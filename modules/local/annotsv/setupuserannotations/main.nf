process ANNOTSV_SETUPUSERANNO {
    tag "${genome_build}"
    label 'process_low'

    conda null
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/363f212881f1b2f5c3395a6c7d1270694392e3a6f886e46e091e83527fed9b6b/data' :
        'community.wave.seqera.io/library/annotsv:3.5.3--71a461cb86d570b7' }"

    input:
    tuple val(meta), path(annotations)
    path ft_included_in_sv
    path sv_included_in_ft
    path any_overlap
    val genome_build
    // Optional one-gene-per-line candidate-gene panel. `[]` sentinel from the
    // workflow when params.candidateGenesFile is not set. Staged under
    // candidate_genes_in/ so the shell's target dir doesn't collide with the
    // input path when we `cp` into ./candidate_genes/. The rest of the module
    // is agnostic to this file — it's only staged so downstream consumers
    // (IGV_REPORTS_SV, SV_REPORT_INDEX, ANNOTATIONS_SV, ANNOTATIONS_CNA) can
    // reach it through a single canonical channel: SETUPUSERANNO.out.genes.
    path candidate_genes, stageAs: 'candidate_genes_in/*'

    output:
    tuple val(meta), path("${prefix}"), emit: annotations
    // Optional emit — absent when `candidate_genes` was `[]`. Downstream
    // consumers use `.ifEmpty([])` to keep their channel non-blocking.
    path "candidate_genes/*", optional: true, emit: genes
    tuple val("${task.process}"), val('annotsv'), eval("AnnotSV --version | sed 's/AnnotSV //'"), emit: versions_annotsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}.annotsv_ref" : "annotsv_user_annotations")
    def args = task.ext.args ?: ''

    """
    # Copy annotations into a real writable directory.
    # ${annotations} is staged as a symlink; appending /. forces cp to copy the
    # directory CONTENTS rather than re-creating the symlink itself (macOS + Linux).
    mkdir -p ${prefix}
    cp -rL ${annotations}/. ${prefix}/ || true

    mkdir -p ${prefix}/Annotations_Human/Users/${genome_build}/FtIncludedInSV
    mkdir -p ${prefix}/Annotations_Human/Users/${genome_build}/SVincludedInFt
    mkdir -p ${prefix}/Annotations_Human/Users/${genome_build}/AnyOverlap

    for f in ${ft_included_in_sv}; do
        [ -f "\$f" ] && cp -fL "\$f" ${prefix}/Annotations_Human/Users/${genome_build}/FtIncludedInSV/ || true
    done
    for f in ${sv_included_in_ft}; do
        [ -f "\$f" ] && cp -fL "\$f" ${prefix}/Annotations_Human/Users/${genome_build}/SVincludedInFt/ || true
    done
    for f in ${any_overlap}; do
        [ -f "\$f" ] && cp -fL "\$f" ${prefix}/Annotations_Human/Users/${genome_build}/AnyOverlap/ || true
    done

    # Inject a #chrom header into any BED file that lacks one; AnnotSV requires it
    # to produce named annotation columns in the VCF output.
    ensure_bed_headers.sh ${prefix} ${genome_build}

    # Remove any pre-existing .formatted.sorted.bed files so checkUsersBED processes
    # the original .bed files and generates fresh .formatted.sorted.header.tsv files. Only delete when a source-side <base>.bed exists 
    find ${prefix}/Annotations_Human/Users/${genome_build} \
        -name "*.formatted.sorted.bed" \
        -exec sh -c 'base="\${1%.formatted.sorted.bed}"  -f "\${base}.bed" ] && rm "\$1" ' _ {} \\; 2>/dev/null || true

    # log that the files are created in the user annotations dir and use same permissions as the other annotations 
    ls -alhR ${prefix}/Annotations_Human/Users/${genome_build}
    chmod -R 775 ${prefix}/Annotations_Human/Users/${genome_build}


    # ANNOTSV_INSTALLANNOTATIONS normally supplies these, but a pre-built bundle passed via params.annotsv_annotations may be missing them
    # if so, download the missing file before the annotsv command to avoid an error.
    #
    # Only GRCh37/GRCh38 use ncbiRefSeq.txt.gz; CHM13 uses hs1_curGene*.txt.gz and
    # is left untouched (its bundle always ships the formatted BED).
    genes_dir="${prefix}/Annotations_Human/Genes/${genome_build}"
    mkdir -p "\$genes_dir"

    case "${genome_build}" in
        GRCh38) ucsc_db="hg38" ;;
        GRCh37) ucsc_db="hg19" ;;
        *)      ucsc_db=""     ;;
    esac

    if [ ! -f "\$genes_dir/ncbiRefSeq.txt.gz" ] && [ ! -f "\$genes_dir/genes.RefSeq.sorted.bed" ]; then
        if [ -z "\$ucsc_db" ]; then
            echo "ERROR: ${genome_build} has no RefSeq gene file and no UCSC fallback is defined." >&2
            exit 1
        fi
        ncbi_url="https://hgdownload.soe.ucsc.edu/goldenPath/\${ucsc_db}/database/ncbiRefSeq.txt.gz"
        echo "INFO: no ncbiRefSeq.txt.gz or genes.RefSeq.sorted.bed in \$genes_dir; downloading \$ncbi_url"
        # Download to a temp name so a truncated transfer never looks complete to AnnotSV.
        if command -v curl >/dev/null 2>&1; then
            curl -fsSL --retry 3 --retry-delay 5 -o "\$genes_dir/ncbiRefSeq.txt.gz.part" "\$ncbi_url"
        elif command -v wget >/dev/null 2>&1; then
            wget -q --tries=3 -O "\$genes_dir/ncbiRefSeq.txt.gz.part" "\$ncbi_url"
        else
            echo "ERROR: neither curl nor wget is available to fetch \$ncbi_url" >&2
            exit 1
        fi
        # Verify it is a readable gzip before committing the final name.
        gzip -t "\$genes_dir/ncbiRefSeq.txt.gz.part"
        mv "\$genes_dir/ncbiRefSeq.txt.gz.part" "\$genes_dir/ncbiRefSeq.txt.gz"
        echo "INFO: staged \$(du -h "\$genes_dir/ncbiRefSeq.txt.gz" | cut -f1) ncbiRefSeq.txt.gz for ${genome_build}"
    else
        echo "INFO: RefSeq gene source already present in \$genes_dir; skipping UCSC download"
    fi

    # Run AnnotSV on a dummy VCF to trigger checkUsersBED, which generates
    # .formatted.sorted.bed and .formatted.sorted.header.tsv for each user BED file.
    cat > dummy_sv.vcf << 'DUMMYVCF'
##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1000000	.	N	<DEL>	.	PASS	SVTYPE=DEL;END=1001000	GT	0/1
DUMMYVCF

    AnnotSV \\
        -annotationsDir ${prefix} \\
        -SVinputFile dummy_sv.vcf \\
        -outputFile dummy_out \\
        -genomeBuild ${genome_build}

    rm -rf dummy_out* dummy_sv.vcf

    # Append each BED file's stem to its AnnotSV-generated .header.tsv column names
    # (so AnnotSV INFO fields are distinguishable per BED source), then remove the
    # original .bed files, keeping only .formatted.sorted.bed and .header.tsv.
    fix_user_beds.sh ${prefix} ${genome_build}

    # Stage the candidate-gene panel (when provided) into ./candidate_genes/
    # so the `genes` output emit catches it and the module's publishDir
    # exposes it alongside the annotation-setup artefacts. Empty when
    # candidate_genes was `[]` — the optional output emit stays absent.
    if [ -n "${candidate_genes}" ] && [ -f "${candidate_genes}" ]; then
        mkdir -p candidate_genes
        cp -fL ${candidate_genes} candidate_genes/
    fi
    """

    stub:
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}.annotsv_ref" : "annotsv_user_annotations")
    def args = task.ext.args ?: ''

    """
    mkdir -p ${prefix}/Annotations_Human/Users/${genome_build}/FtIncludedInSV
    mkdir -p ${prefix}/Annotations_Human/Users/${genome_build}/SVincludedInFt
    mkdir -p ${prefix}/Annotations_Human/Users/${genome_build}/AnyOverlap
    touch ${prefix}/Annotations_Human/Users/${genome_build}/FtIncludedInSV/stub.formatted.sorted.bed
    touch ${prefix}/Annotations_Human/Users/${genome_build}/FtIncludedInSV/stub.header.tsv
    touch ${prefix}/Annotations_Human/Users/${genome_build}/SVincludedInFt/stub.formatted.sorted.bed
    touch ${prefix}/Annotations_Human/Users/${genome_build}/SVincludedInFt/stub.header.tsv
    touch ${prefix}/Annotations_Human/Users/${genome_build}/AnyOverlap/stub.formatted.sorted.bed
    touch ${prefix}/Annotations_Human/Users/${genome_build}/AnyOverlap/stub.header.tsv

    # Mirror the Genes/<build> layout the real script guarantees, so stub runs
    # exercise the same directory shape without touching the network.
    mkdir -p ${prefix}/Annotations_Human/Genes/${genome_build}
    touch ${prefix}/Annotations_Human/Genes/${genome_build}/genes.RefSeq.sorted.bed

    # Stub-mirror the candidate-gene emit so downstream stub tests still see
    # the file when it was passed in.
    if [ -n "${candidate_genes}" ] && [ -f "${candidate_genes}" ]; then
        mkdir -p candidate_genes
        cp -fL ${candidate_genes} candidate_genes/
    fi
    """
}
