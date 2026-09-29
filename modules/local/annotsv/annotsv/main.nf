process ANNOTSV_ANNOTSV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/363f212881f1b2f5c3395a6c7d1270694392e3a6f886e46e091e83527fed9b6b/data' :
        'community.wave.seqera.io/library/annotsv:3.5.3--71a461cb86d570b7' }"

    // Singularity uses --writable-tmpfs so variantconvert can write its cached
    // combined.local.json config. Docker handles this by copying variantconvert to
    // the writable task work dir in the script block (see below).
    containerOptions "${ workflow.containerEngine in ['singularity', 'apptainer'] ? '--writable-tmpfs' : ''}"

    input:
    tuple val(meta), path(sv_vcf), path(sv_vcf_index), path(candidate_small_variants)
    tuple val(meta2), path(annotations)
    tuple val(meta3), path(candidate_genes)
    tuple val(meta4), path(false_positive_snv)
    tuple val(meta5), path(gene_transcripts)

    output:
    tuple val(meta), path("*.annotated.tsv")    , emit: tsv
    tuple val(meta), path("*.unannotated.tsv")  , emit: unannotated_tsv, optional: true
    tuple val(meta), path("*.reheader.vcf")     , emit: vcf            , optional: true
    tuple val(meta), path("*.vc_config.json")   , emit: vc_config      , optional: true
    tuple val("${task.process}"), val('annotsv'), eval("AnnotSV --version | sed 's/AnnotSV //'"), emit: versions_annotsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def genome_build = params.genome ?: 'GRCh38'

    def cand_genes     = candidate_genes          ? "-candidateGenesFile ${candidate_genes}"              : ""
    def small_variants = candidate_small_variants ? "-candidateSnvIndelFiles ${candidate_small_variants}" : ""
    def fp_snv         = false_positive_snv       ? "-snvIndelFiles ${false_positive_snv}"                : ""
    def transcripts    = gene_transcripts         ? "-txFile ${gene_transcripts}"                         : ""

    // For Docker: copy variantconvert to the task work dir so the combined.local.json
    // cache file can be written. Singularity is handled by --writable-tmpfs above.
    // PYTHONPATH must include src/ so `import variantconvert` resolves when AnnotSV
    def vc_setup = workflow.containerEngine == 'docker'
        ? "cp -r /opt/conda/share/python3/variantconvert \${PWD}/variantconvert_rw\nexport PYTHONPATH=\"\${PWD}/variantconvert_rw/src:\${PYTHONPATH:-}\""
        : ""
    def vc_dir = workflow.containerEngine == 'docker'
        ? "-variantconvertDir \${PWD}/variantconvert_rw"
        : ""

    """
    ${vc_setup}

    # Inject original VCF ##INFO field metadata into all variantconvert base *.json
    # configs before AnnotSV runs. variantconvert reconstructs ##INFO headers from its
    # COLUMNS_DESCRIPTION.INFO dict; fields declared there use the config's Number/Type/
    # Description, while unknown fields fall back to Number=./Type=String and a generic
    # "Imported from the INFO field" description. Pre-populating the config with the
    # original VCF headers ensures CALLER, PRECISE, IMPRECISE, MATEID, etc. are emitted
    # with their correct metadata in the variantconvert output VCF.
    # After AnnotSV runs, the generated .local.json is copied as ${prefix}.vc_config.json
    # for auditability — it records exactly which INFO field declarations were in effect.
    python3 << 'VCCONFIG'
import gzip, json, pathlib, re
import variantconvert

def open_vcf(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, 'r')

# Parse ##INFO header lines from the input SV VCF.
info_fields = {}
with open_vcf('${sv_vcf}') as fh:
    for line in fh:
        if not line.startswith('#'):
            break
        m = re.match(r'##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="(.*)">', line.rstrip())
        if m:
            info_fields[m.group(1)] = {'Number': m.group(2), 'Type': m.group(3), 'Description': m.group(4)}

# Locate variantconvert configs directory relative to the installed package.
vc_configs = pathlib.Path(variantconvert.__file__).parent / 'configs'

# Inject any INFO field not already declared in COLUMNS_DESCRIPTION.INFO across
# all base *.json configs (shipped with container; .local.json files are generated
# by variantconvert at runtime and do not exist before AnnotSV runs).
for cfg_path in vc_configs.rglob('*.json'):
    if '.local.' in cfg_path.name:
        continue
    cfg = json.loads(cfg_path.read_text())
    cd_info = cfg.setdefault('COLUMNS_DESCRIPTION', {}).setdefault('INFO', {})
    changed = False
    for fid, fmeta in info_fields.items():
        if fid not in cd_info:
            cd_info[fid] = fmeta
            changed = True
    if changed:
        cfg_path.write_text(json.dumps(cfg, indent=4))
VCCONFIG

    AnnotSV \\
        -annotationsDir ${annotations} \\
        ${cand_genes} \\
        ${small_variants} \\
        ${fp_snv} \\
        ${transcripts} \\
        -outputFile ${prefix}.tsv \\
        -SVinputFile ${sv_vcf} \\
        ${vc_dir} \\
        ${args}

    mv *_AnnotSV/* .

    # Copy the runtime-generated variantconvert config as an audit artifact.
    # .local.json is created by variantconvert when AnnotSV calls it, so it exists here.
    python3 << 'VCCOPY'
import pathlib, shutil, variantconvert
vc = pathlib.Path(variantconvert.__file__).parent / 'configs' / '${genome_build}' / 'annotsv3_from_vcf.combined.local.json'
if vc.exists():
    shutil.copy(str(vc), '${prefix}.vc_config.json')
VCCOPY

    # AnnotSV emits INFO field names with parentheses, slashes, and quotes that
    # violate the VCF spec (only alphanumeric + underscore allowed). Write a
    # corrected copy as ${prefix}.reheader.vcf; the original variantconvert
    # output (${prefix}.vcf) is kept intact for debugging.
    #
    # ACMG_class is forced to Number=. (not Number=1) because AnnotSV's
    # split-mode rows flatten into a single VCF record with comma-separated
    # ACMG_class values like ".,full=NA,full=NA,...". vembrane strictly enforces
    # Number=1 declarations and rejects multi-value records, so any record with
    # multi-gene annotation overlap would be silently dropped. Declaring
    # Number=. lets vembrane accept the records as InfoTuple; downstream filter
    # expressions in FILTER_SV normalize InfoTuple → scalar via
    # str(x).split(",")[0] before comparison.
    #
    # SVTYPE and AnnotSV_ranking_score are forced to Number=1 because they are
    # genuinely single-valued per variant.
    #
    # User-BED key normalization: when a BED file passed to -FtIncludedInSV /
    # -SVincludedInFt / -AnyOverlap lacks the standard column-header line
    # (#chrom chromStart chromEnd name score strand), AnnotSV emits INFO keys
    # using positional indices: '_<stem>' (col 4), '_<stem>.1' (col 5),
    # '_<stem>.2' (col 6) instead of 'name_<stem>', 'score_<stem>',
    # 'strand_<stem>'. The six rules below rename both header (##INFO=<ID=...>)
    # and data-line (;_stem=val) occurrences to the standard names so the
    # filter / table / qmd / yaml consumers downstream see uniform keys.
    # Order matters: .2 → strand_, then .1 → score_, then bare → name_.
    # No other AnnotSV / Severus / SAVANA / MINDA INFO keys begin with a
    # leading underscore, so the generic '_<stem>' pattern is collision-safe.
    if [ -f "${prefix}.vcf" ]; then
        sed -E \\
            -e "s|'SampleID'|SampleID|g" \\
            -e 's|Compound_htz\\(sample\\)|Compound_htz_sample|g' \\
            -e 's|Count_hom\\(sample\\)|Count_hom_sample|g' \\
            -e 's|Count_htz\\(sample\\)|Count_htz_sample|g' \\
            -e 's|Count_htz/allHom\\(sample\\)|Count_htz_allHom_sample|g' \\
            -e 's|Count_htz/total\\(cohort\\)|Count_htz_total_cohort|g' \\
            -e 's|Count_total\\(cohort\\)|Count_total_cohort|g' \\
            -e 's|<ID=SVTYPE,Number=[^,]*,|<ID=SVTYPE,Number=1,|g' \\
            -e 's|<ID=ACMG_class,Number=[^,]*,|<ID=ACMG_class,Number=.,|g' \\
            -e 's|<ID=AnnotSV_ranking_score,Number=[^,]*,|<ID=AnnotSV_ranking_score,Number=1,|g' \\
            -e 's|<ID=_([^,]+)\\.2,|<ID=strand_\\1,|g' \\
            -e 's|<ID=_([^,]+)\\.1,|<ID=score_\\1,|g' \\
            -e 's|<ID=_([^,]+),|<ID=name_\\1,|g' \\
            -e 's|([[:space:];])_([^;=]+)\\.2=|\\1strand_\\2=|g' \\
            -e 's|([[:space:];])_([^;=]+)\\.1=|\\1score_\\2=|g' \\
            -e 's|([[:space:];])_([^;=]+)=|\\1name_\\2=|g' \\
            "${prefix}.vcf" > "${prefix}.reheader.vcf"
    fi
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def create_vcf = args.contains("-vcf 1") ? "touch ${prefix}.vcf\ntouch ${prefix}.reheader.vcf\ntouch ${prefix}.vc_config.json" : ""

    """
    touch ${prefix}.tsv
    touch ${prefix}.unannotated.tsv
    ${create_vcf}
    """
}
