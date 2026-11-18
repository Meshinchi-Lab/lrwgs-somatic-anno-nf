// SVINFT_ANNOTATE -- same-type reference DEL/DUP matches for the pass-2 SV tables, in the
// AnnotSV "SVincludedInFt" style with a coverage floor.
//
// A reference BED given to AnnotSV as FtIncludedInSV attaches an entry only when it lies
// ENTIRELY INSIDE the called SV, so a precise long-read call that sits slightly inside its
// published counterpart is missed. This step attaches, to every DEL or DUP record, each entry
// of the same type in the reference DEL/DUP BEDs that the SV covers by at least
// `min_ref_coverage` of the entry's length; the floor stops megabase-scale entries from
// matching every small SV inside them. The BEDs use the AnnotSV user-BED format (the T-ALL
// defaults are the ST17 DEL/DUP sets).
//
// Output: the pass-2 table with every original column written back byte for byte, plus
// SVinFt_DEL and SVinFt_DUP ("name [score]" entries joined by "|", names encoded as AnnotSV
// writes them). SV_REPORT_INDEX merges those two columns into AnnotSV's matches for the same
// BEDs before its cleaning step. Same container and R stack as SV_REPORT_INDEX; the EnsDb
// supplies chromosome names and lengths. The workflow skips this step when either BED or the
// EnsDb is unset or not found.
process SVINFT_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/meshinchi-lab/sv_report_index:2026.11' :
        'ghcr.io/meshinchi-lab/sv_report_index:2026.11' }"

    input:
    tuple val(meta), path(tsv)
    path del_bed
    path dup_bed
    path ensembl_db
    val  min_ref_coverage

    output:
    // A NEW file name (never the input's): the input is a symlink into FILTER_SV's output,
    // and writing to the same name would overwrite it. The `.filtered.pass2.tsv` suffix is
    // what the report's `\.filtered\.pass2\.tsv$` glob matches.
    tuple val(meta), path("${prefix}.filtered.pass2.tsv"), emit: tsv
    tuple val("${task.process}"), val('r-base'),
          eval("Rscript -e 'cat(R.version[[\"major\"]], \".\", R.version[[\"minor\"]], sep=\"\")'"),
          topic: versions, emit: versions_r
    tuple val("${task.process}"), val('plyranges'),
          eval("Rscript -e 'cat(as.character(packageVersion(\"plyranges\")))'"),
          topic: versions, emit: versions_plyranges

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.svinft"
    """
    svinft_annotate.R ${tsv} ${prefix}.filtered.pass2.tsv ${del_bed} ${dup_bed} ${ensembl_db} ${min_ref_coverage}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.svinft"
    """
    awk 'BEGIN { OFS = "\\t" } NR == 1 { print \$0, "SVinFt_DEL", "SVinFt_DUP"; next } { print \$0, "", "" }' ${tsv} > ${prefix}.filtered.pass2.tsv
    """
}
