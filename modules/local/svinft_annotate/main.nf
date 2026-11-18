// ST17_SVINFT -- same-type ST17 DEL/DUP matches for the pass-2 SV tables, in the AnnotSV
// "SVincludedInFt" style with a coverage floor.
//
// AnnotSV receives the ST17 DEL/DUP BEDs as FtIncludedInSV, so a reference entry is attached
// only when it lies ENTIRELY INSIDE the called SV. A precise long-read call that sits slightly
// inside its published counterpart is missed (e.g. the BCL11B-enhancer tandem duplication
// called as 3H_PBM68877 severus_DUP2115). This step attaches, to every DEL or DUP record, each
// ST17 entry of the same type that the SV covers by at least `min_ref_coverage` of the entry's
// length; the floor stops megabase-scale entries from matching every small SV inside them.
//
// Output: the pass-2 table with every original column written back byte for byte, plus
// ST17_DEL_SVinFt and ST17_DUP_SVinFt ("name [score]" entries joined by "|", names encoded as
// AnnotSV writes them). SV_REPORT_INDEX merges those two columns into name_ST17_DEL/DUP and
// score_ST17_DEL/DUP before its cleaning step. Same container and R stack as SV_REPORT_INDEX;
// the EnsDb supplies chromosome names and lengths. The workflow skips this step when either
// BED or the EnsDb is unset or not found.
process ST17_SVINFT {
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
    prefix = task.ext.prefix ?: "${meta.id}.st17_svinft"
    """
    st17_svinft.R ${tsv} ${prefix}.filtered.pass2.tsv ${del_bed} ${dup_bed} ${ensembl_db} ${min_ref_coverage}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.st17_svinft"
    """
    awk 'BEGIN { OFS = "\\t" } NR == 1 { print \$0, "ST17_DEL_SVinFt", "ST17_DUP_SVinFt"; next } { print \$0, "", "" }' ${tsv} > ${prefix}.filtered.pass2.tsv
    """
}
