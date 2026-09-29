"""Unit tests for bin/vcf_to_bedpe.py.

Tests use stdlib-only imports (pytest handles the discovery).
The bin/ dir is added to sys.path so the script can be imported as a module.
"""

import subprocess
import sys
from pathlib import Path

BIN_DIR = Path(__file__).parent.parent.parent / 'bin'
sys.path.insert(0, str(BIN_DIR))


def test_import_placeholder():
    """Sanity check: pytest discovers this file and can run."""
    assert BIN_DIR.exists()


def test_parse_bnd_partner_t_bracket_p_bracket():
    """t[p[ form — REF followed by [partner:pos["""
    from vcf_to_bedpe import parse_bnd_partner
    assert parse_bnd_partner('N[chr21:12345[') == ('chr21', 12345)


def test_parse_bnd_partner_bracket_p_bracket_t():
    """[p[t form — [partner:pos[ followed by REF"""
    from vcf_to_bedpe import parse_bnd_partner
    assert parse_bnd_partner('[chr21:12345[N') == ('chr21', 12345)


def test_parse_bnd_partner_t_squarebracket_p_squarebracket():
    """t]p] form — REF followed by ]partner:pos]"""
    from vcf_to_bedpe import parse_bnd_partner
    assert parse_bnd_partner('N]chr21:12345]') == ('chr21', 12345)


def test_parse_bnd_partner_squarebracket_p_squarebracket_t():
    """]p]t form — ]partner:pos] followed by REF"""
    from vcf_to_bedpe import parse_bnd_partner
    assert parse_bnd_partner(']chr21:12345]N') == ('chr21', 12345)


def test_parse_bnd_partner_returns_none_for_symbolic_alt():
    """Symbolic ALT like <DEL> should not match the BND regex."""
    from vcf_to_bedpe import parse_bnd_partner
    assert parse_bnd_partner('<DEL>') == (None, None)


def test_parse_bnd_partner_handles_chr_prefix_and_no_prefix():
    from vcf_to_bedpe import parse_bnd_partner
    assert parse_bnd_partner('N[21:12345[') == ('21', 12345)


def test_format_name_del_with_gene():
    from vcf_to_bedpe import format_name
    assert format_name('sev_DEL42', 'DEL', 100, 'RUNX1') == 'sev_DEL42|DEL|RUNX1'


def test_format_name_del_no_gene_uses_dot_sentinel():
    from vcf_to_bedpe import format_name
    assert format_name('sev_DEL42', 'DEL', 100, '.') == 'sev_DEL42|DEL'


def test_format_name_del_no_gene_empty_string():
    from vcf_to_bedpe import format_name
    assert format_name('sev_DEL42', 'DEL', 100, '') == 'sev_DEL42|DEL'


def test_format_name_ins_encodes_svlen_bp_suffix():
    """INS records substitute verbose ALT with +<SVLEN>bp for readability."""
    from vcf_to_bedpe import format_name
    assert format_name('sev_INS17', 'INS', 250, 'LEF1') == 'sev_INS17|INS|+250bp|LEF1'


def test_format_name_ins_absolute_svlen_for_negative_values():
    """Some callers emit negative SVLEN for INS — use abs()."""
    from vcf_to_bedpe import format_name
    assert format_name('sev_INS17', 'INS', -250, 'LEF1') == 'sev_INS17|INS|+250bp|LEF1'


def test_format_name_bnd():
    from vcf_to_bedpe import format_name
    assert format_name('sev_BND3977', 'BND', 0, 'RUNX1') == 'sev_BND3977|BND|RUNX1'


def test_format_name_bnd_no_gene():
    from vcf_to_bedpe import format_name
    assert format_name('sev_BND3977', 'BND', 0, '.') == 'sev_BND3977|BND'


def test_classify_bucket_snv_always_forces_small():
    """SNV rows land in SMALL regardless of computed distance."""
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 100, 101, 'chr21', 100, 101, is_snv=True) == 'small'


def test_classify_bucket_small_sv_same_chrom_under_1kb():
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 100, 101, 'chr21', 500, 501, is_snv=False) == 'small'


def test_classify_bucket_small_boundary_exactly_999bp():
    """999 bp inter-anchor is SMALL (strict <1000)."""
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 0, 1, 'chr21', 998, 999, is_snv=False) == 'small'


def test_classify_bucket_medium_boundary_exactly_1000bp():
    """1000 bp is MEDIUM (inclusive lower bound)."""
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 0, 1, 'chr21', 999, 1000, is_snv=False) == 'medium'


def test_classify_bucket_medium_sv():
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 100, 101, 'chr21', 5000, 5001, is_snv=False) == 'medium'


def test_classify_bucket_large_boundary_exactly_10000bp():
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 0, 1, 'chr21', 9999, 10000, is_snv=False) == 'large'


def test_classify_bucket_large_sv_over_10kb():
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 100, 101, 'chr21', 50000, 50001, is_snv=False) == 'large'


def test_classify_bucket_bnd_different_chroms_is_large():
    """BND translocation between different chromosomes always LARGE."""
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr21', 100, 101, 'chr5', 12000000, 12000001, is_snv=False) == 'large'


def test_classify_bucket_snv_same_position_is_small_not_large():
    """Regression guard."""
    from vcf_to_bedpe import classify_bucket
    assert classify_bucket('chr1', 999, 1000, 'chr1', 999, 1000, is_snv=True) == 'small'


def test_sv_row_del_emits_span_bedpe():
    """DEL: chrom1 == chrom2, anchors from POS to END (both 0-based half-open)."""
    from vcf_to_bedpe import sv_row_to_bedpe
    row = sv_row_to_bedpe('chr21\t1000\t2000\tDEL\t1000\t<DEL>\tsev_DEL1\tGENE1')
    assert row.chrom1 == 'chr21'
    assert row.start1 == 999   # POS - 1 (0-based BED)
    assert row.end1 == 1000    # POS
    assert row.chrom2 == 'chr21'
    assert row.start2 == 1999
    assert row.end2 == 2000
    assert row.svtype == 'DEL'
    assert row.name == 'sev_DEL1|DEL|GENE1'


def test_sv_row_bnd_uses_partner_from_alt():
    from vcf_to_bedpe import sv_row_to_bedpe
    row = sv_row_to_bedpe(
        'chr21\t34797509\t34797509\tBND\t0\tN]chr5:171274715]\tsev_BND1\tRUNX1'
    )
    assert row.chrom1 == 'chr21'
    assert row.start1 == 34797508
    assert row.end1 == 34797509
    assert row.chrom2 == 'chr5'
    assert row.start2 == 171274714
    assert row.end2 == 171274715
    assert row.svtype == 'BND'
    assert row.name == 'sev_BND1|BND|RUNX1'


def test_sv_row_bnd_unparseable_alt_returns_none():
    """Malformed BND ALT — must not crash, just return None."""
    from vcf_to_bedpe import sv_row_to_bedpe
    row = sv_row_to_bedpe('chr21\t100\t100\tBND\t0\t<DEL>\tsev_BAD\t.')
    assert row is None


def test_sv_row_ins_name_uses_svlen():
    from vcf_to_bedpe import sv_row_to_bedpe
    row = sv_row_to_bedpe('chr21\t500\t500\tINS\t42\tACGT_50bp_seq\tsev_INS1\tLEF1')
    assert row.svtype == 'INS'
    assert row.name == 'sev_INS1|INS|+42bp|LEF1'


def test_sv_row_missing_end_defaults_to_pos():
    """END='.' or empty should default to POS."""
    from vcf_to_bedpe import sv_row_to_bedpe
    row = sv_row_to_bedpe('chr21\t500\t.\tINS\t42\tACGT\tsev_INS2\t.')
    assert row.start1 == 499
    assert row.end1 == 500
    assert row.start2 == 499
    assert row.end2 == 500


def test_sv_row_dot_gene_omitted_from_name():
    from vcf_to_bedpe import sv_row_to_bedpe
    row = sv_row_to_bedpe('chr21\t1000\t2000\tDEL\t1000\t<DEL>\tsev_X\t.')
    assert row.name == 'sev_X|DEL'


def test_snv_row_emits_same_anchor_bedpe():
    from vcf_to_bedpe import snv_row_to_bedpe
    row = snv_row_to_bedpe('chr21\t6116503\trs123\tA|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript')
    assert row.chrom1 == 'chr21'
    assert row.chrom2 == 'chr21'
    assert row.start1 == row.start2 == 6116502
    assert row.end1 == row.end2 == 6116503
    assert row.svtype == 'SNV'
    assert row.name == 'rs123|SNV|TP53'


def test_snv_row_no_csq_uses_dot_gene():
    from vcf_to_bedpe import snv_row_to_bedpe
    row = snv_row_to_bedpe('chr21\t6116503\trs123\t.')
    assert row.name == 'rs123|SNV'


def test_snv_row_empty_csq_uses_dot_gene():
    from vcf_to_bedpe import snv_row_to_bedpe
    row = snv_row_to_bedpe('chr21\t6116503\trs123\t')
    assert row.name == 'rs123|SNV'


def test_snv_row_csq_with_missing_symbol_field():
    """VEP CSQ where SYMBOL field is empty pipe-delimited (Consequence||)."""
    from vcf_to_bedpe import snv_row_to_bedpe
    row = snv_row_to_bedpe('chr21\t6116503\trs123\tA|intergenic_variant|MODIFIER||||')
    assert row.name == 'rs123|SNV'


def test_snv_row_csq_multiple_transcripts_takes_first():
    """VEP emits comma-separated CSQ per transcript; use the first (--pick default)."""
    from vcf_to_bedpe import snv_row_to_bedpe
    row = snv_row_to_bedpe(
        'chr21\t6116503\trs123\t'
        'A|missense_variant|MODERATE|TP53|ENSG0|Transcript,'
        'A|missense_variant|MODERATE|SOMETHING_ELSE|ENSG1|Transcript'
    )
    assert 'TP53' in row.name
    assert 'SOMETHING_ELSE' not in row.name


def test_cli_end_to_end_writes_three_bucket_files(tmp_path):
    """Run vcf_to_bedpe.py as a subprocess with mixed input and verify:
    - all three output files created
    - rows routed to correct buckets
    - SNV lands in small
    - BND lands in large (different chroms)
    """
    sv_tsv = tmp_path / 'sv.tsv'
    sv_tsv.write_text(
        'chr21\t1000\t1500\tDEL\t500\t<DEL>\tsv_small\tGENE_A\n'
        'chr21\t10000\t15000\tDUP\t5000\t<DUP>\tsv_med\tGENE_B\n'
        'chr21\t100000\t150000\tINV\t50000\t<INV>\tsv_large\tGENE_C\n'
        'chr21\t200000\t200000\tBND\t0\tN[chr5:5000000[\tsv_bnd\tGENE_D\n'
    )
    snv_tsv = tmp_path / 'snv.tsv'
    snv_tsv.write_text(
        'chr21\t6116503\trs_snv1\tA|missense_variant|MODERATE|TP53|ENSG0|Transcript\n'
    )
    out_small = tmp_path / 'out.small.bedpe'
    out_med = tmp_path / 'out.medium.bedpe'
    out_large = tmp_path / 'out.large.bedpe'
    script = BIN_DIR / 'vcf_to_bedpe.py'
    result = subprocess.run(
        [
            sys.executable, str(script),
            '--sv-in', str(sv_tsv),
            '--snv-in', str(snv_tsv),
            '--out-small', str(out_small),
            '--out-medium', str(out_med),
            '--out-large', str(out_large),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out_small.exists() and out_med.exists() and out_large.exists()
    small_lines = out_small.read_text().splitlines()
    med_lines = out_med.read_text().splitlines()
    large_lines = out_large.read_text().splitlines()
    assert len(small_lines) == 2  # DEL + SNV
    assert any('sv_small' in ln for ln in small_lines)
    assert any('rs_snv1' in ln for ln in small_lines)
    assert len(med_lines) == 1  # DUP
    assert 'sv_med' in med_lines[0]
    assert len(large_lines) == 2  # INV + BND
    assert any('sv_large' in ln for ln in large_lines)
    assert any('sv_bnd' in ln for ln in large_lines)
    bnd_line = next(ln for ln in large_lines if 'sv_bnd' in ln)
    bnd_fields = bnd_line.split('\t')
    assert bnd_fields[0] == 'chr21'
    assert bnd_fields[3] == 'chr5'


def test_cli_snv_optional_when_no_snv_in_flag(tmp_path):
    """--snv-in is optional; small bucket only has SVs when omitted."""
    sv_tsv = tmp_path / 'sv.tsv'
    sv_tsv.write_text('chr21\t1000\t1500\tDEL\t500\t<DEL>\tsv_small\tGENE_A\n')
    out_small = tmp_path / 'out.small.bedpe'
    out_med = tmp_path / 'out.medium.bedpe'
    out_large = tmp_path / 'out.large.bedpe'
    script = BIN_DIR / 'vcf_to_bedpe.py'
    result = subprocess.run(
        [
            sys.executable, str(script),
            '--sv-in', str(sv_tsv),
            '--out-small', str(out_small),
            '--out-medium', str(out_med),
            '--out-large', str(out_large),
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out_small.read_text().strip() == '\t'.join([
        'chr21', '999', '1000', 'chr21', '1499', '1500',
        'sv_small|DEL|GENE_A', '.', '+', '+', 'DEL',
    ])
    assert out_med.read_text().strip() == ''
    assert out_large.read_text().strip() == ''
