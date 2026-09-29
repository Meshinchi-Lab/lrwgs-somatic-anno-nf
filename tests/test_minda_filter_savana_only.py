"""Tests for bin/minda_filter_savana_only.py"""
import sys
import os
import subprocess
from pathlib import Path

SAVANA_VCF = """\
##fileformat=VCFv4.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate BND ID">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcolo829
chr21\t10000\tID_1_1\tN\t<DEL>\t.\t.\tSVTYPE=BND;END=20000\t.\t.
chr21\t20000\tID_1_2\tN\t]chr21:10000]\t.\t.\tSVTYPE=BND;END=20001\t.\t.
chr22\t50000\tID_2_1\tN\t<INS>\t.\t.\tSVTYPE=BND;END=50001\t.\t.
chr22\t70000\tID_3_1\tN\t<DEL>\t.\t.\tSVTYPE=BND;END=80000\t.\t.
"""

SUPPORTED_IDS = "ID_1_1\nID_1_2\n"

HEADER_TXT = """\
##INFO=<ID=CALLER,Number=1,Type=String,Description="Originating SV caller">
##INFO=<ID=SAVANA_ID,Number=1,Type=String,Description="SAVANA variant ID from MINDA SUPP_VEC">
##INFO=<ID=SAVANA_SVTYPE,Number=1,Type=String,Description="SAVANA SV type">
##INFO=<ID=SAVANA_POS,Number=1,Type=Integer,Description="SAVANA ALT start position">
##INFO=<ID=SAVANA_END,Number=1,Type=Integer,Description="SAVANA END position">
##INFO=<ID=SAVANA_CHR2,Number=1,Type=String,Description="SAVANA second chromosome (BND only)">
##INFO=<ID=SAVANA_POS2,Number=1,Type=Integer,Description="SAVANA second breakpoint position (BND only)">
##INFO=<ID=SAVANA_MATEID,Number=1,Type=String,Description="SAVANA BND mate ID">
"""


def _write(path, content):
    Path(path).write_text(content)


def run_script(savana_vcf, supported_ids, header_txt, sample_id="test"):
    script = Path(__file__).parent.parent / "bin" / "minda_filter_savana_only.py"
    cmd = [
        sys.executable, str(script),
        "--savana_vcf",    savana_vcf,
        "--supported_ids", supported_ids,
        "--header_txt",    header_txt,
        "--sample_id",     sample_id,
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


class TestMindaFilterSavanaOnly:

    def test_produces_output_vcf(self, tmp_path):
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        _write(tmp_path / "supported.txt", SUPPORTED_IDS)
        _write(tmp_path / "header.txt", HEADER_TXT)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "savana.vcf"),
            str(tmp_path / "supported.txt"),
            str(tmp_path / "header.txt"),
        )
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "test.savana_only.vcf").exists()

    def test_supported_ids_are_excluded(self, tmp_path):
        """ID_1_1 and ID_1_2 are in supported_ids — they must not appear in output."""
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        _write(tmp_path / "supported.txt", SUPPORTED_IDS)
        _write(tmp_path / "header.txt", HEADER_TXT)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "savana.vcf"),
            str(tmp_path / "supported.txt"),
            str(tmp_path / "header.txt"),
        )
        assert r.returncode == 0, r.stderr
        content = (tmp_path / "test.savana_only.vcf").read_text()
        data_lines = [l for l in content.splitlines() if not l.startswith('#')]
        ids = [l.split('\t')[2] for l in data_lines]
        assert "ID_1_1" not in ids
        assert "ID_1_2" not in ids

    def test_unsupported_ids_are_included(self, tmp_path):
        """ID_2_1 and ID_3_1 are not in supported_ids — they must appear in output."""
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        _write(tmp_path / "supported.txt", SUPPORTED_IDS)
        _write(tmp_path / "header.txt", HEADER_TXT)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "savana.vcf"),
            str(tmp_path / "supported.txt"),
            str(tmp_path / "header.txt"),
        )
        assert r.returncode == 0, r.stderr
        content = (tmp_path / "test.savana_only.vcf").read_text()
        data_lines = [l for l in content.splitlines() if not l.startswith('#')]
        ids = [l.split('\t')[2] for l in data_lines]
        assert "ID_2_1" in ids
        assert "ID_3_1" in ids

    def test_output_contains_new_info_headers(self, tmp_path):
        """Output VCF header must include the 8 SAVANA_* INFO lines from header.txt."""
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        _write(tmp_path / "supported.txt", SUPPORTED_IDS)
        _write(tmp_path / "header.txt", HEADER_TXT)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "savana.vcf"),
            str(tmp_path / "supported.txt"),
            str(tmp_path / "header.txt"),
        )
        assert r.returncode == 0, r.stderr
        content = (tmp_path / "test.savana_only.vcf").read_text()
        assert "SAVANA_ID" in content
        assert "CALLER" in content

    def test_sample_column_renamed(self, tmp_path):
        """The sample column in #CHROM line must be replaced with sample_id."""
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        _write(tmp_path / "supported.txt", SUPPORTED_IDS)
        _write(tmp_path / "header.txt", HEADER_TXT)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "savana.vcf"),
            str(tmp_path / "supported.txt"),
            str(tmp_path / "header.txt"),
            sample_id="mysample",
        )
        assert r.returncode == 0, r.stderr
        content = (tmp_path / "mysample.savana_only.vcf").read_text()
        chrom_line = next(l for l in content.splitlines() if l.startswith('#CHROM'))
        assert chrom_line.endswith("mysample")
        assert "colo829" not in chrom_line
