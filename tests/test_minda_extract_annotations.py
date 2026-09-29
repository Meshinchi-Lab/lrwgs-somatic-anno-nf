"""Tests for bin/minda_extract_annotations.py"""
import sys
import os
import gzip
import subprocess
from pathlib import Path

SAVANA_VCF = """\
##fileformat=VCFv4.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=MATEID,Number=1,Type=String,Description="Mate BND ID">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcolo829
chr21\t10000\tID_1_1\tN\t<DEL>\t.\t.\tSVTYPE=BND;END=20000;MATEID=ID_1_2\t.\t.
chr21\t20000\tID_1_2\tN\t]chr21:10000]\t.\t.\tSVTYPE=BND;END=20001;MATEID=ID_1_1\t.\t.
chr22\t50000\tID_2_1\tN\t<INS>\t.\t.\tSVTYPE=BND;END=50001\t.\t.
"""

SEVERUS_VCF = """\
##fileformat=VCFv4.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=MATE_ID,Number=1,Type=String,Description="Mate BND ID">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcolo829
chr21\t9990\tseverus_BND1_1\tN\t<BND>\t.\t.\tSVTYPE=BND;END=9990;MATE_ID=severus_BND1_2\t.\t.
chr21\t19990\tseverus_BND1_2\tN\t<BND>\t.\t.\tSVTYPE=BND;END=19990;MATE_ID=severus_BND1_1\t.\t.
chr22\t50100\tseverus_INS3\tN\t<INS>\t.\t.\tSVTYPE=INS;END=50101\t.\t.
"""

MINDA_VCF = """\
##fileformat=VCFv4.1
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Supporting callers">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr21\t9990\t.\tN\t<BND>\t.\t.\tSVTYPE=BND;SUPP_VEC=severus_BND1_1,ID_1_1
chr22\t50100\t.\tN\t<INS>\t.\t.\tSVTYPE=INS;SUPP_VEC=severus_INS3,ID_2_1
"""


def _write(path, content):
    Path(path).write_text(content)


def run_script(minda_vcf, severus_vcf, savana_vcf, sample_id="test"):
    script = Path(__file__).parent.parent / "bin" / "minda_extract_annotations.py"
    cmd = [
        sys.executable, str(script),
        "--minda_vcf",   minda_vcf,
        "--severus_vcf", severus_vcf,
        "--savana_vcf",  savana_vcf,
        "--sample_id",   sample_id,
    ]
    return subprocess.run(cmd, capture_output=True, text=True)


class TestMindaExtractAnnotations:

    def test_produces_three_output_files(self, tmp_path):
        _write(tmp_path / "minda.vcf", MINDA_VCF)
        _write(tmp_path / "severus.vcf", SEVERUS_VCF)
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "minda.vcf"),
            str(tmp_path / "severus.vcf"),
            str(tmp_path / "savana.vcf"),
        )
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "test.annotation.vcf.gz").exists()
        assert (tmp_path / "test.annotation.vcf.gz.tbi").exists()
        assert (tmp_path / "test.header.txt").exists()
        assert (tmp_path / "test.savana_supported_ids.txt").exists()

    def test_annotation_vcf_contains_bnd_partner(self, tmp_path):
        """Both BND partners of a SEVERUS BND pair should appear in the annotation VCF."""
        _write(tmp_path / "minda.vcf", MINDA_VCF)
        _write(tmp_path / "severus.vcf", SEVERUS_VCF)
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "minda.vcf"),
            str(tmp_path / "severus.vcf"),
            str(tmp_path / "savana.vcf"),
        )
        assert r.returncode == 0, r.stderr
        content = gzip.open(tmp_path / "test.annotation.vcf.gz", "rt").read()
        assert "severus_BND1_1" in content
        assert "severus_BND1_2" in content

    def test_supported_ids_contains_both_savana_mates(self, tmp_path):
        """supported_ids.txt must include both BND partners of a SAVANA pair."""
        _write(tmp_path / "minda.vcf", MINDA_VCF)
        _write(tmp_path / "severus.vcf", SEVERUS_VCF)
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "minda.vcf"),
            str(tmp_path / "severus.vcf"),
            str(tmp_path / "savana.vcf"),
        )
        assert r.returncode == 0, r.stderr
        ids = (tmp_path / "test.savana_supported_ids.txt").read_text().splitlines()
        assert "ID_1_1" in ids
        assert "ID_1_2" in ids

    def test_header_contains_all_info_fields(self, tmp_path):
        _write(tmp_path / "minda.vcf", MINDA_VCF)
        _write(tmp_path / "severus.vcf", SEVERUS_VCF)
        _write(tmp_path / "savana.vcf", SAVANA_VCF)
        os.chdir(tmp_path)
        r = run_script(
            str(tmp_path / "minda.vcf"),
            str(tmp_path / "severus.vcf"),
            str(tmp_path / "savana.vcf"),
        )
        assert r.returncode == 0, r.stderr
        header = (tmp_path / "test.header.txt").read_text()
        for field in ["CALLER", "SAVANA_ID", "SAVANA_SVTYPE", "SAVANA_POS",
                      "SAVANA_END", "SAVANA_CHR2", "SAVANA_POS2", "SAVANA_MATEID"]:
            assert field in header, f"Missing INFO field: {field}"
