import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from exomiser_vcf_path import is_gzip_file  # noqa: E402


def test_gzip_detected_by_content_not_extension(tmp_path):
    p = tmp_path / "input_for_exomiser"  # no extension, like the real staging file
    p.write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n"))
    assert is_gzip_file(str(p))


def test_plain_text_is_not_gzip(tmp_path):
    p = tmp_path / "x.vcf.gz"  # misleading extension
    p.write_text("##fileformat=VCFv4.2\n")
    assert not is_gzip_file(str(p))


def test_missing_file_is_not_gzip(tmp_path):
    assert not is_gzip_file(str(tmp_path / "nope"))
