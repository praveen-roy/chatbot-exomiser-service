"""VCF path helpers for the Exomiser worker.

Reconstructed from the service's usage contract (exomiser_service_app.py): the VCF is
streamed from S3 to a staging file with no real extension
(`input_for_exomiser.download`), then normalized here to a path whose extension matches
its actual content, because Exomiser's CLI picks its VCF reader from the file extension
(`.vcf` vs `.vcf.gz`). A mislabeled extension makes Exomiser fail to parse the VCF.

Two functions, matching the imports in exomiser_service_app.py:
  - is_gzip_file(path) -> bool
  - normalize_exomiser_vcf_path(path, tmpdir) -> str  (returns the corrected path)
"""
from __future__ import annotations

import os
import shutil

# gzip magic number: first two bytes 0x1f 0x8b
_GZIP_MAGIC = b"\x1f\x8b"


def is_gzip_file(path: str) -> bool:
    """True if the file begins with the gzip magic bytes (content-based, not extension)."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == _GZIP_MAGIC
    except (OSError, IOError):
        return False


def normalize_exomiser_vcf_path(path: str, tmpdir: str) -> str:
    """Give the VCF a correct extension for Exomiser based on its ACTUAL content.

    The downloaded staging file has no meaningful extension. Exomiser selects its reader
    from the extension, so:
      - gzip content  -> ensure the path ends in `.vcf.gz`
      - plain content -> ensure the path ends in `.vcf`
    Returns the (possibly renamed) path. Idempotent if already correctly named.
    """
    gz = is_gzip_file(path)
    target_ext = ".vcf.gz" if gz else ".vcf"

    lower = path.lower()
    already_correct = (gz and lower.endswith(".vcf.gz")) or (
        not gz and lower.endswith(".vcf") and not lower.endswith(".vcf.gz")
    )
    if already_correct:
        return path

    target = os.path.join(tmpdir, "input_for_exomiser" + target_ext)
    if os.path.abspath(target) == os.path.abspath(path):
        return path
    # Move (not copy) to avoid doubling disk use on large VCFs; fall back to copy across fs.
    try:
        os.replace(path, target)
    except OSError:
        shutil.copyfile(path, target)
    return target
