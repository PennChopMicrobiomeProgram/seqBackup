from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import stat
import tempfile
import warnings
from pathlib import Path

from seqBackupLib.nanopore import NanoporeDir, NanoporeFastq, natural_sort_key

# A real ONT run is always well over this; the check exists to catch a run that
# was pointed at the wrong folder or copied while still in progress.
DEFAULT_MIN_TOTAL_SIZE = 500000000  # 500MB

FASTQ_PASS_DIRNAME = "fastq_pass"

# Where runs are almost always archived.
DEFAULT_DESTINATION_DIR = Path("/mnt/isilon/microbiome/raw_data")

# Top-level run files worth keeping alongside the reads.  ``*.txt`` is
# intentionally narrow: we want ``final_summary_*.txt`` but not the (huge)
# ``sequencing_summary_*.txt``.
REPORT_GLOBS = (
    "*.csv",
    "*.json",
    "*.tsv",
    "*.md",
    "report_*.html",
    "final_summary_*.txt",
)

READ_ONLY = stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH


def return_md5(fp: Path) -> str:
    # from https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    with open(fp, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def find_fastq_groups(fastq_pass_dir: Path) -> dict[str, list[Path]]:
    """Map an output name to the ordered list of fastq.gz chunks to concatenate.

    A multiplexed run has ``fastq_pass/barcode01/``, ``fastq_pass/unclassified/``
    etc., each holding many chunk files -> one output file per subdirectory.  A
    non-multiplexed run drops the chunks straight into ``fastq_pass/`` -> a
    single ``fastq_pass.fastq.gz`` output.
    """
    groups: dict[str, list[Path]] = {}

    for subdir in sorted(d for d in fastq_pass_dir.iterdir() if d.is_dir()):
        chunks = sorted(subdir.glob("*.fastq.gz"), key=natural_sort_key)
        if chunks:
            groups[subdir.name] = chunks

    top_level = sorted(fastq_pass_dir.glob("*.fastq.gz"), key=natural_sort_key)
    if top_level:
        groups.setdefault(FASTQ_PASS_DIRNAME, top_level)

    return groups


def concatenate_gzips(chunks: list[Path], dest: Path) -> str:
    """Concatenate gzip members byte-for-byte (a valid multi-member gzip).

    Returns the md5 of the written file, computed in the same pass so a large
    archive is not read back a second time.
    """
    hash_md5 = hashlib.md5()
    with open(dest, "wb") as out:
        for chunk in chunks:
            with open(chunk, "rb") as handle:
                while data := handle.read(1024 * 1024):
                    hash_md5.update(data)
                    out.write(data)
    return hash_md5.hexdigest()


def _expected_barcode(stem: str) -> str | None:
    if stem == "unclassified" or stem.startswith("barcode"):
        return stem
    return None


def backup_nanopore(
    run_dir: Path,
    dest_dir: Path,
    min_total_size: int = DEFAULT_MIN_TOTAL_SIZE,
    allow_check_failures: bool = False,
) -> Path:
    run_dir = Path(run_dir)
    dest_dir = Path(dest_dir)

    # Parsing the folder name determines the archive location, so this always
    # hard-fails regardless of --allow-check-failures.
    nd = NanoporeDir(run_dir.name)

    fastq_pass = run_dir / FASTQ_PASS_DIRNAME
    if not fastq_pass.is_dir():
        raise IOError("fastq_pass directory not found", str(fastq_pass))

    groups = find_fastq_groups(fastq_pass)
    if not groups:
        raise IOError("No .fastq.gz files found under fastq_pass", str(fastq_pass))

    write_dir = dest_dir / nd.build_archive_dir()
    if write_dir.exists():
        raise FileExistsError(f"Archive directory already exists: {write_dir}")

    dest_dir.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(dir=dest_dir, prefix=f".{nd.build_archive_dir()}.staging.")
    )
    try:
        md5s: list[tuple[str, str]] = []
        header_failures: list[tuple[str, dict]] = []
        total_size = 0

        for stem, chunks in groups.items():
            out_fp = staging / f"{stem}.fastq.gz"
            digest = concatenate_gzips(chunks, out_fp)
            total_size += out_fp.stat().st_size
            md5s.append((out_fp.name, digest))

            with gzip.open(out_fp, "rt") as handle:
                nf = NanoporeFastq(
                    handle,
                    expected_flowcell=nd.folder_info["flowcell_id"],
                    expected_barcode=_expected_barcode(stem),
                )
            ok, problems = nf.check_fp_vs_content()
            if not ok:
                header_failures.append((out_fp.name, problems))

        if total_size < min_total_size:
            message = (
                f"Concatenated fastq files total {total_size} bytes, below the minimum "
                f"of {min_total_size}. Check the run folder or lower --min-total-size."
            )
            if allow_check_failures:
                warnings.warn(message)
            else:
                raise ValueError(message)

        if header_failures:
            message = (
                "FASTQ header info does not match the run folder",
                header_failures,
            )
            if allow_check_failures:
                warnings.warn(f"{message[0]}: {message[1]}")
            else:
                raise ValueError(*message)

        for pattern in REPORT_GLOBS:
            for src in sorted(run_dir.glob(pattern)):
                if src.is_file():
                    shutil.copyfile(src, staging / src.name)

        md5_out_fp = staging / ".".join([nd.build_archive_dir(), "md5"])
        with open(md5_out_fp, "w") as md5_out:
            for name, digest in md5s:
                md5_out.write("\t".join([name, digest]) + "\n")

        for child in staging.iterdir():
            child.chmod(READ_ONLY)

        # All checks passed: publish the staged archive.
        write_dir.mkdir(parents=True, exist_ok=False)
        for child in staging.iterdir():
            shutil.move(str(child), str(write_dir / child.name))
    except BaseException:
        if write_dir.is_dir():
            shutil.rmtree(write_dir, ignore_errors=True)
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    return write_dir


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Concatenate and back up Oxford Nanopore fastq_pass reads"
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="The MinKNOW run folder (contains fastq_pass/ and the run reports).",
    )
    parser.add_argument(
        "--destination-dir",
        required=False,
        type=Path,
        default=DEFAULT_DESTINATION_DIR,
        help=f"Destination folder to write the archive to (default: {DEFAULT_DESTINATION_DIR}).",
    )
    parser.add_argument(
        "--min-total-size",
        required=False,
        type=int,
        default=DEFAULT_MIN_TOTAL_SIZE,
        help="Minimum combined size (bytes) of all concatenated fastq files.",
    )
    parser.add_argument(
        "--allow-check-failures",
        action="store_true",
        help="Continue archiving even if validation checks fail",
    )
    args = parser.parse_args(argv)
    return backup_nanopore(
        args.run_dir,
        args.destination_dir,
        args.min_total_size,
        args.allow_check_failures,
    )
