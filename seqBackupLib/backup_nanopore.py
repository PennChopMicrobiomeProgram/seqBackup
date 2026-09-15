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
DEFAULT_MIN_FILE_SIZE = 500000000  # 500MB

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


def find_fastq_chunks(fastq_pass_dir: Path) -> list[Path]:
    """All fastq.gz chunks under fastq_pass/, in a stable concatenation order.

    A multiplexed run keeps its chunks under ``fastq_pass/<barcode>/``; a
    non-multiplexed run drops them straight into ``fastq_pass/``. Either way,
    every chunk is merged into one archived file (matching the Illumina
    convention of archiving the undemultiplexed reads and splitting later), so
    sorting is by (containing folder, natural chunk order) purely to make the
    result reproducible and easy to reason about -- not to group output files.
    """
    return sorted(
        fastq_pass_dir.rglob("*.fastq.gz"),
        key=lambda p: (p.parent.name, natural_sort_key(p)),
    )


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


def backup_nanopore(
    run_dir: Path,
    dest_dir: Path,
    sample_sheet: Path,
    min_file_size: int = DEFAULT_MIN_FILE_SIZE,
    allow_check_failures: bool = False,
) -> Path:
    run_dir = Path(run_dir)
    dest_dir = Path(dest_dir)
    sample_sheet = Path(sample_sheet)

    # Parsing the folder name determines the archive location, so this always
    # hard-fails regardless of --allow-check-failures.
    nd = NanoporeDir(run_dir.name)

    # Required so the run's metadata is always recorded; a transfer run with no
    # metadata is expected to point this at a dummy placeholder file.
    if not sample_sheet.is_file():
        raise IOError("Sample sheet does not exist", str(sample_sheet))

    fastq_pass = run_dir / FASTQ_PASS_DIRNAME
    if not fastq_pass.is_dir():
        raise IOError("fastq_pass directory not found", str(fastq_pass))

    chunks = find_fastq_chunks(fastq_pass)
    if not chunks:
        raise IOError("No .fastq.gz files found under fastq_pass", str(fastq_pass))

    write_dir = dest_dir / nd.build_archive_dir()
    if write_dir.exists():
        raise FileExistsError(f"Archive directory already exists: {write_dir}")

    dest_dir.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(dir=dest_dir, prefix=f".{nd.build_archive_dir()}.staging.")
    )
    try:
        flowcell_id = nd.folder_info["flowcell_id"]
        out_fp = staging / f"{flowcell_id}.fastq.gz"

        digest = concatenate_gzips(chunks, out_fp)
        file_size = out_fp.stat().st_size

        with gzip.open(out_fp, "rt") as handle:
            nf = NanoporeFastq(handle, expected_flowcell=flowcell_id)
        ok, problems = nf.check_fp_vs_content()

        if file_size < min_file_size:
            message = (
                f"Concatenated fastq is {file_size} bytes, below the minimum of "
                f"{min_file_size}. Check the run folder or lower --min-file-size."
            )
            if allow_check_failures:
                warnings.warn(message)
            else:
                raise ValueError(message)

        if not ok:
            message = ("FASTQ header info does not match the run folder", problems)
            if allow_check_failures:
                warnings.warn(f"{message[0]}: {message[1]}")
            else:
                raise ValueError(*message)

        for pattern in REPORT_GLOBS:
            for src in sorted(run_dir.glob(pattern)):
                if src.is_file():
                    shutil.copyfile(src, staging / src.name)

        # Copied last so the supplied sheet wins over the run's own
        # sample_sheet_*.csv if they happen to share a name.
        shutil.copyfile(sample_sheet, staging / sample_sheet.name)

        md5_out_fp = staging / ".".join([nd.build_archive_dir(), "md5"])
        with open(md5_out_fp, "w") as md5_out:
            md5_out.write("\t".join([out_fp.name, digest]) + "\n")

        for path in staging.iterdir():
            if path.is_file():
                path.chmod(READ_ONLY)

        # All checks passed: publish the staged archive.
        write_dir.mkdir(parents=True, exist_ok=False)
        for child in sorted(staging.iterdir()):
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
        "--sample-sheet",
        required=True,
        type=Path,
        help="The sample sheet associated with the run.",
    )
    parser.add_argument(
        "--min-file-size",
        required=False,
        type=int,
        default=DEFAULT_MIN_FILE_SIZE,
        help="Minimum size (bytes) of the concatenated fastq file.",
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
        args.sample_sheet,
        args.min_file_size,
        args.allow_check_failures,
    )
