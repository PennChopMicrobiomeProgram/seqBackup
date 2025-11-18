import argparse
import os
import shutil
import stat
import re
import gzip
import hashlib
import warnings
from pathlib import Path
from seqBackupLib.illumina import IlluminaFastq


DEFAULT_MIN_FILE_SIZE = 500000000  # 500MB


def build_fp_to_archive(fp: Path, has_index: bool, lane: str) -> list[Path]:

    if re.search("R1_001.fastq", fp.name) is None:
        raise IOError("The file doesn't look like an R1 file: {}".format(fp))

    label = ["R2"]
    if has_index:
        label.extend(["I1", "I2"])

    if "_L" in fp.name:
        rexp = "".join(["(L00", lane, "_)(R1)(_001.fastq.gz)$"])
        modified_fp = [
            re.sub(rexp, "".join(["\\1", lab, "\\3"]), fp.name) for lab in label
        ]
    else:
        modified_fp = [fp.name.replace("R1", lab) for lab in label]
    return [fp] + [fp.parent / n for n in modified_fp]


def return_md5(fp: Path) -> str:
    # from https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    with open(fp, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def backup_fastq(
    forward_reads: Path,
    dest_dir: Path,
    sample_sheet_fp: Path,
    has_index: bool,
    min_file_size: int,
):

    R1 = IlluminaFastq(gzip.open(forward_reads, mode="rt"))

    # build the strings for the required files
    RI_fps = build_fp_to_archive(forward_reads, has_index, R1.lane)

    # create the Illumina objects and check the files
    illumina_fastqs = [IlluminaFastq(gzip.open(fp, mode="rt")) for fp in RI_fps]
    r1 = illumina_fastqs[0]

    if not all([ifq.check_fp_vs_content()[0] for ifq in illumina_fastqs]):
        [ifq.check_fp_vs_content(verbose=True) for ifq in illumina_fastqs]
        raise ValueError(
            "The file path and header information don't match",
            [str(ifq) for ifq in illumina_fastqs if not ifq.check_fp_vs_content()[0]],
        )
    if not all([ifq.check_file_size(min_file_size) for ifq in illumina_fastqs]):
        raise ValueError(
            "File seems suspiciously small. Please check if you have the correct file or lower the minimum file size threshold",
            [ifq.check_file_size(min_file_size) for ifq in illumina_fastqs],
        )
    if not all([ifq.check_index_read_exists() for ifq in illumina_fastqs]):
        warnings.warn(
            "No barcodes in headers. Were the fastq files generated properly?"
        )

    # parse the info from the headers in EACH file and check they are consistent within each other
    if not all([fastq.is_same_run(illumina_fastqs[0]) for fastq in illumina_fastqs]):
        raise ValueError("The files are not from the same run.")

    ## Archiving steps

    # make sure the sample sheet exists
    if not sample_sheet_fp.is_file():
        raise IOError("Sample sheet does not exist", str(sample_sheet_fp))

    # create the folder to write to
    write_dir = dest_dir / r1.build_archive_dir()
    write_dir.mkdir(parents=True, exist_ok=False)

    ### All the checks are done and the files are safe to archive!

    # move the files to the archive location and set readable permissions
    # keep the files writable by the owner to allow intentional updates or tests
    # that simulate corruption
    permission = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH
    md5s = []
    for fp in RI_fps:
        if "_L" in fp.name:
            dest_name = fp.name
        else:
            dest_name = fp.name.replace("_S0_", f"_S0_L{r1.lane.zfill(3)}_")
        output_fp = write_dir / dest_name
        shutil.copyfile(fp, output_fp)
        output_fp.chmod(permission)
        md5s.append((dest_name, return_md5(fp)))

    # copy the sample sheet to destination folder
    shutil.copyfile(sample_sheet_fp, write_dir / sample_sheet_fp.name)

    # write md5sums to a file
    md5_out_fp = write_dir / ".".join([r1.build_archive_dir(), "md5"])
    with open(md5_out_fp, "w") as md5_out:
        [md5_out.write("\t".join(md5) + "\n") for md5 in md5s]

    return write_dir


def verify_archive(archive_dir: Path) -> bool:
    md5_files = list(archive_dir.glob("*.md5"))
    if not md5_files:
        raise FileNotFoundError(f"No md5 file found in {archive_dir}")

    if len(md5_files) > 1:
        warnings.warn(
            f"Multiple md5 files found in {archive_dir}. Using {md5_files[0].name}."
        )

    md5_fp = md5_files[0]
    missing_files = []
    mismatched_hashes = []

    with open(md5_fp) as md5_file:
        for line in md5_file:
            expected = line.strip().split("\t")
            if len(expected) != 2:
                raise ValueError(f"Invalid md5 line in {md5_fp}: {line}")
            filename, expected_md5 = expected
            file_fp = archive_dir / filename

            if not file_fp.is_file():
                missing_files.append(filename)
                continue

            computed_md5 = return_md5(file_fp)
            if computed_md5 != expected_md5:
                mismatched_hashes.append((filename, expected_md5, computed_md5))

    if missing_files or mismatched_hashes:
        raise ValueError(
            "MD5 verification failed",
            {
                "missing_files": missing_files,
                "mismatched_hashes": mismatched_hashes,
            },
        )

    return True


def main(argv=None):
    parser = argparse.ArgumentParser(description="Backs up fastq files")

    parser.add_argument(
        "--forward-reads", required=True, type=Path, help="Gzipped R1 fastq file"
    )
    parser.add_argument(
        "--destination-dir",
        required=True,
        type=Path,
        help="Destination folder to copy the files to.",
    )
    parser.add_argument(
        "--sample-sheet",
        required=True,
        type=Path,
        help="The sample sheet associated with the run.",
    )
    parser.add_argument(
        "--no-index",
        action="store_true",
        help="Skip index reads (I1/I2) during backup",
    )
    parser.add_argument(
        "--min-file-size",
        required=False,
        type=int,
        default=DEFAULT_MIN_FILE_SIZE,
        help="Minimum file size to register in bytes",
    )
    args = parser.parse_args(argv)
    return backup_fastq(
        args.forward_reads,
        args.destination_dir,
        args.sample_sheet,
        not args.no_index,
        args.min_file_size,
    )

    # maybe also ask for single or double reads


def verify_main(argv=None):
    parser = argparse.ArgumentParser(description="Verify md5 sums for an archived run")
    parser.add_argument(
        "--archive-dir",
        required=True,
        type=Path,
        help="Archive directory containing the md5 checksum file and reads.",
    )
    args = parser.parse_args(argv)

    return verify_archive(args.archive_dir)
