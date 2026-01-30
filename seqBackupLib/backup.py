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
    allow_check_failures: bool = False,
):

    R1 = IlluminaFastq(gzip.open(forward_reads, mode="rt"))

    # build the strings for the required files
    RI_fps = build_fp_to_archive(forward_reads, has_index, R1.lane)

    # create the Illumina objects and check the files
    illumina_fastqs = [IlluminaFastq(gzip.open(fp, mode="rt")) for fp in RI_fps]
    r1 = illumina_fastqs[0]

    fp_vs_content_results = [ifq.check_fp_vs_content()[0] for ifq in illumina_fastqs]
    if not all(fp_vs_content_results):
        [ifq.check_fp_vs_content(verbose=True) for ifq in illumina_fastqs]
        message = (
            "The file path and header information don't match",
            [
                str(ifq)
                for ifq, ok in zip(illumina_fastqs, fp_vs_content_results)
                if not ok
            ],
        )
        if allow_check_failures:
            warnings.warn(f"{message[0]}: {message[1]}")
        else:
            raise ValueError(*message)
    file_size_results = [ifq.check_file_size(min_file_size) for ifq in illumina_fastqs]
    if not all(file_size_results):
        message = (
            "File seems suspiciously small. Please check if you have the correct file or"
            " lower the minimum file size threshold",
            file_size_results,
        )
        if allow_check_failures:
            warnings.warn(f"{message[0]}: {message[1]}")
        else:
            raise ValueError(*message)
    if not all([ifq.check_index_read_exists() for ifq in illumina_fastqs]):
        warnings.warn(
            "No barcodes in headers. Were the fastq files generated properly?"
        )

    # parse the info from the headers in EACH file and check they are consistent within each other
    same_run_results = [
        fastq.is_same_run(illumina_fastqs[0]) for fastq in illumina_fastqs
    ]
    if not all(same_run_results):
        message = "The files are not from the same run."
        if allow_check_failures:
            warnings.warn(message)
        else:
            raise ValueError(message)

    ## Archiving steps

    # make sure the sample sheet exists
    if not sample_sheet_fp.is_file():
        raise IOError("Sample sheet does not exist", str(sample_sheet_fp))

    # create the folder to write to
    write_dir = dest_dir / r1.build_archive_dir()
    write_dir.mkdir(parents=True, exist_ok=False)

    ### All the checks are done and the files are safe to archive!

    # move the files to the archive location and remove permission
    permission = stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
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
    parser.add_argument(
        "--allow-check-failures",
        action="store_true",
        help="Continue archiving even if validation checks fail",
    )
    args = parser.parse_args(argv)
    return backup_fastq(
        args.forward_reads,
        args.destination_dir,
        args.sample_sheet,
        not args.no_index,
        args.min_file_size,
        args.allow_check_failures,
    )

    # maybe also ask for single or double reads
