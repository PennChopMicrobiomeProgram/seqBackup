import argparse
import os
import shutil
import stat
import re
import gzip
import hashlib
import warnings

from seqBackupLib.illumina import IlluminaFastq

def build_fp_to_archive(file_name, has_index, lane):

    if re.search("R1_001.fastq", file_name) is None:
        raise IOError("The file doesn't look like an R1 file: {}".format(file_name))

    label = ["R2"]
    if has_index:
        label.extend(["I1", "I2"])

    rexp = "".join(["(L00", lane, "_)(R1)(_001.fastq.gz)$"])
    modified_fp = [re.sub(rexp, "".join(["\\1", lab, "\\3"]), file_name) for lab in label]
    return [file_name] + modified_fp

def return_md5(fname):
    # from https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def backup_fastq(forward_reads, dest_dir, sample_sheet_fp, has_index, min_file_size):
    
    R1 = IlluminaFastq(gzip.open(forward_reads, mode = 'rt'))    

    # build the strings for the required files    
    file_names_RI = build_fp_to_archive(forward_reads, has_index, R1.lane)

    # create the Illumina objects and check the files
    illumina_fastqs = []
    for fp in file_names_RI:
        illumina_temp = IlluminaFastq(gzip.open(fp, mode = 'rt'))
        if not illumina_temp.check_fp_vs_content()[0]:
            print(illumina_temp.check_fp_vs_content()[1:])
            raise ValueError("The file path and header information don't match")
        if not illumina_temp.check_file_size(min_file_size):
            raise ValueError("File {0} seems suspiciously small. Plese check if you have the correct file or lower the minimum file size threshold".format(fp))
        if not illumina_temp.check_index_read_exists():
            warnings.warn("No barcodes in headers. Were the fastq files generated properly?: {0}".format(fp))
        illumina_fastqs.append(illumina_temp)

    # parse the info from the headers in EACH file and check they are consistent within each other
    if not all([fastq.is_same_run(illumina_fastqs[0]) for fastq in illumina_fastqs]):
        raise ValueError("The files are not from the same run.")

    ## Archiving steps

    # make sure the sample sheet exists
    if not os.path.isfile(sample_sheet_fp):
        raise IOError("Sample sheet does not exist: {}".format(sample_sheet_fp))

    # create the folder to write to
    write_dir = os.path.join(dest_dir, illumina_temp.build_archive_dir())

    # create the folder. If it exists exit
    if os.path.isdir(write_dir):
        raise IOError("The folder already exists: {}".format(write_dir))
    os.mkdir(write_dir)

    ### All the checks are done and the files are safe to archive!

    # move the files to the archive location and remove permission
    permission = stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH
    for fp in file_names_RI:
        shutil.copyfile(fp, os.path.join(write_dir, os.path.basename(fp)))
        os.chmod(os.path.join(write_dir, os.path.basename(fp)), permission) #this doesn't work on isilon

    # copy the sample sheet to destination folder
    shutil.copyfile(sample_sheet_fp, os.path.join(write_dir, os.path.basename(sample_sheet_fp)))

    # write md5sums to a file
    md5s = [(os.path.basename(fp), return_md5(fp)) for fp in file_names_RI]
    md5out_fp = os.path.join(write_dir, ".".join([illumina_temp.build_archive_dir(), "md5"]))
    with open(md5out_fp, "w") as md5_out:
        [md5_out.write("\t".join(md5) + "\n") for md5 in md5s]

def main(argv=None):
    parser = argparse.ArgumentParser(description="Backs up fastq files")

    parser.add_argument(
        "--forward-reads", required=True,
        type=str,
        help="R1.fastq")
    parser.add_argument(
        "--destination-dir", required=True,
        type=str,
        help="Destination folder to copy the files to.")
    parser.add_argument(
        "--sample-sheet", required=True,
        type=str,
        help="The sample sheet associated with the run.")
    parser.add_argument(
        "--has-index", required=False,
        type=bool, default=True,
        help="Are index reads generated")
    parser.add_argument(
        "--min-file-size", required=False,
        type=int, default=500000000,
        help="Minimum file size to register in bytes")
    args = parser.parse_args(argv)

    backup_fastq(args.forward_reads, args.destination_dir, args.sample_sheet, args.has_index, args.min_file_size)

    # maybe also ask for single or double reads
