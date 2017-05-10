import argparse
import os
import re
import gzip

from seqBackupLib.illumina import IlluminaFastq

def build_fp_to_archive(file_name, has_index, lane):
    label = ["R2"]
    if has_index:
        label.extend(["I1", "I2"])

    rexp = "".join(["(L00", lane, "_)(R1)(_001.fastq.gz)$"])
    modified_fp = [re.sub(rexp, "".join(["\\1", lab, "\\3"]), file_name) for lab in label]
    return [file_name] + modified_fp

def backup_fastq(forward_reads, dest_dir, has_index, min_file_size):
    
    R1 = IlluminaFastq(gzip.GzipFile(forward_reads))    

    # build the strings for the required files    
    file_name_RI = build_fp_to_archive(forward_reads, has_index, R1.lane)
    
    # create the Illumina objects and check the files
    illumina_fastqs = []
    for fp in file_name_RI:
        illumina_temp = IlluminaFastq(gzip.GzipFile(fp))
        illumina_temp.check_fp_vs_content()
        illumina_temp.check_file_size(min_file_size)
        illumina_temp.check_index_read_exists()
        illumina_fastqs.append(str(illumina_temp))

    # parse the info from the headers in EACH file and check they are consistent within each other
    if not all([fastq == illumina_fastqs[0] for fastq in illumina_fastqs]):
        raise ValueError("The files are not from the same run.")


    ## Everything to do with the destination folder

    # create the folder. If it exists exit


    # move the files there



    # move the incoming folder to a subfolder called "imported"


    # write md5sums to a file


    # remove write permission from file



def main(argv=None):
    parser = argparse.ArgumentParser(description="Backs up fastq files")

    parser.add_argument(
        "--forward-reads", required=True,
        type=str,
        help="R1.fastq")
    #parser.add_argument(
    #    "--lane-num", required=True,
    #    type=int,
    #    help="Number of lanes sequenced")
    #parser.add_argument(
    #    "--run-name", required=False,
    #    type=str,
    #    help="The run name to backup as. Defaults to the Illumina file.")
    parser.add_argument(
        "--destination-dir", required=True,
        type=str,
        help="Destination folder to copy the files to.")
    parser.add_argument(
        "--has-index", required=False,
        type=bool, default=True,
        help="Are index reads generated")
    parser.add_argument(
        "--min-file-size", required=False,
        type=int, default=500000000,
        help="Minimum file size to register in bytes")
    args = parser.parse_args(argv)

    backup_fastq(args.forward_reads, args.destination_dir, args.has_index, args.min_file_size)

    # maybe also ask for single or double reads
