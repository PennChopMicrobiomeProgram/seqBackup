import argparse
import os
import re

def find_run_name(forward_reads):
    dir_split = forward_reads.split(os.sep)
    matches = [re.search("\\d{6}_[DM]\\d{5}_\\d{4}", d) for d in dir_split]
    matches = [dir_split[i] for i, m in enumerate(matches) if m]
    if len(matches) != 1:
        raise ValueError("Could not find run name in directory: {0}".format(forward_reads))
    return matches[0]


def backup_fastq(forward_reads, num_lanes, run_name, dest_dir, has_index):
    

    # parse out the run name
    if run_name is None:
        run_name = find_run_name(forward_reads)

    
    # parse out the file name
    file_name = forward_reads.split(os.sep)[-1]


    # build the strings for the required files
    


    # check if the files are actually there

    # check their size. If less than 500 Mb quit with error

    
    # parse the info from the headers in EACH file and check they are consistent within each other

    # check to make sure the run name matches the folder name


    # check to make sure they have barcode info





    # build the destination folder

    # create the folder. If it exists exit


    # move the files there



    # move the incoming folder to a subfolder called "imported"





def main(argv=None):
    parser = argparse.ArgumentParser(description="Backs up fastq files")

    parser.add_argument(
        "--forward-reads", required=True,
        type=argparse.FileType("r"),
        help="R1.fastq")
    parser.add_argument(
        "--number-of-lanes", required=True,
        type=int,
        help="Number of lanes sequenced")
    parser.add_argument(
        "--run-name", required=False,
        type=str,
        help="The run name to backup as. Defaults to the Illumina file.")
    parser.add_argument(
        "--destination-dir", required=False,
        type=str,
        help="Destination folder to copy the files to.")
    parser.add_argument(
        "--has-index", required=False,
        type=bool,
        help="Are index reads generated")
    args = parser.parse_args(argv)

    # Check if the R1 file exists
    fwd_fp = args.forward_reads.name
    args.forward_reads.close()

    backup_fastq(args.forward_reads.name, args.number_of_lanes, args.run_name, args.destination_dir, args.has_index)

    # maybe also ask for single or double reads
