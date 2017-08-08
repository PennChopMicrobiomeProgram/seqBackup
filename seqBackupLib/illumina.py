import gzip
import os.path
import re
import warnings


class IlluminaFastq(object):
    machine_types = {"D": "Illumina-HiSeq", "M": "Illumina-MiSeq"}

    def __init__(self, f):
        self.file = f
        self.fastq_info = self._parse_header()
        self.folder_info = self._parse_folder()
        
    def __str__(self):
        return "_".join([self.fastq_info["instrument"], 
                         self.fastq_info["run_number"], 
                         self.fastq_info["flowcell_id"],
                         self.fastq_info["lane"]])

    def is_same_run(self, other):
        run_check = self.fastq_info["run_number"] == other.fastq_info["run_number"]
        instrument_check = self.fastq_info["instrument"] == other.fastq_info["instrument"]
        flowcell_check = self.fastq_info["flowcell_id"] == other.fastq_info["flowcell_id"]
        return (run_check and instrument_check and flowcell_check)

    def _parse_header(self):
        line = next(self.file).strip()
        if not line.startswith("@"):
            raise ValueError("Not a FASTQ header line")
        # Remove first character, @
        line = line[1:]
        word1, _, word2 = line.partition(" ")
        
        keys1 = ("instrument", "run_number", "flowcell_id", "lane")
        vals1 = dict((k, v) for k, v in zip(keys1, word1.split(":")))

        keys2 = ("read", "is_filtered", "control_number", "index_reads")
        vals2 = dict((k, v) for k, v in zip(keys2, word2.split(":")))

        vals1.update(vals2)
        return vals1

    def _parse_folder(self):
        matches = re.match("(\\d{6})_([DM]\\d{5})_0*(\\d{1,4})_(.*)", self.run_name)
        keys1 = ("date", "instrument", "run_number", "flowcell_id")
        vals1 = dict((k, v) for k, v in zip(keys1, matches.groups()))

        if self.machine_type == "Illumina-HiSeq":
            vals1["flowcell_id"] =  vals1["flowcell_id"][1:]

        matches = re.match("Undetermined_S0_L00([1-8])_([RI])([12])_001.fastq.gz", os.path.basename(self.filepath))
        keys2 = ("lane", "read_or_index", "read")
        vals2 = dict((k, v) for k, v in zip(keys2, matches.groups()))
        
        vals1.update(vals2)
        return vals1

    @property
    def machine_type(self):
        instrument_code = self.fastq_info["instrument"][0]
        return self.machine_types[instrument_code]

    @property
    def date(self):
        year = self.run_name[0:2]
        month = self.run_name[2:4]
        day = self.run_name[4:6]
        return "20{0}-{1}-{2}".format(year, month, day)

    @property
    def lane(self):
        return self.fastq_info["lane"]

    @property
    def filepath(self):
        return self.file.name


    @property
    def run_name(self):
        dir_split = self.filepath.split(os.sep)
        matches = [re.match("\\d{6}_[DM]\\d{5}_\\d{4}", d) for d in dir_split]
        matches = [dir_split[i] for i, m in enumerate(matches) if m]
        if len(matches) != 1:
            raise ValueError("Could not find run name in directory: {0}".format(self.filepath))
        return matches[0]

    def build_archive_dir(self):
        return '_'.join([self.run_name, 'L{:0>3}'.format(self.lane)])

    def check_fp_vs_content(self):
        run_check = self.fastq_info["run_number"] == self.folder_info["run_number"]
        instrument_check = self.fastq_info["instrument"] == self.folder_info["instrument"]
        flowcell_check = self.fastq_info["flowcell_id"] == self.folder_info["flowcell_id"]
        lane_check = self.lane == self.folder_info["lane"]
        read_check = self.fastq_info["read"] == self.folder_info["read"]
        return (run_check and instrument_check and flowcell_check and lane_check and read_check)
    
    def check_file_size(self, min_file_size):
        return os.path.getsize(self.filepath) > min_file_size

    def check_index_read_exists(self):
        return len(self.fastq_info["index_reads"]) > 2
