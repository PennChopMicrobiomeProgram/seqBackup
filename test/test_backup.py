import unittest
import gzip
import tempfile
from io import StringIO

from seqBackupLib.illumina import IlluminaFastq
#from seqBackupLib.backup import (build_fp_to_archive, return_md5, backup_fastq)
from seqBackupLib.backup import *

class BackupTests(unittest.TestCase):
    def setUp(self):
        self.curr_dir = os.path.dirname(os.path.abspath(__file__))
        self.fastq_filepath = fastq_filepath = os.path.join(self.curr_dir, "170323_M04734_0028_000000000-B2MVT/small_Undetermined_S0_L001_R1_001.fastq.gz")
        self.temp_out_dir = tempfile.mkdtemp(dir=self.curr_dir)

    #def tearDown(self):
    #    shutil.rmtree(self.temp_out_dir)

    def test_build_fp_to_archive(self):
        list1 = build_fp_to_archive("Undetermined_S0_L001_R1_001.fastq.gz", True, "1")
        self.assertItemsEqual(list1, ["Undetermined_S0_L001_R1_001.fastq.gz", "Undetermined_S0_L001_R2_001.fastq.gz", "Undetermined_S0_L001_I1_001.fastq.gz", "Undetermined_S0_L001_I2_001.fastq.gz"])
        
        list1 = build_fp_to_archive("Undetermined_S0_L001_R1_001.fastq.gz", False, "1")
        self.assertItemsEqual(list1, ["Undetermined_S0_L001_R1_001.fastq.gz", "Undetermined_S0_L001_R2_001.fastq.gz"])

    def test_backup_fastq(self):
        has_index = True
        min_file_size = 5
        backup_fastq(self.fastq_filepath, self.temp_out_dir, has_index, min_file_size)
        # check the md5sums of the first fastq is the same
        md5_orj = return_md5(self.fastq_filepath)
        md5_trans = return_md5(os.path.join(self.temp_out_dir, os.path.basename(self.fastq_filepath)))
        self.assertEqual(md5_orj, md5_trans)
        
        # check write permissions of the files
        
    
    def test_return_md5(self):
        self.assertEqual(return_md5(self.fastq_filepath), "13695e47114c02536ae3ca6823a42261")
