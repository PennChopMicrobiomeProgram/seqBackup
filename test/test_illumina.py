import unittest
import os
import gzip
from io import StringIO

from seqBackupLib.illumina import IlluminaFastq

class IlluminaTests(unittest.TestCase):
    
    def test_illuminafastq(self):
        fastq_file = StringIO(
            u"@M03543:47:C8LJ2ANXX:1:2209:1084:2044 1:N:0:NNNNNNNN+NNNNNNNN")
        fastq_filepath = (
            "Miseq/160511_M03543_0047_000000000-APE6Y/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        folder_info = {"date":"160511", "instrument":"M03543", "run_number":"47", "flowcell_id":"000000000-APE6Y", "lane":"1", "read_or_index":"R", "read":"1"}
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)

        self.assertEqual(fq.machine_type, "Illumina-MiSeq")
        self.assertEqual(fq.date, "2016-05-11")
        self.assertEqual(fq.lane, "1")
        self.assertEqual(fq.filepath, fastq_filepath)
        self.assertEqual(fq.run_name, "160511_M03543_0047_000000000-APE6Y")

        self.assertDictEqual(fq.folder_info, folder_info)

        fastq_file = StringIO(
            u"@D00727:27:CA7HHANXX:1:1105:1243:1992 1:N:0:NGATCAGT+NNAAGGAG")
        fastq_filepath = (
            "Hiseq/170330_D00727_0027_ACA7HHANXX/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        folder_info = {"date":"170330", "instrument":"D00727", "run_number":"27", "flowcell_id":"CA7HHANXX", "lane":"1", "read_or_index":"R", "read":"1"}
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)

        self.assertEqual(fq.machine_type, "Illumina-HiSeq")
        self.assertEqual(fq.date, "2017-03-30")
        self.assertEqual(fq.lane, "1")
        self.assertEqual(fq.filepath, fastq_filepath)
        self.assertEqual(fq.run_name, "170330_D00727_0027_ACA7HHANXX")

        self.assertDictEqual(fq.folder_info, folder_info)

    def test_fp_vs_content(self):
        # check correct case for Miseq data
        fastq_file = StringIO(
            u"@M04734:28:000000000-B2MVT:1:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertTrue(fq.check_fp_vs_content())

        # check correct case for Hiseq data
        fastq_file = StringIO(
            u"@D00727:27:CA7HHANXX:1:1105:1243:1992 1:N:0:NGATCAGT+NNAAGGAG")
        fastq_filepath = (
            "Hiseq/170330_D00727_0027_ACA7HHANXX/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertTrue(fq.check_fp_vs_content())

        # case when the lane number doesn't match
        fastq_file = StringIO(
            u"@M04734:28:000000000-B2MVT:3:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertFalse(fq.check_fp_vs_content())

        # case when the flow cell ID doesn't match
        fastq_file = StringIO(
            u"@M04734:28:000000000-BBBBB:1:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertFalse(fq.check_fp_vs_content())

        # case when the machine doesn't match
        fastq_file = StringIO(
            u"@D04734:28:000000000-BBBBB:1:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertFalse(fq.check_fp_vs_content())

        # case when the read doesn't match
        ### important: It won't distinguish between R1 and I1.
        fastq_file = StringIO(
            u"@M04734:28:000000000-B2MVT:1:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertFalse(fq.check_fp_vs_content())

    def test_check_index_read_exists(self):
        # test passing
        fastq_file = StringIO(
            u"@M04734:28:000000000-B2MVT:1:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertTrue(fq.check_index_read_exists())
        
        # test failing
        fastq_file = StringIO(
            u"@M04734:28:000000000-B2MVT:1:2106:17605:1940 1:N:0:0")
        fastq_filepath = (
            "Miseq/170323_M04734_0028_000000000-B2MVT/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertFalse(fq.check_index_read_exists())

    def test_build_archive_dir(self):
        # for MiSeq
        fastq_file = StringIO(
            u"@M03543:47:C8LJ2ANXX:1:2209:1084:2044 1:N:0:NNNNNNNN+NNNNNNNN")
        fastq_filepath = (
            "Miseq/160511_M03543_0047_000000000-APE6Y/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertEqual(fq.build_archive_dir(), "160511_M03543_0047_000000000-APE6Y_L001")

        # for HiSeq
        fastq_file = StringIO(
            u"@D00727:27:CA7HHANXX:1:1105:1243:1992 1:N:0:NGATCAGT+NNAAGGAG")
        fastq_filepath = (
            "Hiseq/170330_D00727_0027_ACA7HHANXX/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        self.assertEqual(fq.build_archive_dir(), "170330_D00727_0027_ACA7HHANXX_L001")

    def test_check_file_size(self):
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        fastq_filepath = os.path.join(curr_dir, "170323_M04734_0028_000000000-B2MVT/Undetermined_S0_L001_R1_001.fastq.gz")
        fq = IlluminaFastq(gzip.open(fastq_filepath, mode = 'rt'))
        self.assertTrue(fq.check_file_size(50))
        self.assertFalse(fq.check_file_size(50000))

    def test_is_same_run(self):
        fastq_file = StringIO(
            u"@M03543:47:C8LJ2ANXX:1:2209:1084:2044 1:N:0:NNNNNNNN+NNNNNNNN")
        fastq_filepath = (
            "Miseq/160511_M03543_0047_000000000-APE6Y/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq = IlluminaFastq(fastq_file)
        fastq_file.seek(0)
        fq1 = IlluminaFastq(fastq_file)

        fastq_file = StringIO(
            u"@D00727:27:CA7HHANXX:1:1105:1243:1992 1:N:0:NGATCAGT+NNAAGGAG")
        fastq_filepath = (
            "Hiseq/170330_D00727_0027_ACA7HHANXX/Data/Intensities/"
            "BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz")
        fastq_file.name = fastq_filepath
        fq2 = IlluminaFastq(fastq_file)

        self.assertTrue(fq.is_same_run(fq1))
        self.assertFalse(fq.is_same_run(fq2))
