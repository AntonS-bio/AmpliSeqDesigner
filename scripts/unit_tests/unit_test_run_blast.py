from os.path import expanduser, realpath, dirname
import unittest
from typing import List
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration, BlastResult
import metadata_utils as metadata_utils
from run_blast import BlastRunner

class TestBlastRunner(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    temp_dir=expanduser("~/HandyAmpliconTool/unit_test_data/temp_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    ref_fasta=f'{valid_data}/GCF_000195995.1_short.fna'
    long_fasta=f'{valid_data}/GCF_000195995.1_for_vcf.fna'
    fasta_seq="GAAATAGCCTGCTGATAGAGACTTTCATTCTCGGTTCCAGAGCGTTGTTGCAGTGCAGGATAAATAAAGGAGTAAAG"

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        return super().setUp()
    
    def test_db_from_file(self):
        blast_runner=BlastRunner()
        self.assertTrue(blast_runner.db_from_file(self.ref_fasta, self.temp_dir))
        self.assertTrue(blast_runner.db_from_file(self.ref_fasta, self.temp_dir))
        self.assertRaises(OSError,blast_runner.db_from_file, file_name=self.config_file, db_dir=self.temp_dir)
        self.assertRaises(ValueError,blast_runner.db_from_file, file_name="NoneSuch", db_dir=self.temp_dir)
        self.assertRaises(ValueError,blast_runner.db_from_file, file_name="NoneSuch", db_dir=self.temp_dir)

    def test_db_from_string(self):
        blast_runner=BlastRunner()
        fasta_header="RandomSeq"
        self.assertTrue(blast_runner.db_from_string(fasta_header, self.fasta_seq, self.temp_dir))
        self.assertTrue(blast_runner.db_from_string(fasta_header, self.fasta_seq, self.temp_dir))
        

    def test_run_from_string(self):
        blast_runner=BlastRunner()
        blast_runner.db_from_file(self.long_fasta, db_dir=self.config_data.temp_blast_db)
        hits: List[BlastResults] = blast_runner.run_from_string("test_seq", self.fasta_seq, self.temp_dir )
        self.assertTrue(len(hits)==1)
        self.assertTrue(hits[0].qstart==1)
        self.assertTrue(hits[0].qend==77)
        self.assertTrue(hits[0].pident==100.0)
        self.assertTrue(hits[0].is_flipped==False)
        print(hits)


if __name__ == '__main__':
    unittest.main(verbosity=2)