from os.path import expanduser, realpath, dirname
unit_test_dir = dirname(realpath(__file__))
from sys import path
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)

import unittest
from data_classes import InputConfiguration, FlankingAmplicon, Amplicon, BlastResult

class TestFileValidators(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        return super().setUp()
    
    
    def test_snp(self) -> None:
        #amplicon=Amplicon.from_bed_line("AL513382_143N_pHCM1_120N_pHCM2\t100\t200\tTest Amplicon\n", self.config_data.reference_fasta)
        self.assertRaises(ValueError,Amplicon.from_bed_line, 
                          bed_line="AL513382_143N_pHCM1_120N_pHCM2\t100\t10000000000\tTest Amplicon\n",
                          ref_fasta_file=self.config_data.reference_fasta)
        self.assertRaises(ValueError,Amplicon.from_bed_line, 
                          bed_line="NoneSuch\t100\t1000\tTest Amplicon\n",
                          ref_fasta_file=self.config_data.reference_fasta)
        
    def test_blastresult(self) -> None:
        result=BlastResult.from_blast_line("3.1.1_1_669587_R\t2\t20\tNC_003198.1\t1663171\t1663189\t89.474\t0.16\tTCTGGTACAAAGCGGCGAA\n")
        self.assertEqual(result.value, "3.1.1_1_669587_R 2 20 NC_003198.1 1663171 1663189")

    def test_from_parent_bed_line(self) -> None:
        amplicon=Amplicon.from_bed_line("AL513382_143N_pHCM1_120N_pHCM2\t100\t200\tTest Amplicon\n", self.config_data.reference_fasta)
        flanking_amplicon=FlankingAmplicon.from_parent_bed_line(self.config_data.reference_fasta,True,self.config_data.flank_len_to_check, amplicon)
        self.assertTrue(flanking_amplicon.seq=="AGAGATTACGTCTGGTTGCAAGAGATCATAACAGGGGAAATTGATTGAAAATAAATATATCGCCAGCAGCACATGAACAAGTTTCGGAATGTGATCAATT")
        self.assertTrue(flanking_amplicon.name=="AL513382_143N_pHCM1_120N_pHCM2_100_200_Test Amplicon_left")

if __name__ == '__main__':
    unittest.main(verbosity=2)