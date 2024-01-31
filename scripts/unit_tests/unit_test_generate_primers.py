from os.path import expanduser, realpath, dirname
from typing import List
import unittest
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration, Amplicon, BlastResult
from primers_generator import PrimersGenerator



class TestPrimersGenerator(unittest.TestCase):

    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    long_fasta=f'{valid_data}/GCF_000195995.1_for_vcf.fna'
    existing_primers_bed=f'{valid_data}/test_primer.bed'
    invalid_existing_primers_bed=f'{invalid_data}/test_primer.bed'
    existing_primers_fasta=f'{valid_data}/existing_primers.fasta'
    primer_seq_in_ref_direct="GCGGGCTTTACTGCCGGTAATGAAAAGG"
    primer_seq_in_ref_revcomp="CCTTTTCATTACCGGCAGTAAAGCCCGC"
    species_snps=f'{valid_data}/species_snps.pkl'
    gt_snps=f'{valid_data}/gt_snps.pkl'

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        self.generator=PrimersGenerator(self.config_data)
        return super().setUp()

    def test_get_seq_coordinates_in_ref(self):
        self.assertEqual(self.generator._get_seq_coordinates_in_ref(self.primer_seq_in_ref_direct), 897)
        self.assertEqual(self.generator._get_seq_coordinates_in_ref(self.primer_seq_in_ref_revcomp), 897)
        self.assertEqual(self.generator._get_seq_coordinates_in_ref("NoneSuch"), -1)

    def test_sequence_to_primer(self):
        new_primer=self.generator._sequence_to_primer(self.primer_seq_in_ref_direct, False)
        self.assertEqual(new_primer.length, 28)
        self.assertEqual(new_primer.ref_start, 897)
        self.assertEqual(new_primer.seq, self.primer_seq_in_ref_direct)

    def test_load_existing_primer_from_fasta(self):
        self.generator.load_existing_primer_from_fasta(self.existing_primers_fasta)
        self.assertTrue( len(self.generator.existing_primers)==4 )


    def test_load_existing_primer_from_bed(self):
        self.assertRaises(ValueError, self.generator.load_existing_primer_from_bed, existing_primers_bed=self.invalid_existing_primers_bed)
        self.generator.load_existing_primer_from_bed(self.existing_primers_bed)
        self.assertTrue( len(self.generator.existing_primers)==4 )

    def test_has_homodimers(self):
        self.assertTrue(self.generator.forms_homodimers("GCGGGCTTTACTGCCGGTAATGAAAAGG"))
        self.assertFalse( self.generator.forms_homodimers("TGCAACATGAAGGTGACGATG") )

    def test_count_heterodimers(self):
        self.generator.load_existing_primer_from_bed(self.existing_primers_bed)
        self.generator.count_heterodimers(self.primer_seq_in_ref_direct,"Forward")
        #### ADD Further checks here

    def test_for_testing_load_gts(self):
        self.generator._for_testing_load_gts(self.species_snps, self.gt_snps)

    def test_find_candidate_primers(self):
        self.generator._for_testing_load_gts(self.species_snps, self.gt_snps)
        self.generator.find_candidate_primers()


if __name__ == '__main__':
    unittest.main(verbosity=2)