from os.path import expanduser, realpath, dirname
from typing import List
import unittest
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from Bio.Seq import Seq
from data_classes import InputConfiguration, Amplicon, BlastResult, PrimerPair, Primer
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

    # def test_for_testing_load_gts(self):
    #     self.generator._for_testing_load_gts(self.species_snps, self.gt_snps)


    def get_dummy_primer_pairs(self) -> List[PrimerPair]:
        pairs: List[PrimerPair] = []
        for i, sequences in enumerate(zip( ["GGCAT", "GATC", "GATC"], ["CTAG","CTAG", "CTAG"] )):
            forward_primer=Primer(sequences[0],50.0, 50.0, False)
            reverse_primer=Primer(sequences[1],50.0, 50.0, False)
            pairs.append( PrimerPair(f'Pair_{i+1}', forward_primer, reverse_primer ) )
        return pairs

    def test_remove_duplicate_primer_pairs(self):
        pairs=self.get_dummy_primer_pairs()
        self.generator._remove_duplicate_primer_pairs(pairs)
        self.assertEqual(len(pairs),2)

    def test_both_primers_given(self):
        pass

    def test_right_primers_given(self):
        pass

    def test_left_primers_given(self):
        pass

    def test_process_p3_output(self):
        pass

    def test_primer_pairs_list_to_string(self):
        pass
        # pairs=self.get_dummy_primer_pairs()
        # primers_string=self.generator._primer_pairs_list_to_string(pairs)
        # expected_result=">Pair_1_Forward_\nGGCAT\n>Pair_1_Reverse\nCTAG\n>Pair_2_Forward_\nGATC\n>Pair_2_Reverse\nCTAG\n>Pair_3_Forward_\nGATC\n>Pair_3_Reverse\nCTAG\n"
        # self.assertEqual(primers_string, expected_result)

    def test_primers_list_to_string(self):
        pairs=self.get_dummy_primer_pairs()
        a, b= list(zip(*[(f.forward, f.reverse) for f in pairs ]))
        primers_string=self.generator._primers_list_to_string(a+b)
        expected_result='>Primer_1\nGGCAT\n>Primer_2\nGATC\n>Primer_3\nGATC\n>Primer_4\nCTAG\n>Primer_5\nCTAG\n>Primer_6\nCTAG\n'
        self.assertEqual(primers_string, expected_result)
        
    def test_remove_interfering_primers(self):
        self.generator.config.flank_len_to_check=2000
        self.generator.load_existing_primer_from_bed(self.existing_primers_bed)
        self.generator.new_primer_pairs.clear()
        primer_f=Primer("CATTTCGGGTGATTCTGTTATCTGTGTCACACTTTT",50.0,50.0,False)
        primer_r=Primer("CCACAACCAACAGATTGATAAA",50.0,50.0,False)
        interfering_pair=PrimerPair("interfering", primer_f,primer_r)
        self.generator.new_primer_pairs.append(interfering_pair)
        self.assertEqual(len(self.generator.new_primer_pairs) , 1 )
        self.generator._remove_interfering_primers(self.generator.new_primer_pairs)
        self.assertEqual(len(self.generator.new_primer_pairs) , 0 )

    # def test_remove_primers_in_repeat_regions(self):
    #     pass
    #     self.generator.load_existing_primer_from_bed(self.existing_primers_bed)
    #     self.generator._for_testing_load_gts(self.species_snps, self.gt_snps)
    #     self.generator.find_candidate_primers()
    #     self.generator._remove_primers_in_repeat_regions(self.generator.new_primer_pairs)
        #Add comparison later

    def test_find_candidate_primers(self):
        self.generator._for_testing_load_gts(self.species_snps, self.gt_snps)
        a=self.generator.find_candidate_primers()
        with open(self.config_data.temp_blast_db+"temp.tsv","w") as output_file:
            header="\t".join(["Name", "Penalty", "Contig", "Start","End","Length",
                            "Forward","Forward Tm", "Forward GC",
                            "Reverse","Reverse Tm", "Reverse GC"])+"\n"
            output_file.write(header)
            for pair in a:
                output_file.write(pair.to_string()+"\n")
        for primer in self.generator.new_primer_pairs:
            self.assertEqual(primer.forward.seq, self.generator.ref_seq[primer.ref_contig][primer.forward.ref_start:primer.forward.ref_end])
            rev_comp=Seq(self.generator.ref_seq[primer.ref_contig][primer.reverse.ref_start:primer.reverse.ref_end]).reverse_complement()
            self.assertEqual(primer.reverse.seq, str(rev_comp))



if __name__ == '__main__':
    unittest.main(verbosity=2)