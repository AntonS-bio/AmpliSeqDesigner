from os.path import expanduser, realpath, dirname
from typing import List
import unittest
from sys import path
unit_test_dir = dirname(realpath(__file__))
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)
from data_classes import InputConfiguration, Amplicon, BlastResult
from generate_msa import MsaGenerator

class TestMsaGenenerator(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")
    ref_fasta=f'{valid_data}/GCF_000195995.1_short.fna'
    long_fasta=f'{valid_data}/GCF_000195995.1_for_vcf.fna'
    acr_seq="CAACCTTGTTTTTTTCGCCTGGACGATCGGCCCAGTCTTTCAACGACACAAATGCAATACCGGTATTCTGACCGC"

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        return super().setUp()
    
    @property
    def dummy_amplicon(self) -> Amplicon:
        return Amplicon("Test Amplicon",self.acr_seq)
    
    @property
    def dummy_blast_results(self) -> List[BlastResult]:
        first_result=BlastResult.from_blast_line("qseq_1\t523326\t523400\tsseqID\t1\t75\t100.0\t2.52e-34\tCAACCTTGTTTTTTTCGCCTGGACGATCGGCCCAGTCTTTCAACGACACAAATGCAATACCGGTATTCTGACCGC")
        return [first_result]
    
    def test_get_fasta_files(self):
        generator=MsaGenerator(self.config_data.temp_blast_db)
        generator._get_fasta_files(self.valid_data)
        self.assertEqual([f'{self.valid_data}/GCF_000195995.1_short.fna',
                          f'{self.valid_data}/GCF_000195995.1_for_vcf.fna',
                          f'{self.valid_data}/existing_primers.fasta'],generator.file_to_search)
        #self.assertEqual()

    def test_run_blast(self):
        generator=MsaGenerator(self.config_data.temp_blast_db)
        results: List[BlastResult]=generator._run_blast([self.dummy_amplicon],[self.ref_fasta,self.long_fasta])
        expected_results=self.dummy_blast_results
        #has to be done in two steps because the Pool running results doesn't return them in specific order
        self.assertTrue(len(results)==1)
        self.assertTrue( results[0].coordinates_match(expected_results[0])  )


    def test_process_blast_results(self):
        generator=MsaGenerator(self.config_data.temp_blast_db)
        amplicon=self.dummy_amplicon
        expected_results=self.dummy_blast_results
        for result in expected_results:
            result.sseqid=amplicon._uuid
        valid_amplicon_hits=generator._process_blast_results(expected_results, [amplicon])
        for key, values in valid_amplicon_hits.items():
            for value in values:
                self.assertTrue(value.coordinates_match(expected_results[0]) or value.coordinates_match(expected_results[1]))

    def test_generate_msa(self):
        generator=MsaGenerator(self.config_data.temp_blast_db)
        amplicon=self.dummy_amplicon
        a=generator.generate_msa( [amplicon], genomes_dir=self.config_data.negative_genomes)
        self.assertEqual("".join([f for f in a[amplicon._uuid].iloc[0]]), amplicon.seq)

    def test_align_results_helper(self):
        generator=MsaGenerator(self.config_data.temp_blast_db)
        dummy_results=self.dummy_blast_results
        dummy_amplicon=self.dummy_amplicon
        result=generator._align_results_helper([dummy_results, dummy_amplicon._uuid, dummy_amplicon.seq])
        self.assertTrue(dummy_amplicon._uuid in result.keys())
        self.assertTrue(str.lower(dummy_amplicon.seq) == result[dummy_amplicon._uuid])
        for blast_result in dummy_results:
            self.assertTrue(blast_result.qseqid in result.keys())
            self.assertTrue(str.lower(blast_result.qseq) == result[blast_result.qseqid])
        print(result)

    def test_msa_to_dataframe(self):
        generator=MsaGenerator(self.config_data.temp_blast_db)
        amplicon=self.dummy_amplicon
        result=generator._msa_to_dataframe({amplicon._uuid: amplicon.seq})
        self.assertEqual( "".join([f for f in result.iloc[0]]) , amplicon.seq)

if __name__ == '__main__':
    unittest.main(verbosity=2)