from os.path import expanduser, realpath, dirname, exists
from os import mkdir
unit_test_dir = dirname(realpath(__file__))
from sys import path
path.append( f'{unit_test_dir}/..')
print(unit_test_dir)

import unittest
from identify_species_snps import IdentifySpeciesSnps
import name_converters
from data_classes import InputConfiguration, Genotype, Genotypes
from snp_optimiser import SnpOptimiser
import pandas as pd
from os import listdir

import pickle

class TestIdentifySpeciesSnps(unittest.TestCase):
    valid_data=expanduser("~/HandyAmpliconTool/unit_test_data/valid_data/")
    invalid_data=expanduser("~/HandyAmpliconTool/unit_test_data/invalid_data/")
    config_file=expanduser("~/HandyAmpliconTool/unit_test_data/unittest.json")

    def setUp(self) -> None:
        self.config_data = InputConfiguration(self.config_file)
        self.config_data.specificity_limit=0.5
        self.config_data.sensitivity_limit=0.5
        self.snp_identifier=IdentifySpeciesSnps(ref_fasta=self.config_data.reference_fasta,
                                msa_dir=self.config_data.msa_dir,
                                negative_genomes_dir=self.config_data.negative_genomes,
                                temp_blast_db_dir=self.config_data.temp_blast_db,
                                amplicons_bed=self.config_data.multi_gt_intervals)
        
        if not exists(self.config_data.msa_dir):
            mkdir(self.config_data.msa_dir)
        
        # with open(self.config_data.genotypes_data, "rb") as pickled_file:
        #     gt_snps: Genotypes = pickle.load(pickled_file)

        # snp_opimiser=SnpOptimiser()
        # max_iterval_len=1000
        # amplicon_intervals=snp_opimiser.optimise(max_iterval_len,gt_snps, rare_gts=self.config_data.gts_with_few_snps)
        # with open(self.config_data.multi_gt_intervals, "w") as output_bed:
        #     for interval in amplicon_intervals:
        #         interval_start=min([f[1] for f in interval["snps"]])
        #         interval_end=max([f[1] for f in interval["snps"]])
        #         interval_len=interval_end-interval_start
        #         contig_id=interval["snps"][0][0]
        #         gts="_".join(interval["genotypes"])
        #         output_bed.write('\t'.join([contig_id, str(interval_start),str(interval_end),gts])+"\n")

        return super().setUp()
    

    def test_generate_flanking_amplicons(self):
        snp_identifier=IdentifySpeciesSnps.from_config(self.config_data)
        species_genotype: Genotype=snp_identifier.generate_flanking_amplicons()
        for amplicon in species_genotype.amplicons:
            if amplicon.has_flanking:
                has_left, has_right=[False, False]                
                #check that flanking regions are present and have correct coordinates
                for other_amplicon in species_genotype.amplicons:
                    if amplicon.left_flanking_id==other_amplicon.id:
                        self.assertEqual(amplicon.ref_start-self.config_data.flank_len_to_check,other_amplicon.ref_start)
                        has_left=True
                    elif amplicon.right_flanking_id==other_amplicon.id:
                        self.assertEqual(amplicon.ref_start+self.config_data.flank_len_to_check,other_amplicon.ref_end)
                        has_right=True
                self.assertTrue(has_left and has_right)
    
    def test_from_config(self):
        #nothing to check, just verify that function completes without errors
        self.assertTrue(IdentifySpeciesSnps.from_config(self.config_data))
    
    def test_get_bifurcating_snps(self):
        snp_identifier: IdentifySpeciesSnps=IdentifySpeciesSnps.from_config(self.config_data)
        species_genotype: Genotype=snp_identifier.generate_flanking_amplicons()
        self.assertEqual(species_genotype.name,"species")
        #self.assertTrue(len(species_genotype.amplicons)==6)
        species_genotype.amplicons=species_genotype.amplicons[0:3]
        species_genotype: Genotype=snp_identifier.get_bifurcating_snps(species_genotype)
        self.assertTrue( len(species_genotype.defining_snps)==45)
        for snp in species_genotype.defining_snp_coordinates:
            self.assertTrue( True in [f.coord_in_amplicon(snp) for f in species_genotype.amplicons] )

    def test_msa_df_to_msa_file(self):
        snp_identifier: IdentifySpeciesSnps=IdentifySpeciesSnps.from_config(self.config_data)
        dummy_df:pd.DataFrame=pd.DataFrame( data = {'col1': [1, 2], 'col2': [3, 4]} )
        self.assertTrue(  snp_identifier.msa_df_to_msa_file( dummy_df , "test") )

    def test_map_msa_to_ref_coordinates(self):
        #Already tested via test_get_bifurcating_snps
        pass



if __name__ == '__main__':
    unittest.main(verbosity=2)