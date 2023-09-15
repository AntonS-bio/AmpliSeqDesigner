#### At the moment this is placeholder. In future, this will be expanded to accomodate parsing of phylo trees to determine genotypes hierarchy
from typing import Dict, List, Tuple, Set
import pandas as pd
import metadata_utils  as mu
import name_converters as nc
from collections import Counter
from multiprocessing import cpu_count, Pool
import numpy as np
import warnings
from data_classes import Genotype, Genotypes, SNP, Sample
from tqdm import tqdm

class HierarchyUtilities:

    def __init__(self, sensitivity_limit: float, specificity_limit: float) -> None:
        self.genotype_hierarchy: Dict[str,Genotype]={}
        self.sensitivity_limit=sensitivity_limit
        self.specificity_limit=specificity_limit


    def load_hierarchy(self, filename) -> None:
        self.genotype_hierarchy: Dict[str,Genotype]={}
        with open(filename) as input_file:
            for i, line in enumerate(input_file):
                values=line.strip().split("\t")
                if values[0]=="":
                    raise ValueError(f'Line {i} in hierarchy file has empty first column: {line}')
                if values[0] in self.genotype_hierarchy:
                    raise ValueError(f'Genotype {values[0]} appears more than once in first column')
                genotype=Genotype(values[0])
                for value in values[1:]:
                    if value in genotype.subgenotypes:
                        raise ValueError(f'Genotype {value} appears more than once on line {i}: {line}')
                    genotype.subgenotypes.append(value)
                self.genotype_hierarchy[ values[0] ]=genotype

    _snp_data: pd.DataFrame=pd.DataFrame()
    _column_to_gt: List[str]
    _genotype_snps: pd.DataFrame
    
   
    def _count_gt_allele_freq(self, target_gt: Genotype ) -> Genotype:
        target_columns=[f for f in self._snp_data.columns if mu.meta_data.loc[nc.get_sample(f),mu.genotype_column] in target_gt.subgenotypes]
        non_target_columns=[f for f in self._snp_data.columns if f not in target_columns]

        total_target_samples=len(target_columns)
        total_non_target_samples=len(non_target_columns)
        for target_col_alleles, non_target_col_alleles, position in zip([Counter(f).most_common()[0] for f in self._snp_data[ target_columns ].values],
                                                                  [Counter(f) for f in self._snp_data[ non_target_columns ].values],
                                                                  self._snp_data.index):
            snp=SNP(ref_contig_id=position[0], position=position[1],
                    ref_base="REF", alt_base=target_col_alleles[0], )
            total_target_in_non_target=0 if target_col_alleles[0] not in non_target_col_alleles else non_target_col_alleles[target_col_alleles[0]]
            snp.sensitivity=target_col_alleles[1]/total_target_samples
            snp.specificity=1-total_target_in_non_target/total_non_target_samples
            snp.passes_filters = snp.sensitivity>self.specificity_limit and snp.specificity>self.sensitivity_limit
            if snp.passes_filters:
                target_gt.defining_snps.append(snp)

        return target_gt


    def find_defining_snps(self, samples: List[Sample]) -> Genotypes:
        genotypes=Genotypes()
        for gt_name, genotype in self.genotype_hierarchy.items():
            gt_samples=[f for f in samples if f.genotype in genotype.subgenotypes]
            non_gt_samples=[f for f in samples if f.genotype not in genotype.subgenotypes]
            gt_snps=Counter([snp for sample in gt_samples for snp in sample.snps])
            non_gt_snps=Counter([snp for sample in non_gt_samples for snp in sample.snps])
            for gt_snp, gt_snp_count in gt_snps.items():
                gt_snp.sensitivity=gt_snp_count/len(gt_samples)
                gt_snp.specificity=1-non_gt_snps[gt_snp]/len(non_gt_samples)
                gt_snp.passes_filters = gt_snp.sensitivity>self.specificity_limit and gt_snp.specificity>self.sensitivity_limit
                gt_snp.is_genotype_snp=gt_snp.passes_filters
            genotype.defining_snps=[gt for  gt in gt_snps.keys() if gt.passes_filters]
            genotypes.genotypes.append(genotype)
            print(len(genotype.defining_snps))
        return genotypes

    # def find_defining_snps(self, snp_data: pd.DataFrame) -> Genotypes:
    #     self._column_to_gt: List[str]=[mu.meta_data.loc[nc.get_sample(f), mu.genotype_column] for f in snp_data.columns]
    #     self._snp_data=snp_data

    #     if __name__ == 'hierarchy_utils':
    #         present_gts: Set[Genotype]=set()
    #         for gt_name, genotype in self.genotype_hierarchy.items():
    #             if gt_name not in self._column_to_gt:
    #                 if len([ f for f in  genotype.subgenotypes if f in self._column_to_gt ]):
    #                     warnings.warn(f'Genotype {gt_name} is not present in samples, but some of its subgenotypes are')
    #                     present_gts.add(genotype)
    #                 else:
    #                     warnings.warn(f'Neither genotype {gt_name}, nor any of its subgenotypes are present among samples')
    #             else:
    #                 present_gts.add(genotype)
        
    #         print("Identifying genotype SNPs")
    #         pool = Pool(processes= min(  max(cpu_count()-1,1) , len(present_gts) ) )
    #         results = list(tqdm( pool.imap(func=self._count_gt_allele_freq, iterable=list(present_gts)), total=len(present_gts) ))
    #         return Genotypes(genotypes=results)
    #         #pool = Pool(processes= 1 )
    #         results=pool.map( self._count_gt_allele_freq, list(present_gts) ) 
    #         pool.close()
    #         pool.join()
    #         return results



