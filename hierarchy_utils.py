#### At the moment this is placeholder. In future, this will be expanded to accomodate parsing of phylo trees to determine genotypes hierarchy
from typing import Dict, List, Tuple, Set
import pandas as pd
import metadata_utils  as mu
import name_converters as nc
from collections import Counter
from multiprocessing import cpu_count, Pool
import numpy as np
import warnings

class HierarchyUtilities:

    def __init__(self, sensitivity_limit: float, specificity_limit: float) -> None:
        self.genotype_hierarchy: Dict[str,List[str]]={}
        self.sensitivity_limit=sensitivity_limit
        self.specificity_limit=specificity_limit


    def load_hierarchy(self, filename) -> None:
        self.genotype_hierarchy={}
        with open(filename) as input_file:
            for i, line in enumerate(input_file):
                values=line.strip().split("\t")
                if values[0]=="":
                    raise ValueError(f'Line {i} in hierarchy file has empty first column: {line}')
                if values[0] in self.genotype_hierarchy:
                    raise ValueError(f'Genotype {values[0]} appears more than once in first column')
                self.genotype_hierarchy[ values[0] ]=[values[0]]
                for value in values[1:]:
                    if value in self.genotype_hierarchy[values[0]]:
                        raise ValueError(f'Genotype {value} appears more than once on line {i}: {line}')
                    self.genotype_hierarchy[values[0]].append(value)
       

    def get_subgenotypes(self, genotype) -> List[str]:
        if genotype in self.genotype_hierarchy:
            return self.genotype_hierarchy[genotype]
        raise ValueError(f'Genotype {genotype} not found in hierarchy')


    _snp_data: pd.DataFrame=pd.DataFrame()
    _column_to_gt: List[str]
    _genotype_snps: pd.DataFrame
    
    def _count_gt_allele_freq(self, target_gt: str ) -> Tuple[str, pd.DataFrame]:
        print(target_gt)
        target_columns=[f for f in self._snp_data.columns if mu.meta_data.loc[nc.get_sample(f),mu.genotype_column] in self.genotype_hierarchy[target_gt]]
        non_target_columns=[f for f in self._snp_data.columns if f not in target_columns]
        most_common_alleles=pd.DataFrame(index=self._snp_data.index, columns=["Target_top_allele","Target_top_count", "all_non_target_alleles",
                                                                              "non_Target_top_count","Specificity","Sensitivity","Pass"])
        most_common_alleles[ ["Target_top_allele","Target_top_count"] ]=[Counter(f).most_common()[0] for f in self._snp_data[ target_columns ].values]
        most_common_alleles[ "all_non_target_alleles" ]=[Counter(f) for f in self._snp_data[ non_target_columns ].values]
        most_common_alleles[ "non_Target_top_count" ]=[0 if most_common_target not in non_target_alleles else non_target_alleles[most_common_target] for non_target_alleles, most_common_target in most_common_alleles[ ["all_non_target_alleles","Target_top_allele"] ].values]
        total_target_samples=len(target_columns)
        total_non_target_samples=len(non_target_columns)
        most_common_alleles["Sensitivity"]=most_common_alleles["Target_top_count"]/total_target_samples
        most_common_alleles["Specificity"]=most_common_alleles["non_Target_top_count"]/total_non_target_samples
        most_common_alleles["Pass"]=[ sen>=self.sensitivity_limit and spe<=1-self.specificity_limit for  spe, sen in zip(most_common_alleles["Specificity"],most_common_alleles["Sensitivity"])]
        return (target_gt, most_common_alleles)


    def find_defining_snps(self, snp_data: pd.DataFrame) -> Dict[str, pd.DataFrame]: 
        self._column_to_gt: List[str]=[mu.meta_data.loc[nc.get_sample(f), mu.genotype_column] for f in snp_data.columns]
        self._snp_data=snp_data

        if __name__ == 'hierarchy_utils':
            present_gts=set()
            for gt in self.genotype_hierarchy.keys():
                if gt not in self._column_to_gt:
                    if len([ f for f in  self.genotype_hierarchy[gt] if f in self._column_to_gt ]):
                        warnings.warn(f'Genotype {gt} is not present in samples, but some of its subgenotypes are')
                        present_gts.add(gt)
                    else:
                        warnings.warn(f'Neither genotype {gt}, nor any of its subgenotypes are present among samples')
                else:
                    present_gts.add(gt)
        

            pool = Pool(processes= min(  max(cpu_count()-1,1) , len(present_gts) ) )
            #pool = Pool(processes= 1 )
            results=pool.map( self._count_gt_allele_freq, list(present_gts) ) 
            pool.close()
            pool.join()
            genotype_alleles:Dict[str, pd.DataFrame]={}
            for result in  results:
                genotype_alleles[result[0]]=result[1]
            
            return genotype_alleles
