#### At the moment this is placeholder. In future, this will be expanded to accomodate parsing of phylo trees to determine genotypes hierarchy
from typing import Dict, List, Tuple
import pandas as pd
import metadata_utils  as mu
import name_converters as nc
from collections import Counter
from multiprocessing import cpu_count, Pool, Lock

_lock=Lock()
class HierarchyUtilities:

    def __init__(self) -> None:
        self.genotype_hierarchy: Dict[str,List[str]]={}


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
    _genotype_snps: Dict[str, List[object]]={}

    def _find_defining_snps(self, unique_gt: str ) -> Tuple:
        print(unique_gt)
        gt_and_subgts=self.genotype_hierarchy[unique_gt]
        target_columns=[self._snp_data.columns[i] for i, f in enumerate(self._column_to_gt) if f in gt_and_subgts]
        non_target_columns=[self._snp_data.columns[i] for i, f in enumerate(self._column_to_gt) if f not in gt_and_subgts]
        unique_target=self._snp_data[target_columns].apply(lambda x: list(Counter(x).keys()), axis=1)
        unique_target.name="target_alleles"
        unique_non_target=self._snp_data[non_target_columns].apply(lambda x: list(Counter(x).keys()), axis=1)
        unique_non_target.name="non_target_alleles"
        merged_alleles=pd.merge(unique_target, unique_non_target, right_index=True, left_index=True)
        #if lenght of merge set of alleles is the same as lenght of both target and non-target alleles, there are no shared allels between two groups
        merged_alleles["segragating_allele"]=merged_alleles.apply(lambda x: len( set( x["target_alleles"]+x["non_target_alleles"] ) ) == len( x["target_alleles"])+len(x["non_target_alleles"]) , axis=1)
        return (unique_gt, list(merged_alleles.index[merged_alleles["segragating_allele"]]) )
        #_lock.aqcuire()
        #_lock.release()
        #if len(self._genotype_snps)>5:

    def find_defining_snps(self, snp_data: pd.DataFrame) -> Dict[str, List[object]]:
        self._column_to_gt: List[str]=[mu.meta_data.loc[nc.get_sample(f), mu.genotype_column] for f in snp_data.columns]
        self._snp_data=snp_data
        self._genotype_snps: Dict[str, List[object]]={}

        if __name__ == 'hierarchy_utils':
            pool = Pool(processes=max(cpu_count()-1,1))
            results=pool.map( self._find_defining_snps, set(self._column_to_gt))
            pool.close()
            pool.join()
        self._genotype_snps=dict(results)
        return self._genotype_snps
