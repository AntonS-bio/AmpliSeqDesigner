#### At the moment this is placeholder. In future, this will be expanded to accomodate parsing of phylo trees to determine genotypes hierarchy
from typing import Dict, List
import pandas as pd
import metadata_utils  as mu
import name_converters as nc
from collections import Counter

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

    def find_defining_snps(self, snp_data: pd.DataFrame) -> Dict[str, List[object]]:
        column_to_gt: List[str]=[mu.meta_data.loc[nc.get_sample(f), mu.genotype_column] for f in snp_data.columns]
        genotype_snps: Dict[str, List[object]]={}
        for unique_gt in set(column_to_gt):
            print(unique_gt)
            gt_and_subgts=self.genotype_hierarchy[unique_gt]
            target_columns=[snp_data.columns[i] for i, f in enumerate(column_to_gt) if f in gt_and_subgts]
            non_target_columns=[snp_data.columns[i] for i, f in enumerate(column_to_gt) if f not in gt_and_subgts]
            unique_target=snp_data[target_columns].apply(lambda x: list(Counter(x).keys()), axis=1)
            unique_target.name="target_alleles"
            unique_non_target=snp_data[non_target_columns].apply(lambda x: list(Counter(x).keys()), axis=1)
            unique_non_target.name="non_target_alleles"
            merged_alleles=pd.merge(unique_target, unique_non_target, right_index=True, left_index=True)
            #if lenght of merge set of alleles is the same as lenght of both target and non-target alleles, there are no shared allels between two groups
            merged_alleles["segragating_allele"]=merged_alleles.apply(lambda x: len( set( x["target_alleles"]+x["non_target_alleles"] ) ) == len( x["target_alleles"])+len(x["non_target_alleles"]) , axis=1)
            genotype_snps[unique_gt]=list(merged_alleles.index[merged_alleles["segragating_allele"]])
        return genotype_snps
