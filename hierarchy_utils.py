#### At the moment this is placeholder. In future, this will be expanded to accomodate parsing of phylo trees to determine genotypes hierarchy
from typing import Dict, List, Tuple, Set
import pandas as pd
import metadata_utils  as mu
import name_converters as nc
from collections import Counter
from multiprocessing import cpu_count, Pool
import numpy as np

class HierarchyUtilities:

    def __init__(self) -> None:
        self.genotype_hierarchy: Dict[str,List[str]]={}
        self.sensitivity_limit=0.98
        self.specificity_limit=0.98


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
    # def _find_defining_snps(self, unique_gt: str ) -> Tuple:
    #     print(unique_gt)
    #     gt_and_subgts=self.genotype_hierarchy[unique_gt]
    #     target_columns=[self._snp_data.columns[i] for i, f in enumerate(self._column_to_gt) if f in gt_and_subgts]
    #     non_target_columns=[self._snp_data.columns[i] for i, f in enumerate(self._column_to_gt) if f not in gt_and_subgts]
    #     unique_target=self._snp_data[target_columns].apply(lambda x: list(Counter(x).keys()), axis=1)
    #     unique_target.name="target_alleles"
    #     unique_non_target=self._snp_data[non_target_columns].apply(lambda x: list(Counter(x).keys()), axis=1)
    #     unique_non_target.name="non_target_alleles"
    #     merged_alleles=pd.merge(unique_target, unique_non_target, right_index=True, left_index=True)
    #     #if lenght of merge set of alleles is the same as lenght of both target and non-target alleles, there are no shared allels between two groups
    #     merged_alleles["segragating_allele"]=merged_alleles.apply(lambda x: len( set( x["target_alleles"]+x["non_target_alleles"] ) ) == len( x["target_alleles"])+len(x["non_target_alleles"]) , axis=1)
    #     return (unique_gt, list(merged_alleles.index[merged_alleles["segragating_allele"]]) )
    
    def _count_gt_allele_freq(self, unique_gt: str ) -> Tuple:
        print(unique_gt)
        columns=self._snp_data.columns[np.asarray(self._column_to_gt)==unique_gt]
        return (unique_gt, [Counter(f) for f in self._snp_data[columns].to_numpy()])

    def _count_multi_gt_allele_freq(self, gts_to_exclude: Set[str] ) -> Tuple:
        columns=self._snp_data.columns[ [f not in gts_to_exclude for f in self._column_to_gt] ]
        gts_to_use=set([f for f in self._column_to_gt if f not in gts_to_exclude])
        #columns=self._snp_data.columns[np.asarray(self._column_to_gt)==gts_to_use]
        return (gts_to_use, [Counter(f) for f in self._snp_data[columns].to_numpy()])        
    
    def _find_defining_snps(self,target_gt: str) -> Tuple[str,pd.DataFrame]:
        target_gt_and_subgts=set([f for f in  self._column_to_gt if f in self.genotype_hierarchy[target_gt]])
        #there are three sets of gts:
        #those that are targeted in this specific call
        #those that are never targeted
        #those that are targeted by another call
        #the last category needs to be carefully managed
        non_target_gts=[f for f in self._genotype_snps.columns if f not in target_gt_and_subgts ]
        allele_split=pd.DataFrame(index=self._genotype_snps.index, columns=["GT", "Pass","Target_Allele", "Count_target","Count_non_target"])
        for index in self._genotype_snps.index:
            target_df=pd.DataFrame([(key,item) for f in self._genotype_snps.loc[ [index], list(target_gt_and_subgts)].values[0] for key, item in f.items()], columns=["Allele","Count"])
            target_alleles=pd.pivot_table(data=target_df, values="Count", index="Allele", aggfunc=sum)
            most_frequent_target_allele=target_alleles.index[np.argmax(target_alleles)]
            most_frequent_count_target=target_alleles.loc[most_frequent_target_allele,"Count"]

            non_target_df=pd.DataFrame([(key,item) for f in self._genotype_snps.loc[ [index], list(non_target_gts)].values[0] for key, item in f.items()], columns=["Allele","Count"])
            non_target_alleles=pd.pivot_table(data=non_target_df, values="Count", index="Allele", aggfunc=sum)
            if most_frequent_target_allele in non_target_alleles.index:
                most_frequent_count_non_target=non_target_alleles.loc[most_frequent_target_allele,"Count"]
            else:
                most_frequent_count_non_target=0

            #sensitivity - most frequent alleles prevalence in >X%
            if most_frequent_count_target >= target_alleles["Count"].sum()*self.sensitivity_limit and \
                most_frequent_count_non_target <= non_target_alleles["Count"].sum()*(1-self.specificity_limit):
                allele_split.loc[ [index], ["GT", "Pass","Target_Allele", "Count_target","Count_non_target"]]=[target_gt,True,most_frequent_target_allele, most_frequent_count_target, most_frequent_count_non_target]
            else:
                allele_split.loc[ [index],  ["GT", "Pass","Target_Allele", "Count_target","Count_non_target"]]=[target_gt,False,most_frequent_target_allele, most_frequent_count_target, most_frequent_count_non_target]
        return (target_gt, allele_split)


    def find_defining_snps(self, snp_data: pd.DataFrame) -> Dict[str, pd.DataFrame]: 
        self._column_to_gt: List[str]=[mu.meta_data.loc[nc.get_sample(f), mu.genotype_column] for f in snp_data.columns]
        self._snp_data=snp_data
        self._genotype_snps: pd.DataFrame

        if __name__ == 'hierarchy_utils':
            gts_to_use=set()
            for target_gt, target_subgts in self.genotype_hierarchy.items():
                #add subgenotypes to set of target gts, but exclude the gts no in dataset samples
                gts_to_use.update( [f for f in target_subgts if f in self._column_to_gt]  )
                gts_to_use.add(target_gt)

            #gts_to_use=list(set([f for f in  self._column_to_gt if f in self.genotype_hierarchy]))
            pool = Pool(processes= min(max(cpu_count()-1,1), len(gts_to_use) ))
            result=pool.map( self._count_gt_allele_freq, list(gts_to_use) ) # dict[key=gt, value=Counter(alleleFreq)]
            pool.close()
            pool.join()
            gts_alles: Dict[str, Counter]= dict(result)
            del result
            non_target_gts=self._count_multi_gt_allele_freq(gts_to_use)
            
        
            gts_alleles_df=pd.DataFrame(index=snp_data.index, columns=list(gts_alles.keys()) )
            gts_alleles_df["Non_target_gts"]=non_target_gts[1]
            for gt in gts_alles:
                gts_alleles_df[gt]=gts_alles[gt]

            self._genotype_snps=gts_alleles_df
            #pool = Pool(processes= min(max(cpu_count()-1,1), len(self.genotype_hierarchy) ))
            pool = Pool(processes=1)
            results=pool.map( self._find_defining_snps, list(self.genotype_hierarchy.keys()) ) # dict[key=gt, value=DataFrame]
            pool.close()
            pool.join()
            genotype_alleles:Dict[str, pd.DataFrame]={}
            for result in  results:
                genotype_alleles[result[0]]=result[1]
            
            return genotype_alleles
