import pandas as pd
import numpy as np
from collections import Counter
from typing import List
from data_classes import Genotype, Genotypes
class SnpOptimiser:

    def __init__(self) -> None:
        pass

    def optimise(self, snp_interval: int, genotypes: Genotypes, rare_gts: List[str]) -> List[object]:
        """Select SNPs based on how many fall within a maximum permitted amplicon range

        :param snp_interval: Maximum length of amplicon 
        :type snp_interval: int

        :param genotypes: genotypes with already identified genotypes defining SNPs
        :type genotypes: Genotypes
        """        
        interval_snps=[]
        sorted_snps=genotypes.all_snps_coord_sorted()
        i=0
        max_reached_index=0
        while i<len(sorted_snps):
            j=i #this means the first snp is automatically added
            while j<len(sorted_snps) and sorted_snps[j][0]==sorted_snps[i][0] and sorted_snps[j][1]-sorted_snps[i][1]<snp_interval:
                j+=1
            if j>max_reached_index: #this account for some intervals containing >2 snps. Without this, these intervals would create multiple entries in interval_snps
                interval_snps.append( {"snps": sorted_snps[i:j], "genotypes":[] } ) # j, not j-1 because python excludes the last element of index
                max_reached_index=j-1
            i+=1
        # check which GTs are captured by which lists of SNPs
        # Remove those that capture same GT multiple times - this is likely due to structural variant
        for interval in interval_snps:
            interval["genotypes"]=set()
            for genotype in genotypes.genotypes:
                temp=set(interval['snps']) & set(genotype.defining_snp_coordinates) #this is probably not the most efficient way.
                if len(temp)>9: #allow for some SNPs in close proximity. 9 is three AA deletion or insertion
                    continue
                elif len(temp)<=9 and len(temp)>0: #SNPs or upto 3 codons as variants
                    interval["genotypes"].add(genotype.name)

        interval_snps=[snp_interval for snp_interval in interval_snps if len(snp_interval["genotypes"])>1 or len(set(snp_interval["genotypes"]) & set(rare_gts))>0 ]
        captured_genotypes=set([f for genotypes in interval_snps for f in genotypes["genotypes"]])
        not_captured_genotypes=[f.name for f in genotypes.genotypes if f.name not in captured_genotypes ]
        if len(not_captured_genotypes)==0:
            print("Genotypes not captured by multigenotype intervals: None")
        else:
            print(f'Genotypes not captured by multigenotype intervals: {", ".join(not_captured_genotypes)}')
        not_captured_due_to_no_snps=[f.name for f in genotypes.genotypes if len(f.defining_snp_coordinates)==0 ]
        if len(not_captured_due_to_no_snps)>0:
            print(f'Of these {", ".join(not_captured_due_to_no_snps)} do not have defining SNPs')
        else:
            print(f'All of them have defining SNPs')
        return interval_snps

