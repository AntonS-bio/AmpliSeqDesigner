from typing import Dict, List
import pandas as pd
from collections import Counter
from data_classes import Genotype, Genotypes, Sample, InputConfiguration

class HierarchyUtilities:
    """Class representing a hierarchy structure of the target organims.
    The data is supplied as tab-delimited file
    each row must contain (in first column) the target genotype
    and MAY also contain in subsequent column subgenotypes of target genotype
    There is no need to create rows for non-target genotypes
    If the genotype has not subgenotypes, the row would only have target genotype
    ex: 
    4.3.1 4.3.1.1.P1 4.3.1.2 4.3.1.2
    3.1.1
    """
    def __init__(self) -> None:
        self.genotype_hierarchy: Dict[str,Genotype]={}


    def load_hierarchy(self, filename: str) -> Dict[str,Genotype]:
        """Loads hierarchy data
        :param filename: Path to tab delimited file containing genotype hierarchy
        :type filename: str
        """
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
        return self.genotype_hierarchy

    _snp_data: pd.DataFrame=pd.DataFrame()
    _column_to_gt: List[str]
    _genotype_snps: pd.DataFrame

    def find_defining_snps(self, samples: List[Sample]) -> Genotypes:
        """Identifies SNPs that are specific to genotypes
        :param samples: Collection of all samples that were loaded from VCF files.
        :type samples: List[Sample]
        """
        genotypes=Genotypes()
        for gt_name, genotype in self.genotype_hierarchy.items():
            gt_samples=[f for f in samples if f.genotype in genotype.subgenotypes]
            non_gt_samples=[f for f in samples if f.genotype not in genotype.subgenotypes]
            gt_snps=Counter([snp for sample in gt_samples for snp in sample.snps])
            non_gt_snps=Counter([snp for sample in non_gt_samples for snp in sample.snps])
            for snps_dict, invert_specificity_sensitivity in zip([gt_snps, non_gt_snps], [False, True]):
                #a genotype can be defined by SNPs not present in it or SNPs present in it
                #for this reason, SNPs in genotype samples are not sufficient and SNPs not in genotype also have to be checked
                #this creates potiential double counting, which needs to be checked
                for gt_snp, gt_snp_count in snps_dict.items():
                    if invert_specificity_sensitivity:
                        specificity=1-gt_snp_count/len(gt_samples)
                        sensitivity=non_gt_snps[gt_snp]/len(non_gt_samples)
                        allele_depth=non_gt_snps[gt_snp]
                    else:
                        sensitivity=gt_snp_count/len(gt_samples)
                        specificity=1-non_gt_snps[gt_snp]/len(non_gt_samples)
                        allele_depth=gt_snp_count 
                    if sensitivity>InputConfiguration.specificity_limit and specificity>InputConfiguration.sensitivity_limit and gt_snp not in genotype.defining_snps:
                        # gt_snp not in genotype.defining_snps check for redundancy
                        snp_copy=gt_snp.copy()
                        snp_copy.sensitivity=sensitivity
                        snp_copy.specificity=specificity
                        snp_copy.passes_filters=True
                        snp_copy.is_genotype_snp=True
                        genotype_allele=gt_snp.ref_base if invert_specificity_sensitivity else gt_snp.alt_base
                        genotype.add_genotype_allele(snp_copy, genotype_allele, allele_depth )

            genotypes.genotypes.append(genotype)
            print(f'{genotype.name} has {str(len(genotype.defining_snps))} SNPs')
        return genotypes

       