from typing import List
from data_classes import  Genotypes, SNP, Genotype
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
        sorted_snps: List[SNP]=genotypes.all_snps_coord_sorted()
        i=0
        max_reached_index=0
        while i<len(sorted_snps):
            j=i #this means the first snp is automatically added
            while j<len(sorted_snps) and sorted_snps[j].ref_contig_id==sorted_snps[i].ref_contig_id and \
                  sorted_snps[j].position-sorted_snps[i].position<snp_interval:
                j+=1
            if j>max_reached_index: #this account for some intervals containing >2 snps. Without this, these intervals would create multiple entries in interval_snps
                interval_snps.append( {"snps": sorted_snps[i:j], "genotypes":[] } ) # j, not j-1 because python excludes the last element of index
                max_reached_index=j-1
            i+=1
        # check which GTs are captured by which lists of SNPs
        # Remove those that capture same GT multiple times - this is likely due to structural variant
        bifurcating_gts=set() #these are GTs that divide dataset into two parts, both of which are targets
        #this means same SNP will capture both of these genotypes. This makes every such SNP 
        #a multi GT SNP, which is not how it is supposed to be. 
        for interval in interval_snps:
            if len(interval["snps"])>9: #allow for some SNPs in close proximity, but limit to 9 SNPs per interval
                continue
            interval["genotypes"]=set()
            for snp in interval["snps"]:
                gts_with_snp: List[Genotype]=genotypes.genotypes_with_snp(snp)
                if len(gts_with_snp)>1:
                    #possibly genotyping scheme consists of two genotypes splitting all samples.
                    #in this case use ALT defined genotypes
                    gt_with_alt=[f for f in gts_with_snp if f.get_genotype_allele(snp)==snp.alt_base][0]
                    interval["genotypes"].add(gt_with_alt.name)
                    gt_with_ref=[f for f in gts_with_snp if f.get_genotype_allele(snp)==snp.ref_base][0]
                    bifurcating_gts.add(gt_with_ref.name)
                else:
                    interval["genotypes"].add(gts_with_snp[0].name)
                
                    

        interval_snps=[snp_interval for snp_interval in interval_snps if len(snp_interval["genotypes"])>1 or len(set(snp_interval["genotypes"]) & set(rare_gts))>0 ]
        captured_genotypes=set([f for genotypes in interval_snps for f in genotypes["genotypes"]])
        for gt in bifurcating_gts: #the reference GT will not captured through about algorithms, so has to added manually.
            captured_genotypes.add(gt)
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

