import primer3
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
import warnings
from os.path import exists
from matplotlib import pyplot as plt
from typing import Dict, List
root_dir="/home/lshas17/"
import pandas as pd
from itertools import product
from sys import path
path.insert(1, f'{root_dir}/HandyAmpliconTool/scripts/')
from data_classes import BlastResult, PrimerPair, Primer, InputConfiguration


class PrimersGenerator():
    def __init__(self, config: InputConfiguration) -> None:
        self.ref_seq: Dict[str, str]={}
        self.config=config
        if not exists(config.reference_fasta):
            raise IOError(f'Reference FASTA {config.reference_fasta} does not exist')
        for record in SeqIO.parse(config.reference_fasta, "fasta"):
            self.ref_seq[record.id]=str(record.seq)
        self.existing_primers: List[Primer]=[]


    def _get_seq_coordinates_in_ref(self, sequence) -> int:
        seq=Seq(sequence)
        seq_rc=str(seq.reverse_complement())
        for ref_contig in self.ref_seq.values():
            direct_start=ref_contig.find(str(seq))
            if direct_start!=-1:
                return direct_start
            else:
                rc_start=ref_contig.find(seq_rc)
                if rc_start!=-1:
                    return rc_start
        return -1

    def _sequence_to_primer(self, primer_seq) -> Primer:
        primer=Primer(primer_seq, gc_fraction(primer_seq), primer3.calc_tm(primer_seq) )
        primer.ref_start=self._get_seq_coordinates_in_ref(primer_seq)
        return primer

    def load_existing_primer_from_fasta(self, existing_primers_fasta: str):
        self.existing_primers.clear()
        for record in SeqIO.parse(existing_primers_fasta, "fasta"):
            primer=self._sequence_to_primer( str(record.seq) )
            self.existing_primers.append(primer)

    def load_existing_primer_from_bed(self, existing_primers_bed: str):
        self.existing_primers.clear()
        with open(existing_primers_bed) as bed_file:
            for i, line in enumerate(bed_file):
                values=line.strip().split("\t")
                if values[0] not in self.ref_seq:
                    raise ValueError(f'Contig {values[0]} from existing amplicons file {existing_primers_bed} not found in reference fasta {self.config.reference_fasta}')
                if int(values[1])+1>len(self.ref_seq[values[0]]):
                    raise ValueError(f'Amplicon {i+1} in existing amplicons file {existing_primers_bed} has coordinate {values[0]} {int(values[1])} greater than lenght of reference sequence: {self.config.reference_fasta}')
                primer_seq=self.ref_seq[values[0]][int(values[1]):int(values[1])+1]
                primer=self._sequence_to_primer( primer_seq )
                self.existing_primers.append(primer)

    # vcf_file_name=f'{root_dir}/HandyAmpliconTool/test_data/outputs/2_3_1/snps.vcf'
    # header_rows=0
    # for line in open(vcf_file_name):
    #     if line[0]=="#":
    #         header_rows+=1
    # vcf_data=pd.read_csv(vcf_file_name, sep="\t", header=0, skiprows=header_rows-1) # -1 because header row also starts with #CHROM



    # def find_primer_candidate(self):
    #     target_lineage="2.3.1"
    #     interval_len=1500
    #     min_interval_len=200
    #     target_lineage_index=vcf_data.index[ [f for f in vcf_data.index if vcf_data.loc[f, "ID"].find(target_lineage)>-1]  ]

    #     new_primer_pairs: List[PrimerPair] = []

    #     for index in target_lineage_index:
    #         snp_position=vcf_data.loc[index, "POS"]
    #         valid_range=range(snp_position-interval_len, snp_position+interval_len)
    #         interval_indices=[ f for f in vcf_data.index if vcf_data.loc[f, "POS"] in valid_range ]

    #         interval_left_serovar_snps=[f for f in interval_indices if vcf_data.loc[f,"ID"].find("Serovar")>-1 and f<index ]
    #         interval_right_serovar_snps=[f for f in interval_indices if vcf_data.loc[f,"ID"].find("Serovar")>-1 and f>index]
    #         if len(interval_left_serovar_snps)>0 and len(interval_right_serovar_snps)>0:
    #             print(f'SNP at {snp_position} has serovar SNPs left and right of genotype SNP')
    #         elif len(interval_right_serovar_snps)>0:
    #             print(f'SNP at {snp_position} has serovar SNPs right of genotype SNP')
    #         elif len(interval_left_serovar_snps):
    #             print(f'SNP at {snp_position} has serovar SNPs left of genotype SNP')
    #         else:
    #             print(f'SNP at {snp_position} has no serovar SNPs within {str(interval_len)} region')


    #         #four three cases: 
    #         # 1 - left and right serovar SNPs anchored primers
    #         # 2 - left serovar anchored primer
    #         # 3 - right serovar anchored primer
    #         # 4 - neither side has serovar primer, ignore this kind of SNP

    #         for left_snp_index, right_snp_index in product(interval_left_serovar_snps, interval_right_serovar_snps):
    #             print("fixed")
    #             #case 1
    #             right_snp_pos=vcf_data.loc[right_snp_index,"POS"]
    #             left_snp_pos=vcf_data.loc[left_snp_index,"POS"]
    #             if right_snp_pos-left_snp_pos>interval_len or right_snp_pos-left_snp_pos<min_interval_len: #distance between SNPs is too long:
    #                 print(f'SNPs too far apart: {left_snp_pos} to {right_snp_pos} vs max len {str(interval_len)}')
    #                 continue
    #             forward=ref_seq[left_snp_pos-1:left_snp_pos+20]
    #             reverse=ref_seq[right_snp_pos-20:right_snp_pos]
    #             seq_args={"SEQUENCE_ID":f'{target_lineage}_{str(snp_position)}', 
    #                     "SEQUENCE_TEMPLATE": str(ref_seq[left_snp_pos-1:right_snp_pos].seq),
    #                     "SEQUENCE_INCLUDED_REGION": [0, right_snp_pos-left_snp_pos+1 ]}
    #             #add "fixed" to name to show that this pair has both sides fixed
    #             global_args=add_fixed_primer(for_seq=str(forward.seq),rev_seq=str(reverse.seq.reverse_complement()), template=seq_args["SEQUENCE_TEMPLATE"])
    #             primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
    #             new_primer_pairs += process_p3_output(primers, id_prefix=f'{target_lineage}_{snp_position}_both_fixed')

    #         for right_snp_index in interval_right_serovar_snps:
    #             print("right")
    #             right_snp_pos=vcf_data.loc[right_snp_index,"POS"]
    #             reverse=ref_seq[right_snp_pos-20:right_snp_pos]
    #             seq_args={"SEQUENCE_ID":f'{target_lineage}_{str(snp_position)}', 
    #                     "SEQUENCE_TEMPLATE": str(ref_seq[right_snp_pos-interval_len:right_snp_pos].seq),
    #                     "SEQUENCE_INCLUDED_REGION": [0, interval_len ]}
    #             global_args=add_fixed_primer(for_seq=None,rev_seq=str(reverse.seq.reverse_complement()),  template=seq_args["SEQUENCE_TEMPLATE"])
    #             primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
    #             new_primer_pairs += process_p3_output(primers, id_prefix=f'{target_lineage}_{snp_position}_right_fixed')

    #         for left_snp_index in interval_left_serovar_snps:
    #             print("left")
    #             left_snp_pos=vcf_data.loc[left_snp_index,"POS"]
    #             forward=ref_seq[left_snp_pos-1:left_snp_pos+20] #-1 due to VCF being 1 indexed
    #             seq_args={"SEQUENCE_ID":f'{target_lineage}_{str(snp_position)}', 
    #                     "SEQUENCE_TEMPLATE": str(ref_seq[left_snp_pos-20:left_snp_pos+interval_len].seq),
    #                     "SEQUENCE_INCLUDED_REGION": [0, interval_len ]}
    #             global_args=add_fixed_primer(for_seq=str(forward.seq), rev_seq=None, template=seq_args["SEQUENCE_TEMPLATE"])
    #             primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
    #             new_primer_pairs += process_p3_output(primers, id_prefix=f'{target_lineage}_{snp_position}_both_fixed')

    #         #sometimes the same pair of F/R sequences is generated multiple times because there are multiple species SNPs around same genotype SNP
    #         unique_new_primer_pairs: List[PrimerPair]=[]
    #         existing_sequences=set()
    #         for pair in new_primer_pairs:
    #             sequences=pair.forward.seq+"_"+pair.reverse.seq
    #             if sequences not in existing_sequences:
    #                 unique_new_primer_pairs.append(pair)
    #                 existing_sequences.add(sequences)



    # def add_fixed_primer(for_seq: str, rev_seq:str, template: str):
    #     global_args={
    #         'PRIMER_OPT_SIZE': 20,
    #         'PRIMER_OPT_TM': 60.0,
    #         'PRIMER_MIN_TM': 50.0,
    #         'PRIMER_MAX_TM': 70.0,
    #         'PRIMER_PRODUCT_SIZE_RANGE': f'100-{len(template)}'}
    #     if for_seq!=None:
    #         global_args['SEQUENCE_PRIMER']=for_seq
    #     if rev_seq!=None:
    #         global_args['SEQUENCE_PRIMER_REVCOMP']=rev_seq
    #     return global_args
    # def has_homodimers(forward_seq, reverse_seq) -> bool:
    #     thermoresult_forward=primer3.bindings.calc_homodimer(forward_seq)
    #     thermoresult_reverse=primer3.bindings.calc_homodimer(reverse_seq)
    #     if (thermoresult_forward.structure_found and thermoresult_forward.tm>10) or \
    #             (thermoresult_reverse.structure_found and thermoresult_reverse.tm>10):
    #         return True
    #     else:
    #         return False 
    # def check_heterodimers(forward_seq, reverse_seq):
    #     global existing_primers
    #     total_structures=0
    #     for primer in existing_primers:
    #         if primer.find("Forward")>-1:
    #             forward_result=primer3.bindings.calc_heterodimer(reverse_seq, existing_primers[primer])
    #             total_structures+=1 if forward_result.structure_found and forward_result.tm>10 else 0
    #         else:
    #             reverse_result=primer3.bindings.calc_heterodimer(forward_seq, existing_primers[primer])
    #             total_structures+=1 if reverse_result.structure_found and reverse_result.tm>10 else 0
    #     return total_structures
    # def process_p3_output(output: Dict, id_prefix: str) -> List[PrimerPair]:
    # primer_pairs=len(output["PRIMER_PAIR"])
    # result: List[PrimerPair]=[]
    # for i in range(0,primer_pairs):
    #     forward=Primer(seq=output[f'PRIMER_LEFT_{str(i)}_SEQUENCE'],
    #                    g_c=output[f'PRIMER_LEFT_{str(i)}_GC_PERCENT'],
    #                    t_m=output[f'PRIMER_LEFT_{str(i)}_TM'])
    #     reverse=Primer(seq=output[f'PRIMER_RIGHT_{str(i)}_SEQUENCE'],
    #                    g_c=output[f'PRIMER_RIGHT_{str(i)}_GC_PERCENT'],
    #                    t_m=output[f'PRIMER_RIGHT_{str(i)}_TM'])
    #     new_pair=PrimerPair(name=id_prefix, forward=forward, reverse=reverse)
    #     new_pair.penalty=output[f'PRIMER_PAIR_{str(i)}_PENALTY']
    #     if not has_homodimers(forward.seq, reverse.seq) and check_heterodimers(forward.seq, reverse.seq)==0:
    #         result.append(new_pair)
    #         #output[f'PRIMER_PAIR_{str(i)}_PRODUCT_SIZE']
    #         print(f'LEFT: {output["PRIMER_LEFT_EXPLAIN"]}')
    #         print(f'RIGHT: {output["PRIMER_RIGHT_EXPLAIN"]}')
    #         print(f'PAIR: {output["PRIMER_PAIR_EXPLAIN"]}')
    #         print("\n")
    # return result