import primer3
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
import warnings
from os.path import exists
from matplotlib import pyplot as plt
from typing import Dict, List, Tuple
root_dir="/home/lshas17/"
import pandas as pd
from itertools import product
from sys import path
path.insert(1, f'{root_dir}/HandyAmpliconTool/scripts/')
from data_classes import BlastResult, PrimerPair, Primer, InputConfiguration, Genotype, Genotypes, SNP
from inputs_validation import ValidateFiles
import pickle


class PrimersGenerator():
    def __init__(self, config: InputConfiguration) -> None:
        self.ref_seq: Dict[str, str]={}
        self.config=config
        if not exists(config.reference_fasta):
            raise IOError(f'Reference FASTA {config.reference_fasta} does not exist')
        for record in SeqIO.parse(config.reference_fasta, "fasta"):
            self.ref_seq[record.id]=str(record.seq)
        self.existing_primers: List[Primer]=[]
        self.new_primer_pairs: List[PrimerPair]=[]


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

    def _sequence_to_primer(self, primer_seq:str, is_reverse:bool) -> Primer:
        """Generates a Primer object from primer sequence.
        :param primer_seq: sequence of the primer
        :type primer_seq: str
        :param is_reverse: Indicates if primer sequence is for reverse primer
        :type is_reverse: bool        
        :return: primer object 
        :rtype: Primer
        """
        primer=Primer(primer_seq, gc_fraction(primer_seq), primer3.calc_tm(primer_seq), is_reverse )
        primer.ref_start=self._get_seq_coordinates_in_ref(primer_seq)
        return primer

    def _id_has_direction_info(self, id_value: str) -> Tuple[bool, str]:
        if str.lower(id_value).find("forward")>-1:
            return (True, "forward")
        elif str.lower(id_value).find("reverse")>-1:
            return (True, "reverse")
        else:
            return (False, "")
        
    def load_existing_primer_from_fasta(self, existing_primers_fasta: str):
        self.existing_primers.clear()
        for record in SeqIO.parse(existing_primers_fasta, "fasta"):
            direction_info=self._id_has_direction_info(record.id)
            if direction_info[0]:
                primer=self._sequence_to_primer( str(record.seq), is_reverse=direction_info[1]=="reverse" )
            else: #this is worst case, add primer as both forward and reverse
                primer=self._sequence_to_primer( str(record.seq), is_reverse=True )
                self.existing_primers.append(primer)
                primer=self._sequence_to_primer( str(record.seq), is_reverse=False)
            self.existing_primers.append(primer)

    def load_existing_primer_from_bed(self, existing_primers_bed: str):
        file_validator=ValidateFiles()
        file_validator.validate_bed(existing_primers_bed)
        self.existing_primers.clear()
        with open(existing_primers_bed) as bed_file:
            for i, line in enumerate(bed_file):
                values=line.strip().split("\t")
                if values[0] not in self.ref_seq:
                    raise ValueError(f'Contig {values[0]} from existing amplicons file {existing_primers_bed} not found in reference fasta {self.config.reference_fasta}')
                if int(values[1])+1>len(self.ref_seq[values[0]]):
                    raise ValueError(f'Amplicon {i+1} in existing amplicons file {existing_primers_bed} has coordinate {values[0]} {int(values[1])} greater than lenght of reference sequence: {self.config.reference_fasta}')
                primer_seq=self.ref_seq[values[0]][int(values[1]):int(values[2])+1] # +1 because python excludes the last item of slice
                direction_info=self._id_has_direction_info(values[0])
                if not direction_info[0] and len(values)>2:
                    direction_info=self._id_has_direction_info(values[3])
                    if not direction_info[0]: #this is worst case, add primer as both forward and reverse
                        primer=self._sequence_to_primer( primer_seq, True )
                        self.existing_primers.append(primer)
                        primer=self._sequence_to_primer( primer_seq, False )
                if direction_info[0]:
                    primer=self._sequence_to_primer( primer_seq, is_reverse= direction_info[1]=="reverse" )
                self.existing_primers.append(primer)

    # vcf_file_name=f'{root_dir}/HandyAmpliconTool/test_data/outputs/2_3_1/snps.vcf'
    # header_rows=0
    # for line in open(vcf_file_name):
    #     if line[0]=="#":
    #         header_rows+=1
    # vcf_data=pd.read_csv(vcf_file_name, sep="\t", header=0, skiprows=header_rows-1) # -1 because header row also starts with #CHROM


    def add_global_primer_args(self, for_seq: str, rev_seq:str, template: str):
        global_args={
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 55.0,
            'PRIMER_MAX_TM': 65.0,
            'PRIMER_PRODUCT_SIZE_RANGE': f'100-{len(template)}'}
        if for_seq!=None:
            global_args['SEQUENCE_PRIMER']=for_seq
        if rev_seq!=None:
            global_args['SEQUENCE_PRIMER_REVCOMP']=rev_seq
        return global_args
    
    def forms_homodimers(self, primer_sequence: str) -> bool:
        """Checks if sequence forms homodimer at above 10C
        :param primer_sequence: sequence of the primers
        :type primer_sequence: str
        :return: True if sequence forms homodimer
        :rtype: bool
        """
        thermoresult_forward=primer3.bindings.calc_homodimer(primer_sequence)
        return (thermoresult_forward.structure_found and thermoresult_forward.tm>10)
        
    def count_heterodimers(self, primer_sequence: str, orientation: str) -> int:
        """Checks if sequence forms heterodimers at above 10C
        with existing primers
        :param primer_sequence: sequence of the primers
        :type primer_sequence: str
        :param orientation: Forward, Reverse or Unknown
        :type orientation: str
        :return: number of heterodimers sequences forms with existing primers
        :rtype: int
        """
        total_structures=0
        for primer in self.existing_primers:
            if (orientation=="Forward" and not primer.is_reverse) or \
                (orientation=="Reverse" and primer.is_reverse) or \
                    (orientation=="Unknown"):
                temp_result=primer3.bindings.calc_heterodimer(primer_sequence, primer.seq)
                total_structures+=1 if temp_result.structure_found and temp_result.tm>10 else 0
        return total_structures
    
    def process_p3_output(self, output: Dict, id_prefix: str) -> List[PrimerPair]:
        """Converts the primer3 output format into list of PrimerPair objects
        with existing primers
        :param output: dictionary object generated by primer3
        :type output: Dict
        :param id_prefix: the name to assign to the each new PrimerPair
        :type id_prefix: str
        :return: number of heterodimers sequences forms with existing primers
        :rtype: int
        """
        primer_pairs=len(output["PRIMER_PAIR"])
        result: List[PrimerPair]=[]
        for i in range(0,primer_pairs):
            forward=Primer(seq=output[f'PRIMER_LEFT_{str(i)}_SEQUENCE'],
                        g_c=output[f'PRIMER_LEFT_{str(i)}_GC_PERCENT'],
                        t_m=output[f'PRIMER_LEFT_{str(i)}_TM'], is_reverse=False)
            reverse=Primer(seq=output[f'PRIMER_RIGHT_{str(i)}_SEQUENCE'],
                        g_c=output[f'PRIMER_RIGHT_{str(i)}_GC_PERCENT'],
                        t_m=output[f'PRIMER_RIGHT_{str(i)}_TM'], is_reverse=True)
            new_pair=PrimerPair(name=id_prefix, forward=forward, reverse=reverse)
            new_pair.penalty=output[f'PRIMER_PAIR_{str(i)}_PENALTY']
            if not self.forms_homodimers(forward.seq) and not self.forms_homodimers(reverse.seq) and self.count_heterodimers(forward.seq, reverse.seq)==0:
                result.append(new_pair)
                #output[f'PRIMER_PAIR_{str(i)}_PRODUCT_SIZE']
                # print(f'LEFT: {output["PRIMER_LEFT_EXPLAIN"]}')
                # print(f'RIGHT: {output["PRIMER_RIGHT_EXPLAIN"]}')
                # print(f'PAIR: {output["PRIMER_PAIR_EXPLAIN"]}')
                # print("\n")
        return result

    def _for_testing_load_gts(self, species_pkl, gts_pkl):
        with open(species_pkl, "rb") as pickled_file:
            species: Genotype = pickle.load(pickled_file)

        with open(gts_pkl, "rb") as pickled_file:
            self.genotypes: Genotypes = pickle.load(pickled_file)

        self.genotypes.genotypes.append(species)

    def _snps_within_interval(self, snps: List[SNP], ref_contig: int, ref_position: int,  interval_len: int) -> List[SNP]:
        result=[snp for snp in snps if snp.ref_contig_id==ref_contig and snp.position>ref_position-interval_len and snp.position<ref_position+interval_len]
        return result
    
    def _get_ref_sequence(self, seq_id, start, end) -> str:
        if seq_id not in self.ref_seq:
            raise ValueError(f'Sequence {seq_id} not found in reference FASTA {self.config.reference_fasta}')
        return str(self.ref_seq[seq_id][start:end])
    
    def _both_primers_given(self, left_snps: List[SNP], right_snps: List[SNP], target_snp:SNP) -> None:
        for left_snp, right_snp in product(left_snps, right_snps):
            #case 1
            distance_between_snps=right_snp.position-left_snp.position
            if distance_between_snps > self.config.flank_len_to_check or distance_between_snps < InputConfiguration.min_amplicon_length: #distance between SNPs is too long:
                #check the left and right separately
                self._left_primers_given(left_snps=[left_snp], target_snp=target_snp)
                self._right_primers_given(right_snps=[right_snp], target_snp=target_snp)
                continue
            forward=self._get_ref_sequence( left_snp.ref_contig_id, left_snp.position-1, left_snp.position+20  )
            reverse=self._get_ref_sequence( right_snp.ref_contig_id, right_snp.position-20,right_snp.position )
            reverse=str(Seq(reverse).reverse_complement())
            seq_args={"SEQUENCE_ID":f'{self.target_gt}_{str(target_snp.coordinate)}',
                    "SEQUENCE_TEMPLATE": self._get_ref_sequence( left_snp.ref_contig_id , left_snp.position-1, right_snp.position),
                    "SEQUENCE_INCLUDED_REGION": [0, right_snp.position-left_snp.position+1 ]}
            #add "fixed" to name to show that this pair has both sides fixed
            global_args=self.add_global_primer_args(for_seq=forward,rev_seq=reverse, template=seq_args["SEQUENCE_TEMPLATE"])
            primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
            self.new_primer_pairs += self.process_p3_output(primers, id_prefix=f'{self.target_gt}_{target_snp.coordinate}_both_fixed')

    def _right_primers_given(self, right_snps:List[SNP], target_snp: SNP) -> None:
        for right_snp in right_snps:
            reverse=self._get_ref_sequence( right_snp.ref_contig_id, right_snp.position-20,right_snp.position )
            reverse=str(Seq(reverse).reverse_complement())
            seq_args={"SEQUENCE_ID":f'{self.target_gt}_{str(target_snp.position)}',
                    "SEQUENCE_TEMPLATE": self._get_ref_sequence( right_snp.ref_contig_id, right_snp.position-self.config.flank_len_to_check ,right_snp.position ),
                    "SEQUENCE_INCLUDED_REGION": [0, self.config.flank_len_to_check ]}
            global_args=self.add_global_primer_args(for_seq=None,rev_seq=reverse,  template=seq_args["SEQUENCE_TEMPLATE"])
            primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
            self.new_primer_pairs += self.process_p3_output(primers, id_prefix=f'{self.target_gt}_{right_snp.position}_right_fixed')

    def _left_primers_given(self, left_snps: List[SNP], target_snp: SNP) -> None:
        for left_snp in left_snps:
            forward=self._get_ref_sequence( left_snp.ref_contig_id, left_snp.position-1, left_snp.position+20  )
            seq_args={"SEQUENCE_ID":f'{self.target_gt}_{str(target_snp.position)}',
                    "SEQUENCE_TEMPLATE":    self._get_ref_sequence(target_snp.ref_contig_id, left_snp.position-20, left_snp.position+self.config.flank_len_to_check),
                    "SEQUENCE_INCLUDED_REGION": [0, self.config.flank_len_to_check ]}
            global_args=self.add_global_primer_args(for_seq=forward, rev_seq=None, template=seq_args["SEQUENCE_TEMPLATE"])
            primers=primer3.bindings.design_primers(seq_args=seq_args, global_args=global_args)
            self.new_primer_pairs += self.process_p3_output(primers, id_prefix=f'{self.target_gt}_{target_snp.position}_both_fixed')


    def find_candidate_primers(self):
        self.target_gt="2.3.1"
        # for every SNP in target_lineage, identify the nearby SNPs 
        all_genotype_snps=[snp for genotype in self.genotypes.genotypes for snp in genotype.defining_snps if genotype.name!="species" and genotype.name!=self.target_gt]
        all_genotype_snps=sorted(all_genotype_snps, key=lambda x: (x.ref_contig_id, x.position) )
        all_species_snps=[snp for genotype in self.genotypes.genotypes for snp in genotype.defining_snps if genotype.name=="species"]
        temp=[genotype for genotype in self.genotypes.genotypes if genotype.name!="species"]
        all_species_snps=sorted(all_species_snps, key=lambda x: (x.ref_contig_id, x.position) )
        interval_len=self.config.flank_len_to_check

        target_genotype=self.genotypes.get_genotype(self.target_gt)
        
        self.new_primer_pairs.clear()
        for i, snp in enumerate(target_genotype.defining_snps):
            #print(i/len(target_genotype.defining_snps))
            other_gt_snps=self._snps_within_interval(all_genotype_snps, snp.ref_contig_id, snp.position, interval_len)
            if len(other_gt_snps)>0:
                #### !!!! Include this im primer selection process to prioritise primers which capture multiple GTs
                print([f.name for f in temp  if other_gt_snps[0] in f.defining_snps]) 
            species_gt_snps=self._snps_within_interval(all_species_snps, snp.ref_contig_id, snp.position, interval_len)
            if len(species_gt_snps)>0 and len(species_gt_snps)<30: # the target SNP has at least one flanking species SNPs, but too many is indicative of problematic region
                left_species_snps=[species_snp for species_snp in species_gt_snps if species_snp.position<snp.position]
                right_species_snps=[species_snp for species_snp in species_gt_snps if species_snp.position>snp.position]
            #four three cases: 
            # 1 - left and right serovar SNPs anchored primers
            # 2 - left serovar anchored primer
            # 3 - right serovar anchored primer
            # 4 - neither side has serovar primer, ignore this kind of SNP
                if len(left_species_snps)>0 and len(right_species_snps)>0:
                    self._both_primers_given(left_snps=left_species_snps, right_snps=right_species_snps, target_snp=snp)
                elif len(left_species_snps):
                    self._left_primers_given(left_snps=left_species_snps, target_snp=snp)
                elif len(right_species_snps):
                    self._right_primers_given(right_snps=right_species_snps, target_snp=snp)

    #         #sometimes the same pair of F/R sequences is generated multiple times because there are multiple species SNPs around same genotype SNP
    #         unique_new_primer_pairs: List[PrimerPair]=[]
    #         existing_sequences=set()
    #         for pair in new_primer_pairs:
    #             sequences=pair.forward.seq+"_"+pair.reverse.seq
    #             if sequences not in existing_sequences:
    #                 unique_new_primer_pairs.append(pair)
    #                 existing_sequences.add(sequences)
