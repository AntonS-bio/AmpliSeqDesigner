from collections import Counter
from typing import List, Dict
import pandas as pd
from generate_msa import MsaGenerator
from data_classes import Amplicon, SNP, FlankingAmplicon, Genotype, InputConfiguration


class IdentifySpeciesSnps:
    """Set of functions to identify SNPs that separate target organism
    from non-target organisms. Works from Multiple Sequence Alignment file
    """
    def __init__(self, ref_fasta: str, negative_genomes_dir: str, msa_dir: str, amplicons_bed:str, temp_blast_db_dir: str) -> None:
        self.ref_fasta=ref_fasta
        self.amplicons_bed=amplicons_bed
        self.negative_genomes_dir=negative_genomes_dir
        self.msa_dir=msa_dir
        self.temp_blast_db_dir=temp_blast_db_dir

    @classmethod
    def from_config(cls, config: InputConfiguration) -> InputConfiguration:
        """Constructor using bedfile lines. 
        :param config: the config input
        :type config: InputConfiguration
        """
        species_snp_finder=cls(config.reference_fasta, config.negative_genomes, config.msa_dir, config.multi_gt_intervals, config.temp_blast_db)
        return species_snp_finder


    def msa_df_to_msa_file(self, msa_df: pd.DataFrame, file_prefix: str) -> bool:
        """Converts pandas DataFrame with SNP data into MSA file
        This is optional and useful for later looking into various SNPs
        """
        with open(f'{self.msa_dir}/{file_prefix}.fasta', "w") as output_file:
            for index in msa_df.index:
                output_file.write(f'>{index}'+"\n")
                output_file.write(f'{"".join( [str(f) for f in msa_df.loc[index]] ) }'+"\n")
        return True


    def get_bifurcating_snps(self, genotype: Genotype) -> Genotype:
        '''Identifies SNPs that separate target and non-target species around amplicon sequences'''
        msa_generator=MsaGenerator(temp_blast_db_dir=self.temp_blast_db_dir)
        
        msa_dfs: Dict[str, pd.DataFrame]=msa_generator.generate_msa(genotype.amplicons, genomes_dir=self.negative_genomes_dir)


        for amplicon_id, msa_df in msa_dfs.items():
            self.msa_df_to_msa_file(msa_df, [f for f in genotype.amplicons if f.id==amplicon_id][0].name) ##this saves MSA files for fasta.
            current_amplicon=[f for f in genotype.amplicons if f.id==amplicon_id][0]
            ampicon_msa_seq: List[str]=[f for f in msa_df.loc[amplicon_id]]
            msa_to_amplicon_coord: Dict[int, int] =self._map_msa_to_ref_coordinates( msa_seq=ampicon_msa_seq )

    
            if msa_df.shape[0]==1:
                current_amplicon.has_homologues=False
            else:
                current_amplicon.has_homologues=True
                for i, col in enumerate(msa_df.columns):
                    target_nucleotide=msa_df.loc[amplicon_id,col]
                    if target_nucleotide=="-":
                        continue #can't target primer to non-existent nucleotide
                    bases_at_position=Counter(msa_df[col])
                    if bases_at_position[target_nucleotide]==1: #i.e. the target strain nucleotide is unique among all other strains
                        #four cases exist (here, T, C and A can be any nucleotide):
                        # 1) A vs TTTT : output T
                        # 2) A vs TTCC : output T and C
                        # 3) A vs ---- : output -
                        # 4) A vs TT-- : output T only, not - 
                        for alt_base, count in bases_at_position.items():
                            if alt_base!=target_nucleotide:
                                if alt_base=="-" and len(bases_at_position)!=2: #case 4
                                    continue #do not output missing base if there are other nucleotides in that position
                                snp=SNP(ref_contig_id=current_amplicon.ref_contig, ref_base=target_nucleotide, alt_base=alt_base,  position=current_amplicon.ref_start+msa_to_amplicon_coord[i])
                                if snp.alt_base=="-":
                                    snp.alt_base="."
                                snp.passes_filters=True
                                snp.specificity=1
                                snp.sensitivity=1
                                snp.is_species_snp=True
                                genotype.add_genotype_allele(snp, snp.ref_base, bases_at_position[alt_base])
        return genotype

    def _map_msa_to_ref_coordinates(self, msa_seq:List[str])-> Dict[int, int]:
        """MSA produces a sequence string with gaps so position of nucleotide in the string
        does not correspond to the position of nucleotide in reference sequence. 
        This function creates a mapping that allows SNP object to have coordinates of the reference
        sequence.
        Where msa position has missing value, the dictionary value is -1
        The full amplicon sequence is included in MSA, see generate_msa.generate_msa and  generate_msa._align_results
        """
        amplicon_position=0
        msa_to_amplicon_position: Dict[int,int]={}
        for msa_position, position_value in enumerate(msa_seq):
            if position_value!="-":
                msa_to_amplicon_position[msa_position]=amplicon_position
                amplicon_position+=1
            else:
                msa_to_amplicon_position[msa_position]=-1
        return msa_to_amplicon_position

    # def identify_insequence_snps(self,  max_blast_length_diff: int,  min_blast_identity:int) -> List[Amplicon]:
    #     """Identifies SNPs within amplicons that distinguish the amplicon from other 
    #     genotypes. 
    #     """        
    #     self.max_blast_length_diff=max_blast_length_diff
    #     self.min_blast_identity=min_blast_identity
    #     genotype: Genotype=Genotype("species")

    #     with open(f'{self.amplicons_bed}') as input_bed_file: #don't keep the file open, hence why load it to memory
    #         for line in input_bed_file:
    #             new_amplicon=Amplicon.from_bed_line(line,self.ref_fasta)
    #             genotype.amplicons.append(new_amplicon)
        
    #     result:Genotype = self.get_bifurcating_snps(genotype)

    #     return result

    def generate_flanking_amplicons(self) -> Genotype:
        """Identifies SNPs within left and right amplicon flanking sequences
        The flanking sequences are max_len-amplicon_len which means that total length
        of flanking plus amplicon sequences may be longer than max_seq_len
        """        
        genotype: Genotype=Genotype("species")

        with open(f'{self.amplicons_bed}') as input_bed_file: #don't keep the file open, hence why load it to memory
            for line in input_bed_file:
                if line.strip()=="":#catches the last line of bed file that can be just \n which user might not realise
                    continue
                new_amplicon=Amplicon.from_bed_line(line,self.ref_fasta)
                genotype.amplicons.append(new_amplicon)
                if InputConfiguration.flank_len_to_check>0:
                    left_flanking=FlankingAmplicon.from_parent_bed_line(self.ref_fasta, True, InputConfiguration.flank_len_to_check, new_amplicon)
                    right_flanking=FlankingAmplicon.from_parent_bed_line(self.ref_fasta, False, InputConfiguration.flank_len_to_check, new_amplicon)
                    genotype.amplicons.append( left_flanking)
                    genotype.amplicons.append( right_flanking)
        
        return genotype
    

    
