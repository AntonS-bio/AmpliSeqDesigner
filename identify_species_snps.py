from generate_msa import MsaGenerator
import pandas as pd
from collections import Counter
from Bio import SeqIO
from typing import List, Dict, Tuple
from data_classes import Amplicon, SNP, FlankingAmplicon


class IdentifySpeciesSnps:

    def __init__(self, ref_fasta: str, negative_genomes_dir: str, msa_dir: str, amplicons_bed:str, temp_blast_db_dir: str) -> None:
        self.ref_fasta=ref_fasta
        self.amplicons_bed=amplicons_bed
        self.negative_genomes_dir=negative_genomes_dir
        self.msa_dir=msa_dir
        self.temp_blast_db_dir=temp_blast_db_dir

    def _get_bifurcating_snps(self, amplicons: List[Amplicon]) -> List[SNP]:
        msa_generator=MsaGenerator(temp_blast_db_dir=self.temp_blast_db_dir)

        msa_dfs: Dict[str, pd.DataFrame]=msa_generator.generate_msa(amplicons, output_dir=self.msa_dir, genomes_dir=self.negative_genomes_dir, 
                                                                        max_blast_length_diff=self.max_blast_length_diff,
                                                                        min_blast_identity=self.min_blast_identity)

        segregating_snps:List[SNP]=[]

        for amplicon_id, msa_df in msa_dfs.items():
            current_amplicon=[f for f in amplicons if f.id==amplicon_id][0]
            ampicon_msa_seq: List[str]=[f for f in msa_df.loc[amplicon_id]]
            msa_to_amplicon_coord: Dict[int, int] =self._map_msa_to_ref_coordinates(amplicon=current_amplicon, msa_seq=ampicon_msa_seq )

      
            if msa_df.shape[0]==1:
                current_amplicon.has_homologues=False
            else:
                current_amplicon.has_homologues=True
                for i, col in enumerate(msa_df.columns):
                    target_nucleotide=msa_df.loc[amplicon_id,col]
                    bases_at_position=Counter(msa_df[col])
                    if bases_at_position[target_nucleotide]==1:
                        snp=SNP(ref_contig_id=current_amplicon.ref_contig, ref_base=target_nucleotide,  position=current_amplicon.ref_start+msa_to_amplicon_coord[i])
                        snp.alt_base=bases_at_position.most_common(1)[0][0] #.most_common() returns list of tuples hence [0][0]
                        snp.passes_filters=True
                        snp.specificity=1
                        snp.sensitivity=1
                        snp.is_species_snp=True
                        segregating_snps.append(snp)
        return segregating_snps

    def _map_msa_to_ref_coordinates(self, amplicon:Amplicon, msa_seq:List[str])-> Dict[int, int]:
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


    def identify_insequence_snps(self,  max_blast_length_diff: int,  min_blast_identity:int) -> List[Amplicon]:
        """Identifies SNPs within amplicons that distinguish the amplicon from other 
        genotypes. 
        """        
        self.max_blast_length_diff=max_blast_length_diff
        self.min_blast_identity=min_blast_identity
        amplicons: List[Amplicon]=[]
        with open(f'{self.amplicons_bed}') as input_bed_file: #don't keep the file open, hence why load it to memory
            for line in input_bed_file:
                new_amplicon=Amplicon.from_bed_line(line,self.ref_fasta)
                amplicons.append(new_amplicon)
        
        result:List[SNP] = self._get_bifurcating_snps(amplicons)

        for snp in result:
            for amplicon in amplicons:
                if amplicon.snp_in_amplicon(snp):
                    amplicon.snps.append(snp)

        return amplicons

    def identify_flanking_snps(self, max_seq_len: int, max_blast_length_diff: int,  min_blast_identity:int) -> List[FlankingAmplicon]:
        """Identifies SNPs within left and right amplicon flanking sequences
        The flanking sequences are max_len-amplicon_len which means that total length
        of flanking plus amplicon sequences may be longer than max_seq_len
        """        
        self.max_blast_length_diff=max_blast_length_diff
        self.min_blast_identity=min_blast_identity
        amplicons: List[FlankingAmplicon]=[]
        with open(f'{self.amplicons_bed}') as input_bed_file: #don't keep the file open, hence why load it to memory
            for line in input_bed_file:
                new_amplicon=Amplicon.from_bed_line(line,self.ref_fasta)
                left_flanking=FlankingAmplicon.from_parent_bed_line(self.ref_fasta, True, max_seq_len, new_amplicon)
                right_flanking=FlankingAmplicon.from_parent_bed_line(self.ref_fasta, False, max_seq_len, new_amplicon)
                amplicons.append( left_flanking)
                amplicons.append( right_flanking)
        
        result:List[SNP] = self._get_bifurcating_snps(amplicons)

        for snp in result:
            for amplicon in amplicons:
                if amplicon.snp_in_amplicon(snp):
                    amplicon.snps.append(snp)

        return amplicons
    
