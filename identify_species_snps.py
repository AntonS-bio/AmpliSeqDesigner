from generate_msa import MsaGenerator
import pandas as pd
from collections import Counter
from Bio import SeqIO
from typing import List, Dict
from data_classes import Amplicon, SNP


class IdentifySpeciesSnps:

    def __init__(self, ref_fasta: str, negative_genomes_dir: str, msa_dir: str, amplicons_bed:str, temp_blast_db_dir: str) -> None:
        self.ref_fasta=ref_fasta
        self.amplicons_bed=amplicons_bed
        self.negative_genomes_dir=negative_genomes_dir
        self.msa_dir=msa_dir
        self.temp_blast_db_dir=temp_blast_db_dir

    def _get_bifurcating_snps(self, amplicons: List[Amplicon]) -> List[SNP]:
        msa_generator=MsaGenerator(temp_blast_db_dir=self.temp_blast_db_dir)

        msa_dfs: Dict[str, pd.DataFrame]=msa_generator.generate_msa(amplicons, output_dir=self.msa_dir, genomes_dir=self.negative_genomes_dir)


        #### To Continue from there -> This needs to produces list of SNPs
        for amplicon_id, msa_df in msa_dfs.items():
            segregating_columns=[]
            print(msa_df.shape)
            if msa_df.shape[0]==1:
                print(f'{amplicon_id} has no homologues')        
            else:
                for col in msa_df.columns:
                    target_nucleotide=msa_df.loc[amplicon_id,col]
                    if Counter(msa_df[col])[target_nucleotide]==1:
                        segregating_columns.append(col)
                print(f'{amplicon_id} has {str(len(segregating_columns))}')

    # def identify_insequence_snps(self):
    #     """Identifies SNPs within amplicons that distinguish the amplicon from other 
    #     genotypes. 
    #     """        
    #     amplicons: List[Amplicon]=[] 
    #     with open(f'{self.amplicons_bed}') as input_bed_file: #don't keep the file open, hence why load it to memory
    #         for line in input_bed_file:
    #             bed_line_values = line.strip().split("\t")
    #             ampl_chr, ampl_start, ampl_end=bed_line_values[0:3]
    #             ampl_start=int(ampl_start)
    #             ampl_end=int(ampl_end)
    #             amplicon_id='_'.join( [str(f) for f in bed_line_values[0:4] ] )
    #             for record in SeqIO.parse(self.ref_fasta,"fasta"):
    #                 if record.id==ampl_chr:
    #                     new_amplicon=Amplicon(amplicon_id, str(record.seq[ampl_start:ampl_end]))
    #                     new_amplicon.middle=True
    #                     amplicons.append( new_amplicon )
    #                     continue
        
    #     self._get_bifurcating_snps(amplicons)

    def identify_flanking_snps(self, max_seq_len: int) -> List[SNP]:
        """Identifies SNPs within left and right amplicon flanking sequences
        The flanking sequences are max_len-amplicon_len which means that total length
        of flanking plus amplicon sequences may be longer than max_seq_len
        """        
        amplicons: List[Amplicon]=[]
        with open(f'{self.amplicons_bed}') as input_bed_file: #don't keep the file open, hence why load it to memory
            for line in input_bed_file:
                bed_line_values = line.strip().split("\t")
                ampl_chr, ampl_start, ampl_end=bed_line_values[0:3]
                ampl_start=int(ampl_start)
                ampl_end=int(ampl_end)
                flank_len=max(max_seq_len-(ampl_end-ampl_start),0)
                if flank_len==0:
                    continue
                amplicon_id='_'.join( [str(f) for f in bed_line_values[0:4] ] )
                for record in SeqIO.parse(self.ref_fasta,"fasta"):
                    if record.id==ampl_chr:
                        #left flanking
                        new_amplicon=Amplicon(amplicon_id, str(record.seq[max(ampl_start-flank_len,0):ampl_start]))
                        new_amplicon.left_flank=True
                        amplicons.append( new_amplicon)
                        new_amplicon=Amplicon(amplicon_id, str(record.seq[ampl_end:min(ampl_end+flank_len,len(record.seq)) ]))
                        new_amplicon.right_flank=True
                        amplicons.append( new_amplicon )
                        continue
        
        result:List[SNP] = self._get_bifurcating_snps(amplicons[0:20]) ####! DEBUG ONLY
        return result