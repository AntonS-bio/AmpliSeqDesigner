from collections import Counter
from typing import List, Dict
from generate_msa import MsaGenerator, MsaResult
from data_classes import Amplicon, SNP, FlankingAmplicon, Genotype, InputConfiguration
from tqdm import tqdm

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


    def msa_df_to_msa_file(self, msa: MsaResult, file_prefix: str) -> bool:
        """Converts pandas DataFrame with SNP data into MSA file
        This is optional and useful for later looking into various SNPs
        """
        with open(f'{self.msa_dir}/{file_prefix}.fasta', "w") as output_file:
            for i in range(0, msa.matrix.shape[0] ):
                output_file.write(f'>{msa.seq_ids[i]}'+"\n")
                output_file.write(msa.row_to_seq(i) + "\n")
        return True


    def get_bifurcating_snps(self, genotype: Genotype) -> Genotype:
        '''Identifies SNPs that separate target and non-target species around amplicon sequences'''
        msa_generator=MsaGenerator(temp_blast_db_dir=self.temp_blast_db_dir)
        
        msa_results: List[MsaResult]=msa_generator.generate_msa(genotype.amplicons, genomes_dir=self.negative_genomes_dir)

        print("Processing MSA data")
        with tqdm(total=len(msa_results)) as progress_meter:
            for msa in msa_results:
                amplicon_id=msa.amplicon_id
                progress_meter.update(1)
                self.msa_df_to_msa_file(msa, [f for f in genotype.amplicons if f.id==amplicon_id][0].name) ##this saves MSA files for fasta.
                current_amplicon=[f for f in genotype.amplicons if f.id==amplicon_id][0]
                ampicon_msa_seq: str= msa.row_to_seq(msa.seq_ids.index(amplicon_id))
                msa_to_amplicon_coord: Dict[int, int] =self._map_msa_to_ref_coordinates( msa_seq= list(ampicon_msa_seq) )

        
                if len(msa.seq_ids)==1:
                    current_amplicon.has_homologues=False
                else:
                    current_amplicon.has_homologues=True
                    for i in range(0, msa.matrix.shape[1]):
                        target_nucleotide=ampicon_msa_seq[i]
                        if target_nucleotide=="-":
                            continue #can't target primer to non-existent nucleotide
                        bases_at_position=Counter( msa.nucleotides_at_col(i) )
                        if bases_at_position[target_nucleotide]<=InputConfiguration.max_matching_negative_genomes: #i.e. the target strain nucleotide is unique among all other strains
                            #four cases exist (here, T, C and A can be any nucleotide):
                            # 1) A vs TTTT : output T
                            # 2) A vs TTCC : output T and C
                            # 3) A vs ---- : output -
                            # 4) A vs TT-- : output T only, not - 
                            for alt_base, count in bases_at_position.items():
                                if alt_base!=target_nucleotide:
                                    if alt_base=="-" and len(bases_at_position)!=2: #case 4
                                        continue #do not output missing base if there are other nucleotides in that position
                                    snp=SNP(ref_contig_id=current_amplicon.ref_seq.refseq_id, ref_base=target_nucleotide, alt_base=alt_base,  position=current_amplicon.ref_seq.ref_start+msa_to_amplicon_coord[i])
                                    if snp.alt_base=="-":
                                        if i==0:
                                            continue #exceptional case where deletion is the first base on amplicon and correct VCF cannot be created
                                        snp.ref_base="".join( [ ampicon_msa_seq[f] for f in range(i-1,i+1) ] ) #For deletions, reference is sequence from previous to deleted base
                                        snp.alt_base=ampicon_msa_seq[i-1]
                                    snp.passes_filters=True
                                    snp.specificity=1
                                    snp.sensitivity=1
                                    snp.is_species_snp=True
                                    genotype.add_genotype_allele(snp, snp.alt_base, bases_at_position[alt_base])
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

    def generate_flanking_amplicons(self) -> Genotype:
        """Identifies SNPs within left and right amplicon flanking sequences
        The flanking sequences are max_len-amplicon_len which means that total length
        of flanking plus amplicon sequences may be longer than max_seq_len
        """        
        genotype: Genotype=Genotype(InputConfiguration.SPECIES_NAME)
        i=0
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
    

    
