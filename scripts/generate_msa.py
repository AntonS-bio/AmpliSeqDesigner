#take all fastas in specified directory and check all for specific gene
from os import listdir, walk, mkdir, remove
from os.path import isfile, join, splitext, exists
import subprocess
from typing import List, Dict, Tuple
from run_blast import BlastRunner
from multiprocessing import Pool
import numpy as np
import numpy.typing as npt
from Bio.Seq import Seq
from data_classes import Amplicon, BlastResult, InputConfiguration
from tqdm import tqdm

class MergedAmplicons:
    """
    Class designed to merge and demerge multiple overlapping amplicons into fewer
    e.g. ____________
         _______
         ____
    Only amplicons which are complete subset of another are merged, because 
    otherwise the amplicon length may become very long
    The merge reduces the number of redundant BLAST searches and 
    MSA alignments.
    """

    def __init__(self) -> None:
        self._source_amplicons: List[Amplicon] = []
        self._source_to_destination: Dict[str, Amplicon]
        self._destination_amplicons: List[Amplicon] = []

    @property
    def source_amplicons(self) -> List[Amplicon]:
        return self._source_amplicons
    
    @property
    def destination_amplicons(self) -> List[Amplicon]:
        return self._destination_amplicons
   
    def get_destination_amplicon(self, amplicon: Amplicon) -> Amplicon:
        if amplicon.id not in self._source_to_destination:
            raise ValueError(f'Amplicon {amplicon.name} not present in merged amplicons')
        return self._source_to_destination[amplicon.id]

    def merge_amplicons(self, amplicons: List[Amplicon]) -> List[Amplicon]:
        """Identifies ovelapping amplicons, checks if one is a subset of another
        and keeps the longer amplicon

        :param amplicon: list of amplicons to merge
        :type amplicon: List[Amplicon]
        """
        #sorting guarantees
        amplicons=sorted(amplicons, key=lambda x: x.len, reverse=True)
        self._source_to_destination={}

        for i, first_amplicon in enumerate(amplicons): #Rewrite later 
            for second_amplicon in amplicons[i+1:]:
                if  second_amplicon.id not in self._source_to_destination.keys() and \
                    first_amplicon.coord_in_amplicon( (second_amplicon.ref_contig, second_amplicon.ref_seq.ref_start) ) and \
                    first_amplicon.coord_in_amplicon( (second_amplicon.ref_contig, second_amplicon.ref_seq.ref_end) ):
                        self._source_to_destination[second_amplicon.id] = first_amplicon.id
        #return those amplicons that are not a subset of another amplicon
        self._destination_amplicons = [f for f in amplicons if f.id not in self._source_to_destination.keys()]
        return self._destination_amplicons

class MsaResult:
    """Result of MSA alignment consiting of two parts:
    Index of sequences IDs and MSA sequences
    """
    def __init__(self, amplicon_id: str, ids: List[str], sequences:List[str]) -> None:
        self._amplicon_id=amplicon_id
        self._ids: List[str]=ids
        self._sequences: npt.NDArray= np.asarray([ self._to_numeric( list(f.upper()) ) for f in sequences ], dtype=int)

    def _to_numeric(self, sequence:List[str]) -> List[int]:
        return [ InputConfiguration.BASE_DIC[f] if f in InputConfiguration.BASE_DIC else InputConfiguration.BASE_DIC["N"] for f in sequence ]
    
    def _to_char(self, sequence: List[int]) -> str:
        return "".join( [ InputConfiguration.NUMBER_DIC[f] for f in sequence ] )

    @property
    def amplicon_id(self) -> str:
        return self._amplicon_id

    @property
    def seq_ids(self) -> List[str]:
        return self._ids

    @property
    def matrix(self) -> npt.NDArray:
        return self._sequences

    def _values_at_col(self, index:int) -> npt.NDArray:
        if index < self._sequences.shape[1]:
            return self._sequences[:,index]
        else:
            raise ValueError(f'Index value {index} exceeds the number of MSA columns')

    def nucleotides_at_col(self, index:int) -> List[str]:
        numeric_values=self._values_at_col(index)
        return  [ InputConfiguration.NUMBER_DIC[f] for f in numeric_values ]

    def row_to_seq(self, index: int) -> str:
        if index < self._sequences.shape[0]:
            return self._to_char(self._sequences[index,:])
        else:
            raise ValueError(f'Index value {index} exceeds the number of MSA rows')


class MsaGenerator:
    #use_sub_dirs=False
    #cpu_threads=8

    def __init__(self, temp_blast_db_dir: str) -> None:
        self.temp_blast_db_dir=temp_blast_db_dir
        self.file_to_search=[]

    def _get_fasta_files(self, dir_to_search: str):
        self.file_to_search=[]
        if InputConfiguration.use_negative_genomes_subdir:
            for path, subdirs, dir_files in walk(dir_to_search):
                for name in dir_files:
                    if isfile(join(path, name)) and (splitext(name)[-1]==".fasta" or splitext(name)[-1]==".fna"):
                        self.file_to_search.append(join(path, name))
        else:
            self.file_to_search = [dir_to_search+"/"+f for f in listdir(dir_to_search) if isfile(join(dir_to_search, f)) and (splitext(f)[-1]==".fasta" or splitext(f)[-1]==".fna")]

    def generate_msa(self, amplicons:List[Amplicon], genomes_dir:str) -> Dict[str, MsaResult]:
        """Takes list of Amplicons and directory of genomes
        blasts amplicon sequences against all genomes in directory 
        and creates multiple sequences alignment file of blast results
        one MSA per amplicon supplied
        """
        if not exists(InputConfiguration.output_dir):
            mkdir(InputConfiguration.output_dir)
        #Collect fasta files against which to run blast
        self._get_fasta_files(genomes_dir)
        if len(self.file_to_search)==0:
            if InputConfiguration.use_negative_genomes_subdir:
                raise ValueError(f'When looking for genomes to BLAST against, no .fna or .fasta files found in {genomes_dir} or sub-directories.')
            else:
                raise ValueError(f'When looking for genomes to BLAST against, no .fna or .fasta files found in {genomes_dir}. Did you mean to include sub-directories?')

        print("Merging amplicons")
        merged_amplicons=MergedAmplicons()
        merged_amplicons.merge_amplicons(amplicons)
        blast_results_raw=self._run_blast( merged_amplicons.destination_amplicons, self.file_to_search )
        blast_results=self._process_blast_results(blast_results_raw, merged_amplicons.destination_amplicons)
        msa_dfs: List[MsaResult] = self._align_blast_results(blast_results, merged_amplicons.destination_amplicons)
        return msa_dfs

    def _align_blast_results(self, blast_results: Dict[str, List[BlastResult]], amplicons: List[Amplicon]) -> List[MsaResult]:
        if __name__ == 'generate_msa':
            print("Generating MSAs")
            pool = Pool(processes= InputConfiguration.cpu_threads )
            aligner_inputs:List=[]
            for amplicon_id, amplicon_blast_results in blast_results.items():
                if len(amplicon_blast_results)==0:
                    continue
                amplicon_seq=[f.seq for f in amplicons if f.id==amplicon_id][0]
                aligner_inputs.append([amplicon_blast_results,amplicon_id,amplicon_seq])

            msa_results = list(tqdm( pool.imap(func=self._align_results_helper, iterable=aligner_inputs), total=len(aligner_inputs) ))
            pool.close()
            return msa_results
            # msa_dfs: Dict[str, pd.DataFrame]={}
            # for amplicon_id, amplicon_blast_results in blast_results.items(): #this is a shortcut and a more robust solution is required in longer-term
            #     for result in msa_results:
            #         if amplicon_id in result:
            #             msa_dfs[amplicon_id]=self._msa_to_dataframe(result)
            #             break

            # return msa_dfs

    def _run_blast(self, subject_sequences: List[Amplicon], query_files: List[str]) -> List[BlastResult] :
        """Runs blast against a single file at a time using Pool
        this is the fastest implemintaiton of all tried
        """
        if not exists(self.temp_blast_db_dir):
            mkdir(self.temp_blast_db_dir)
        with open(self.temp_blast_db_dir+"/temp.fasta", "w") as temp_file:
            for amplicon in subject_sequences:
                temp_file.write(">"+amplicon.id+"\n"+amplicon.seq+"\n")
        subject_seq_file=self.temp_blast_db_dir+"/temp.fasta"

        blast_runner=BlastRunner()           
        blast_runner.db_from_file(subject_seq_file, self.temp_blast_db_dir)
        blast_results: List[ Tuple[str, List[BlastResult]] ]=[]

        if __name__ == 'generate_msa':
            print("Running BLAST against genomes")
            with Pool(processes= InputConfiguration.cpu_threads) as pool: #min(  max(cpu_count()-1,1) , self.cpu_threads ) )
                blast_results = list(tqdm( pool.imap(func=blast_runner.run_from_file, iterable=query_files), total=len(query_files) ))
            return [item for sublist in blast_results for item in sublist]

    def _process_blast_results(self, blast_resuls: List[BlastResult], target_amplicons: List[Amplicon]) -> Dict[str, List[BlastResult] ]:  #str is the name of the amplicon
        """Create a dictionarty which for each amplicon ID lists
        valid blast hits with correctly oriented sequence
        """
        valid_amplicon_hits: Dict[str, List[BlastResult] ]={}
        amplicon_len_delta={}
        amplicon_len={}
        for amplicon in target_amplicons:
            amplicon_len[amplicon.id]=len(amplicon.seq)
            amplicon_len_delta[amplicon.id]=int(len(amplicon.seq) * (InputConfiguration.max_blast_length_diff/100))
            valid_amplicon_hits[amplicon.id]=[]

        for result in blast_resuls:
            #check lenght and identity
            if result.q_hit_len >= (amplicon_len[result.sseqid] - amplicon_len_delta[result.sseqid]) and \
                result.q_hit_len <= (amplicon_len[result.sseqid] + amplicon_len_delta[result.sseqid]) and \
                result.pident >= InputConfiguration.min_blast_identity:
                valid_amplicon_hits[result.sseqid].append(result)

        del blast_resuls
        return valid_amplicon_hits

    def _align_results_helper(self, values:List) -> MsaResult:
        """Take valid blast results and create a file of unalligned hits
        use MSA tool (here Mafft) to align them
        this might be replaced later, but at the moment this is simpler approach
        """
        blast_results: List[BlastResult]; amplicon_id: str; amplicon_seq: str
        blast_results, amplicon_id, amplicon_seq=values
        fasta_file=f'{self.temp_blast_db_dir}/{amplicon_id}.fasta'
        with open(fasta_file, "w") as output:
            output.write(">"+amplicon_id+"\n"+amplicon_seq+"\n")
            for result in blast_results:
                if result.sstart>result.send:
                    seq_to_allign=Seq(result.qseq.replace("-","")).reverse_complement()
                else:
                    seq_to_allign=result.qseq.replace("-","")
                output.write(f'>{result.qseqid}'+"\n")
                output.write(str(seq_to_allign)+"\n")

        outcome=subprocess.run(f'mafft --retree 1 {fasta_file}', \
                               shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        remove(fasta_file) 
        if outcome.returncode==2 or outcome.returncode==1:
            raise OSError(f'Error generating BLAST database: {outcome.stderr}')

        # msa_results: Dict[str,str]={}
        ids=[]
        sequences=[]
        for line in outcome.stdout.decode().strip().split("\n"):
            if line[0]==">":
                ids.append(line[1:])
                sequences.append("")
            else:
                sequences[-1]=sequences[-1]+line

        return MsaResult( amplicon_id, ids, sequences)

    # def _msa_to_dataframe(self, msa_result: Dict[str, str]) -> np.ndarray:
    #     indices=list(msa_result.keys())
    #     columns=len(msa_result[indices[0]])
    #     result=np.zeros( (len(indices), len(columns))  )
    #     #msa_df=pd.DataFrame(index=indices, columns=['pos_'+str(f) for f in range(0,columns)], dtype=str)
    #     for seq_id, sequence in msa_result.items():
    #         result[]
    #         msa_df.loc[seq_id]=list(str(sequence).upper())
    #     return msa_df
