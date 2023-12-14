#take all fastas in specified directory and check all for specific gene
from os import listdir, walk, mkdir, remove
from os.path import isfile, join, splitext, exists
import subprocess
from typing import List, Dict, Tuple
from run_blast import BlastRunner
from multiprocessing import cpu_count, Pool
import pandas as pd
#from Bio import SeqIO
from Bio.Seq import Seq
from data_classes import Amplicon, BlastResult
from tqdm import tqdm

class MsaGenerator:
    use_sub_dirs=False
    cpu_threads=8

    def __init__(self, temp_blast_db_dir: str) -> None:
        self.temp_blast_db_dir=temp_blast_db_dir
        self.file_to_search=[]

    def _get_fasta_files(self, dir_to_search: str):
        self.file_to_search=[]
        if self.use_sub_dirs:
            for path, subdirs, dir_files in walk(dir_to_search):
                for name in dir_files:
                    if isfile(join(path, name)) and (splitext(name)[-1]==".fasta" or splitext(name)[-1]==".fna"):
                        self.file_to_search.append(join(path, name))
        else:
            self.file_to_search = [dir_to_search+"/"+f for f in listdir(dir_to_search) if isfile(join(dir_to_search, f)) and (splitext(f)[-1]==".fasta" or splitext(f)[-1]==".fna")]

    def generate_msa(self, amplicons:List[Amplicon], output_dir: str, genomes_dir:str, max_blast_length_diff: int,  min_blast_identity:int) -> Dict[str, pd.DataFrame]:
        self.max_blast_length_diff=max_blast_length_diff
        self.min_blast_identity=min_blast_identity
        self.output_dir=output_dir
        if not exists(self.output_dir):
            mkdir(self.output_dir)
        #Collect fasta files against which to run blast
        self._get_fasta_files(genomes_dir)
        if len(self.file_to_search)==0:
            if self.use_sub_dirs:
                raise ValueError(f'No .fna or .fasta files found in {genomes_dir} or sub-directories.')
            else:
                raise ValueError(f'No .fna or .fasta files found in {genomes_dir}. Did you mean to include sub-directories?')
        #self.file_to_search=self.file_to_search
        blast_results_raw=self._run_blast(amplicons, self.file_to_search)
        blast_results=self._process_blast_results(blast_results_raw, amplicons)




        if __name__ == 'generate_msa':
            print("Generating MSAs")
            pool = Pool(processes= min(  max(cpu_count()-1,1) , self.cpu_threads ) )
            aligner_inputs:List=[]
            for amplicon_id, amplicon_results in blast_results.items():
                if len(amplicon_results)==0:
                    continue
                amplicon_seq=[f.seq for f in amplicons if f.id==amplicon_id][0]
                aligner_inputs.append([amplicon_results,amplicon_id,amplicon_seq])

            msa_results = list(tqdm( pool.imap(func=self._align_results, iterable=aligner_inputs), total=len(aligner_inputs) ))
            msa_dfs: Dict[str, pd.DataFrame]={}
            for amplicon_id, amplicon_results in blast_results.items(): #this is a shortcut and a more robust solution is required in longer-term
                for result in msa_results:
                    if amplicon_id in result:
                        msa_dfs[amplicon_id]=self._msa_to_dataframe(result)
                        break

            return msa_dfs

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
            pool = Pool(processes= min(  max(cpu_count()-1,1) , self.cpu_threads ) )
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
            amplicon_len_delta[amplicon.id]=int(len(amplicon.seq) * (self.max_blast_length_diff/100))
            valid_amplicon_hits[amplicon.id]=[]

        for result in blast_resuls:
            #check lenght and identity
            if result.q_hit_len >= (amplicon_len[result.sseqid] - amplicon_len_delta[result.sseqid]) and \
                result.q_hit_len <= (amplicon_len[result.sseqid] + amplicon_len_delta[result.sseqid]) and \
                result.pident >= self.min_blast_identity:
                valid_amplicon_hits[result.sseqid].append(result)

        del blast_resuls
        return valid_amplicon_hits

    def _align_results(self, values:List):
        """Take valid blast results and create a file of unalligned hits
        use MSA tool (here Mafft) to align them
        this might be replaced later, but at the moment this is simpler approach
        """
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

        msa_results: Dict[str,str]={}
        current_id=""
        current_seq=""
        for line in outcome.stdout.decode().strip().split("\n"):
            if line[0]==">":
                if current_id!="":
                    msa_results[current_id]=current_seq
                current_id=line[1:]
                current_seq=""
            else:
                current_seq=current_seq+line

        return msa_results

    def _msa_to_dataframe(self, msa_result: Dict[str, str]) -> pd.DataFrame:
        indices=list(msa_result.keys())
        columns=len(msa_result[indices[0]])
        msa_df=pd.DataFrame(index=indices, columns=['pos_'+str(f) for f in range(0,columns)], dtype=str)
        for seq_id, sequence in msa_result.items():
            msa_df.loc[seq_id]=list(str(sequence).upper())
        return msa_df
