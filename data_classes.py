from typing import List, Tuple, Dict
import pandas as pd
import uuid
from io import TextIOWrapper
from Bio import SeqIO
import copy

class SNP:
    def __init__(self, **kwargs) -> None:
        """Constructor

        :param ref_contig_id: The id of the contig which is the reference for this SNP, defaults to empty string
        :type ref_contig_id: str, optional

        :param ref_base: The nucleotide at this position on the reference contig, defaults to empty string
        :type ref_base: str, optional

        :param alt_base: The alternative nucleotide at this position on the reference contig, defaults to empty string
        :type alt_base: str, optional

        :param position: The nucleotide at this position on the reference contig, defaults to empty string
        :type position: int, optional
        
        :param passes_filters: Indicates if SNP passes some filters, defaults to False
        :type passes_filters: bool, optional
        """
        if kwargs.get("ref_contig_id","")!="":
            self._ref_contig_id=kwargs.get("ref_contig_id","")
        if kwargs.get("ref_base","")!="":
            self._ref_base=kwargs.get("ref_base","")
        if kwargs.get("alt_base","")!="":
            self._alt_base=kwargs.get("alt_base","")
        if kwargs.get("position","")!="":
            self._position=kwargs.get("position","")
        self._passes_filters=kwargs.get("passes_filters",False)

    @property
    def coordinate(self) -> Tuple[str, int]:
        return (self._ref_contig_id, self._position)

    @property
    def ref_contig_id(self) -> str:
        return self._ref_contig_id

    @ref_contig_id.setter
    def ref_contig_id(self, value: str):
        self._ref_contig_id = value

    @property
    def position(self) -> int:
        return self._position

    @position.setter
    def position(self, value: int):
        self._position = int(value)

    @property
    def ref_base(self) -> str:
        return self._ref_base

    @ref_base.setter
    def ref_base(self, value: str):
        self._ref_base = value

    @property
    def alt_base(self) -> str:
        return self._alt_base

    @alt_base.setter
    def alt_base(self, value: str):
        self._alt_base = value

    @property
    def sensitivity(self) -> float:
        return self._sensitivity

    @sensitivity.setter
    def sensitivity(self, value: float):
        self._sensitivity = float(value)

    @property
    def specificity(self) -> float:
        return self._specificity

    @specificity.setter
    def specificity(self, value: float):
        self._specificity = float(value)

    @property
    def passes_filters(self) -> bool:
        return self._passes_filters

    @passes_filters.setter
    def passes_filters(self, value: bool):
        self._passes_filters = value

    @property
    def is_genotype_snp(self) -> bool:
        return self._is_genotype_snp

    @is_genotype_snp.setter
    def is_genotype_snp(self, value: bool):
        self._is_genotype_snp = value
        self._is_species_snp=False
    
    @property
    def is_species_snp(self) -> bool:
        return self._is_species_snp

    @is_species_snp.setter
    def is_species_snp(self, value: bool):
        self._is_species_snp = value
        self._is_genotype_snp=False

    def to_file(self, file_handle: TextIOWrapper, **kwargs):
        """Constructor
        :param file_handle: File handle to which the SNP will be written
        :type file_handle: TextIOWrapper

        :param sep: Column separator to use, defaults to tab
        :type sep: str, optional
        """        
        sep=kwargs.get("sep","\t")
        name=f'{self._ref_base}/{self.alt_base}'
        if self._is_species_snp:
            name=name+"/species"
        elif self._is_genotype_snp:
            name=name+"/genotype"

        file_handle.write(sep.join( [str(f) for f in [self._ref_contig_id, self.position, self.position+1, name] ] )+"\n")

    def __eq__(self, other) -> bool:
        return self.coordinate==other.coordinate and self.alt_base==other.alt_base

    def __hash__(self):
        return hash( (self.ref_contig_id, self.position, self.alt_base) )

    def copy(self):
        """Creates deep copy of the SNP instance
        """
        new_snp=copy.deepcopy(self)
        return new_snp

class Sample:
    """Represents a VCF sample and contains SNPs associated with it

    SNPs have to be loaded separately to avoid creating multitude of objects
    """
    def __init__(self, name: str, vcf_file: str) -> None:
        self._name=name
        self._snps=[]
        self._vcf_file=vcf_file
        self._uuid=str(uuid.uuid4())
        self._genotype=""

    
    @property
    def genotype(self) -> str:
        return self._genotype

    @genotype.setter
    def genotype(self, value: str):
        self._genotype = value

    @property
    def vcf_file(self) -> str:
        return self._vcf_file

    @property
    def snps(self) -> List[SNP]:
        return self._snps

    @property
    def name(self) -> str:
        return self._name
    
    @property
    def id(self) -> str:
        return self._uuid

class Amplicon:
    def __init__(self, name: str, seq: str) -> None:
        self._name: str=name
        self.seq: str=seq
        self._snps:List[SNP]=[]
        self._left_flanking_id=""
        self._right_flanking_id=""
        self._has_homologues=False
        self._uuid=str(uuid.uuid4())

    @classmethod
    def from_bed_line(cls, bed_line:str, ref_fasta_file: str):
        """Constructor using bedfile lines. 
        :param bed_line: String from bedfile, if the line has fourth column, this will be included in amplicon name
        :type bed_file: str

        :param ref_fasta_file: Path to fasta file on which the amplicon is based
        :type ref_fasta_file: str

        """
        bed_line_values = bed_line.strip().split("\t")
        ampl_chr, ampl_start, ampl_end=bed_line_values[0:3]
        if len(bed_line_values)>=4:
            name='_'.join( [str(f) for f in bed_line_values[0:4] ] )
        else:
            name='_'.join( [str(f) for f in bed_line_values[0:3] ] )
        ampl_start=int(ampl_start)
        ampl_end=int(ampl_end)
        for record in SeqIO.parse(ref_fasta_file,"fasta"):
            if record.id==ampl_chr:
                    new_amplicon=cls(name, str(record.seq[ampl_start:ampl_end]))
                    new_amplicon.ref_contig=record.id
                    new_amplicon.ref_start=ampl_start
                    new_amplicon.ref_end=ampl_end
                    return new_amplicon
        raise ValueError(f'Contig {ampl_chr} is not found in fasta file: {ref_fasta_file}')

    @property
    def left_flanking_id(self) -> str:
        return self._left_flanking_id

    @left_flanking_id.setter
    def left_flanking_id(self, value: str):
        self._left_flanking_id = value

    @property
    def right_flanking_id(self) -> str:
        return self._right_flanking_id

    @right_flanking_id.setter
    def right_flanking_id(self, value: str):
        self._right_flanking_id = value

    @property
    def has_flanking(self) -> bool:
        return self._left_flanking_id!="" and self._right_flanking_id!=""


    @property
    def name(self) -> str:
        return self._name

    @property
    def id(self) -> str:
        return self._uuid

    @property
    def len(self) -> int:
        return len(self.seq)

    @property
    def has_homologues(self) -> bool:
        return self._has_homologues

    @has_homologues.setter
    def has_homologues(self, value: bool):
        self._has_homologues = value

    @property
    def ref_contig(self) -> str:
        return self._ref_contig

    @ref_contig.setter
    def ref_contig(self, value: str):
        self._ref_contig = value

    @property
    def ref_start(self) -> int:
        return self._ref_start

    @ref_start.setter
    def ref_start(self, value: int):
        self._ref_start = int(value)

    @property
    def ref_end(self) -> int:
        return self._ref_end

    @ref_end.setter
    def ref_end(self, value: int):
        self._ref_end = int(value)

    def snp_in_amplicon(self, snp:SNP) -> bool:
        if snp.ref_contig_id==self.ref_contig and \
        snp.position>=self.ref_start and snp.position<=self.ref_end:
            return True
        else:
            return False

    @property
    def snps(self) -> List[SNP]:
        return self._snps

    @snps.setter
    def snps(self, value: List[SNP]):
        self._snps = value

    def __hash__(self):
        return hash(self.id)

class FlankingAmplicon(Amplicon):
    """Class for defining amplicons flanking a parent amplicon
    """

    def __init__(self, name: str, seq: str, parent: Amplicon, is_left: bool, max_len: int) -> None:
        super().__init__(name, seq)
        self._parent=parent
        if is_left:
            self._parent.left_flanking_id=self.id
        else:
            self._parent.right_flanking_id=self.id

        self._is_left=is_left
        self._max_len=max_len
    
    @property
    def parent(self) -> Amplicon:
        return self._parent

    @property
    def max_len(self) -> int:
        return self._max_len

    @max_len.setter
    def max_len(self, value: int):
        self._max_len = int(value)

    @property
    def is_left(self) -> bool:
        return self._is_left

    @is_left.setter
    def is_left(self, value: bool):
        self._is_left = value

    @classmethod
    def from_parent_bed_line(cls,ref_fasta_file: str, is_left: bool, max_len:int, parent: Amplicon):
        """Constructor basesd on parent sequence and maximum amplicon length. 
        :param ref_fasta_file: Path to fasta file on which the amplicon is based
        :type ref_fasta_file: str

        :param is_left: Boolean indicating if flank is left or right of parent amplicon
        :type is_left: str

        :param max_len: Maximum length of parent and this flanking sequence
        :type max_len: int        

        :param parent: The parent (i.e. actual or central) amplicon to which this one will be flanking
        :type parent: Amplicon
        """
        name = parent._name+"_left" if is_left else parent._name+"_right"
        for record in SeqIO.parse(ref_fasta_file,"fasta"):
            if record.id==parent.ref_contig:
                    new_amplicon=cls(name, "", parent, is_left, max_len )
                    ampl_start, ampl_end = new_amplicon._calculate_flanking_coordinates(record)
                    new_amplicon.seq=str(record.seq[ampl_start:ampl_end])
                    new_amplicon.ref_contig=record.id
                    new_amplicon.ref_start=ampl_start
                    new_amplicon.ref_end=ampl_end
                    return new_amplicon
        raise ValueError(f'Contig {parent.ref_contig} is not found in fasta file: {ref_fasta_file}')

    def _calculate_flanking_coordinates(self, record: SeqIO.SeqRecord) -> Tuple[int,int]:
        """Calculates the strat and end of the flanking amplicon sequences based on amplicon length
        and maximum permitted lenght of the amplicon
        """
        if self.is_left:
            return (max(self.parent.ref_start-self.max_len,0),self.parent.ref_start)
        else:
            return (self.parent.ref_end , min(self.parent.ref_end+self.max_len,len(record.seq)))

class Genotype:

    def __init__(self, name: str) -> None:
        self._name=name
        self._subgenotypes:List[str]=[name] #everygenotype has itself as subgenotypes
        self._defining_snps:List[SNP]=[]
        self._alleles: Dict[SNP, str]={}
        self._amplicons: List[Amplicon]=[]

    @property
    def name(self) -> str:
        return self._name
    
    @property
    def subgenotypes(self) -> List[str]:
        return self._subgenotypes

    @subgenotypes.setter
    def subgenotypes(self, value: List[str]):
        self._subgenotypes = value

    @property
    def amplicons(self) -> List[Amplicon]:
        return self._amplicons

    @amplicons.setter
    def amplicons(self, value: List[Amplicon]):
        self._amplicons = value

    @property
    def defining_snps(self) -> List[SNP]:
        return list(self._alleles.keys())

    def get_genotype_allele(self, snp: SNP) -> str:
        if snp not in self._alleles:
            raise ValueError(f'SNP with coordinates {snp.coordinate} is not present among snps of genotype {self._name}')
        return self._alleles[snp]

    def add_genotype_allele(self, snp: SNP, allele: str):
        self._alleles[snp]=allele

    @property
    def defining_snp_coordinates(self) -> List[Tuple[str,int]]:
        return [f.coordinate for f in self.defining_snps if f.passes_filters]
    
class BlastResult:
    def __init__(self) -> None:
        pass

    @property
    def q_hit_len(self) -> int:
        return len(self._qseq.replace("-",""))

    @property
    def qseq(self) -> str:
        return self._qseq

    @qseq.setter
    def qseq(self, value: str):
        self._qseq = value

    @property
    def evalue(self) -> float:
        return self._evalue

    @evalue.setter
    def evalue(self, value: float):
        self._evalue = float(value)

    @property
    def pident(self) -> float:
        return self._pident

    @pident.setter
    def pident(self, value: float):
        self._pident = float(value)

    @property
    def send(self) -> int:
        return self._send

    @send.setter
    def send(self, value: int):
        self._send = int(value)

    @property
    def sstart(self) -> int:
        return self._sstart

    @sstart.setter
    def sstart(self, value: int):
        self._sstart = int(value)

    @property
    def sseqid(self) -> str:
        return self._sseqid 

    @sseqid.setter
    def sseqid(self, value: str):
        self._sseqid = value

    @property
    def qseqid(self) -> str:
        return self._qseqid

    @qseqid.setter
    def qseqid(self, value: str):
        self._qseqid = value

    @property
    def qstart(self) -> int:
        return self._qstart

    @qstart.setter
    def qstart(self, value: int):
        self._qstart = int(value)

    @property
    def qend(self) -> int:
        return self._qend
    
    @qend.setter
    def qend(self, value: int):
        self._qend = int(value)

    @property
    def query_file_name(self) -> str:
        return self._query_file_name

    @query_file_name.setter
    def query_file_name(self, value: str):
        self._query_file_name = value

class Genotypes:
    def __init__(self, **kwargs) -> None:
        """Constructor

        :param genotypes: List of Genotype objects, defaults to empty list
        :type genotypes: List[Genotype], optional
        """
        self._genotypes=kwargs.get("genotypes",[])

    @property
    def genotypes(self) -> List[Genotype]:
        return self._genotypes

    @genotypes.setter
    def genotypes(self, value: List[Genotype]):
        self._genotypes = list(value)

    def all_snps_coord_sorted(self) -> List[Tuple[str, int]]:
        """Returns a list of unique contig + position pairs 
        sorted by contig and position
        """
        unique_contig_pos: List[Tuple[str, int]]=list(set([snp.coordinate for snps in self.genotypes for snp in snps.defining_snps]))
        unique_contig_pos=sorted(unique_contig_pos, key=lambda x: x)
        return unique_contig_pos

    def genotypes_to_snp_matrix(self) -> pd.DataFrame:
        """Converts list of genotypes into a matrix with each row a position on reference
        and each column a genotype. Two extra columns are Contig and Position
        The values indicate if the SNP passes some filter for that column (usually genotype)
        """
        if len(self._genotypes)==0:
            raise ValueError("The object has no genotypes in it.")
        unique_contig_pos: List[Tuple[str, int]]=self.all_snps_coord_sorted()
        gts: List[str]=[f.name for f in self.genotypes]
        output_df=pd.DataFrame(columns=["Contig","Position"]+gts, index=range(0,len(unique_contig_pos))).fillna(False)
        output_df["Contig"]=[f[0] for f in unique_contig_pos]
        output_df["Position"]=[f[1] for f in unique_contig_pos]
        for gt in self.genotypes:
            for snp in gt.defining_snps:
                index=unique_contig_pos.index((snp.ref_contig_id, snp.position))
                output_df.loc[index, gt.name]=snp.passes_filters
        return output_df
