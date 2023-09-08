from typing import List, Tuple
import pandas as pd
import uuid

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
        self._sensitivity = value

    @property
    def specificity(self) -> str:
        return self._specificity

    @specificity.setter
    def specificity(self, value: str):
        self._specificity = value

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

class Genotype:

    def __init__(self, name: str) -> None:
        self._name=name
        self._subgenotypes:List[str]=[name] #everygenotype has itself as subgenotypes
        self._defining_snps:List[SNP]=[]

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
    def defining_snps(self) -> List[SNP]:
        return self._defining_snps

    @defining_snps.setter
    def defining_snps(self, value: List[SNP]):
        self._defining_snps = value

    @property
    def defining_snp_coordinates(self) -> List[Tuple[str,int]]:
        return [f.coordinate for f in self.defining_snps if f.passes_filters]
    

class Amplicon:
    def __init__(self, name, seq) -> None:
        self._name: str=name
        self.seq: str=seq
        self._uuid=str(uuid.uuid4())

    @property
    def name(self) -> str:
        return self._name

    @property
    def id(self) -> str:
        return self._uuid


    @property
    def left_flank(self) -> bool:
        return self._left_flank

    @left_flank.setter
    def left_flank(self, value: bool):
        self._left_flank = value
        self._right_flank = False
        self._middle = False

    @property
    def right_flank(self) -> bool:
        return self._right_flank

    @right_flank.setter
    def right_flank(self, value: bool):
        self._right_flank = value
        self._left_flank = False
        self._middle = False

    @property
    def is_flanking(self) -> bool:
        return self._right_flank or self._left_flank

    @property
    def middle(self) -> bool:
        return self._middle

    @middle.setter
    def middle(self, value: bool):
        self._middle = value
        self._right_flank = False
        self._left_flank = False

    @property
    def len(self) -> int:
        return len(self.seq)

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
    def position(self, value: List[Genotype]):
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
