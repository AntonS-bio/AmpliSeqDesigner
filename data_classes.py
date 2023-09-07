class Amplicon:
    def __init__(self, id, seq) -> None:
        self.id: str=id
        self.seq: str=seq

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