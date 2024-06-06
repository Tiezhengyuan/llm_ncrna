from Bio import SeqIO
from typing import Iterable



class MyParser:

    def __init__(self, max_len:int = None):
        self.max_len = max_len if max_len and max_len>=18 else None

    def kmer(self, seq:str) -> list:
        if len(seq) >= self.max_len:
            return [seq, ]
        res = []
        for i in range(0, len(seq)-self.max_len+1):
            res.append(seq[i:i+self.max_len])
        return res

    def iterate_fa(self) -> Iterable:
        fa_file = 'mature.fa'
        with open(fa_file, 'r') as f:
            parser = SeqIO.parse(f, 'fasta')
            if self.max_len:
                for rec in parser:
                    seq = str(rec.seq).replace('U', 'T')
                    if len(seq) <= self.max_len:
                        specie = ' '.join(rec.description.split(' ')[2:4])
                        yield seq, specie
            else:
                for rec in parser:
                    seq = str(rec.seq).replace('U', 'T')
                    specie = ' '.join(rec.description.split(' ')[2:4])
                    yield seq, specie

    def other_rna(self, fa_file:str):
        labels, texts, n = [], [], 0
        with open(fa_file, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                _seq = str(rec.seq).replace('U', 'T')
                kmer_seq = self.kmer(_seq)
                texts += kmer_seq
                labels += ['other' for i in kmer_seq]
                n += 1
        print(f"Number of {fa_file}: {n}")
        return texts, labels

    def get_192_test(self, num_rec:int=None):
        fa_file = 'unaligned_192.fa'
        data, n = [], 0
        with open(fa_file, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                _seq = str(rec.seq)
                if len(_seq)>=20:
                    data.append(('192',_seq))
                    n += 1
                if num_rec and n >= num_rec:
                    break
        return data

    def get_192(self, num_rec:int=None):
        fa_file = 'unaligned_192.fa'
        labels, texts, n = [], [], 0
        with open(fa_file, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                _seq = str(rec.seq)
                if len(_seq)>=20:
                    texts.append(_seq)
                    labels.append('192')
                    n += 1
                if num_rec and n >= num_rec:
                    break
        return texts, labels