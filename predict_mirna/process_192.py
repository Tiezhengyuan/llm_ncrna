import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

class Process:

    def __init__(self):
        self.samples = []
        self.min_rc = 3

    def __call__(self):
        self.iter_seq = self.iter_files()
        self.count()
        self.to_fasta()
        self.to_df()

    def iter_files(self):
        min_len = 18
        indir = '/home/yuan/data/192_results'
        for name in os.listdir(indir):
            sample_name = name.split('_')[0]
            if sample_name not in self.samples:
                self.samples.append(sample_name)
            infile = os.path.join(indir, name)
            print(f"scan {sample_name}...")
            with open(infile, 'r') as f:
                parser = SeqIO.parse(f, 'fasta')
                for rec in parser:
                    if len(rec.seq) >= min_len:
                        rc = rec.description.split('_')[-1]
                        yield sample_name, str(rec.seq), rc
                # break

    def count(self):
        self.counts, self.rc = {}, {}
        for sample_name, seq, rc in self.iter_seq:
            if int(rc) >= self.min_rc:
                if seq not in self.counts:
                    self.counts[seq] = []
                self.counts[seq].append({
                    'sample': sample_name,
                    'counts': rc,
                })
                if seq not in self.rc:
                    self.rc[seq] = {}
                self.rc[seq][sample_name] = rc
        print(f"total of sequences: {len(self.counts)}")
    
    def to_fasta(self):
        outfile = '/home/yuan/bio/llm_ncrna/predict_mirna/unaligned_192.fa'
        with open(outfile, 'w') as f:
            n = 1
            for seq, v in self.counts.items():
                sequence = SeqRecord(Seq(seq), id=str(n), description='')
                SeqIO.write(sequence, f, 'fasta')
                n += 1

    def to_df(self):
        self.samples.sort()
        
        outfile = '/home/yuan/bio/llm_ncrna/predict_mirna/unaligned_192.csv'
        with open(outfile, 'w') as f:
            f.write(','.join(['sequence',]+self.samples+['\n',]))
            for seq, v in self.rc.items():
                item = [seq,] + [v.get(str(i), '0') for i in self.samples] + ['\n',]
                f.write(','.join(item))



if __name__ == "__main__":
    Process()()