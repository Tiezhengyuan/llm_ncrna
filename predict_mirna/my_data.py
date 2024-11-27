# retrieve sequences and species as label
from typing import Iterable
import random
import re


from my_dataset import MyDataset
from my_parser import MyParser

class MyData:
    def __init__(self):
        max_len = 30
        self.parser = MyParser(max_len)
        self.train_data, self.species = [], {}
        self.len_range = range(18, max_len+1)

    def get_full_data(self) -> list:
        # true miRNAs
        self.wrap(self.origin_data(), True)
        # trimmed seq
        self.wrap(self.trim_data(), False)
        # mutation
        self.wrap(self.mutation_data(1), False)

        #Negative: shuffle data
        self.wrap(self.shuffle_data(1), False)
        #Negative: other RNAs
        others = ['5srrnadb.fasta', 'gtrnadb.fasta', 'pirbase.fasta']
        for rna_file in others:
            self.wrap(self.parser.other_rna(rna_file), False)


        return self.train_data        

    def wrap(self, data_tuple, new=False) -> list:
        '''
        return data for train
        '''
        if new is True:
            self.train_data = []
        data = []
        n=0
        texts, labels = data_tuple
        for _text, _label in zip(texts, labels):
            # if re.findall(r'[^A|T|G|C|N]', _text):
            #     print('###########',n, _label, _text)
            n+=1
            data.append((_label, _text))
        self.train_data += data
        return data


    def shuffle_data(self, repeat:int):
        '''
        label='shuffle'
        '''
        n = 0
        labels, texts = [], []
        for seq, specie in self.parser.iterate_fa():
            for i in range(repeat):
                seq2 = list(seq)
                random.shuffle(seq2)
                texts.append(''.join(seq2))
                n += 1
        labels = ['shuffle',] * len(texts)
        for l, seq in zip(labels, texts):
            if re.findall(r'[^A|T|G|C]', seq):
                print('###########', l, seq)
        print(f"Number of shuffle seq: {n}")
        return texts, labels

    def trim_data(self):
        n = 0
        labels, texts = [], []
        for sub_len in self.len_range:
            for seq, specie in self.parser.iterate_fa():
                len_seq = len(seq)
                if len_seq > sub_len:            
                    for i in range(0, len_seq-sub_len+1):
                        labels.append(specie)
                        texts.append(seq[i:i+sub_len])
                        n += 1
        print(f"Number of trim seq: {n}")
        return texts, labels

    def random_data(self) -> list:
        '''
        label='random'
        '''
        n, repeat = 0, 1000
        labels, texts = [], []
        for sub_len in self.len_range:
            for i in range(repeat):
                rand_seq = [random.choice(list('ATGC')) for i in range(sub_len)]
                labels.append('random')
                texts.append(''.join(rand_seq))
            n += repeat
        print(f"Number of random seq: {n}")
        return texts, labels

    def mutation_data(self, mut:int):
        texts, labels = [], []
        for seq, specie in self.parser.iterate_fa():
            for i in range(mut):
                pos = random.choice(range(len(seq)))
                seq = seq[:pos] + 'N' + seq[pos+1:]
            texts.append(seq)
            labels.append('mutation')
        return texts, labels


#
    def get_origin_data(self) -> list:
        """
        update self.train_data including all true miRNAs
        """
        self.wrap(self.origin_data(), True)
        return self.train_data
    
    def origin_data(self) -> tuple:
        texts, labels, n = [], [], 0
        for seq, specie in self.parser.iterate_fa():
            texts.append(seq)
            labels.append(specie)
            n += 1
        print(f"Number of origin seq: {n}")
        return texts, labels

    def get_data(self) -> list:
        n = 0
        self.train_data, self.species = [], {}
        for seq, specie in self.mirna_iter:
            self.train_data.append((specie, seq))
            if specie not in self.species:
                self.species[specie] = 0
            self.species[specie] += 1
            n += 1
        print(f"Number of miRNA seq: {n}")
        return self.train_data
    
    def get_length_data(self):
        data = {}
        for seq, specie in self.parser.iterate_fa():
            _len = len(seq)
            if _len not in data:
                data[_len] = {'texts': [], 'labels': []}
            data[_len]['texts'].append(seq)
            data[_len]['labels'].append(specie)
        
        for _len in data:
            _data = data[_len]
            data[_len] = MyDataset(_data['texts'], _data['labels'])
        return data

    def get_random_data(self):
        '''
        label='random'
        '''
        rand_dataset, n, repeat = {}, 0, 10**5
        for sub_len in range(20, 29):
            text = []
            for i in range(repeat):
                rand_seq = [random.choice(list('ATGC')) for i in range(sub_len)]
                text.append(''.join(rand_seq))
            n += repeat
            if text:
                labels = ['random',] * len(text)
                rand_dataset[sub_len] = MyDataset(text, labels)
        print(f"Number of trim seq: {n}")
        return rand_dataset

    def get_trim_data(self):
        '''
        label='trim'
        '''
        rand_dataset, n = {}, 0
        for sub_len in range(18, 29):
            texts, labels = [], []
            for seq, specie in self.parser.iterate_fa():
                len_seq = len(seq)
                if len_seq > sub_len:            
                    for i in range(0, len_seq-sub_len+1):
                        labels.append(specie)
                        texts.append(seq[i:i+sub_len])
                        n += 1
            rand_dataset[sub_len] = MyDataset(texts, labels)
        print(f"Number of trim seq: {n}")
        return rand_dataset

    def get_mutation_data(self, mut:int):
        text, labels = [], []
        for seq, specie in self.parser.iterate_fa():
            if mut > 0:
                for i in range(0, mut):
                    pos = random.choice(range(len(seq)))
                    seq = seq[:pos] + 'A' + seq[pos+1:]
                text.append(seq)
                labels.append('mutation')
            else:
                text.append(seq)
                labels.append(specie)
        return MyDataset(text, labels)

    def top_species(self, num_top=3):
        ordered_species = sorted(self.species.items(), key=lambda x: x[1], reverse=True)
        top_species= []
        for k,v in ordered_species:
            print('\t', k,v)
            top_species.append(k)
            if len(top_species) >= num_top:
                break
        return top_species