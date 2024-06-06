from collections import Counter, OrderedDict
import re
from torchtext.vocab import vocab
from torch.utils.data.dataset import random_split

class MyEmbedding:
    def __init__(self, data:list):
        self.data = data
        self.train_dataset, self.valid_dataset = None, None
    
    
    def split(self):
        num_train = round(len(self.data)*.8)
        num_valid = len(self.data) - num_train
        self.train_dataset, self.valid_dataset = random_split(self.data, [num_train, num_valid])
        print(self.train_dataset[0])
        print(len(self.train_dataset), type(self.train_dataset))
        return self.train_dataset, self.valid_dataset

    def split2(self):
        num_train = round(len(self.data)*.8)
        num_valid = len(self.data) - num_train
        train_dataset, valid_dataset = random_split(self.data, [len(self.data)-1, 1])
        print(train_dataset[0])
        print(len(train_dataset), type(train_dataset))
        return train_dataset


    def tokenizer_single(self, text:str):
        return list(text)
    
    def tokenize(self, tokenizer=None):
        '''
        need self.train_dataset
        '''
        if tokenizer is None:
            tokenizer = self.tokenizer_single
        print("\n## Step 2 tokenization: unique tokens (words)...")
        # count tokens
        input_tokens, label_tokens = Counter(), Counter()
        for label, line in self.train_dataset:
            tokens = tokenizer(line)
            # words in list type
            input_tokens.update(tokens)
            label_tokens.update([label,])

        print('A sentence converted to tokens:', line, tokens)
        print('Vocab-size of input:', len(input_tokens))
        print('Vocab-size of labels:', len(label_tokens))
    
        # sort token couts of input
        sorted_by_freq_tuples = sorted(input_tokens.items(), key=lambda x: x[1], reverse=True)
        self.input_ordered_dict = OrderedDict(sorted_by_freq_tuples)
        print(self.input_ordered_dict)
        
        # sort token couts of labels
        sorted_by_freq_tuples = sorted(label_tokens.items(), key=lambda x: x[1], reverse=True)
        self.label_ordered_dict = OrderedDict(sorted_by_freq_tuples)
        counts = list(self.label_ordered_dict.values())
        print('counts of input:', counts)
        
        return self.input_ordered_dict, self.label_ordered_dict
    
    
    def build_vocab(self):
        '''
        need self.input_ordered_dict
        '''
        print("\n## Step 3 encoding: encoding each unique token into integers...")
        # Convert count value to index value (ranking)
        self.input_vocab = vocab(self.input_ordered_dict)
        self.input_vocab.insert_token("<pad>", 0)
        self.input_vocab.insert_token("<unk>", 1)
        # default token is "<unk>"
        self.input_vocab.set_default_index(1)
        
        # labels
        self.label_vocab = vocab(self.label_ordered_dict)
        # print(label_ordered_dict)
        self.label_vocab.insert_token("<pad>", 0)
        self.label_vocab.insert_token("<unk>", 1)
        # default token is "<unk>"
        self.label_vocab.set_default_index(1)
        
        return self.input_vocab, self.label_vocab