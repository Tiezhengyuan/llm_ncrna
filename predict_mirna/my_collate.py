import torch
import torch.nn as nn

class MyCollate:

    def __init__(self, text_pipeline, label_pipeline):
        self.text_pipeline = text_pipeline
        self.label_pipeline = label_pipeline
    
    def __call__(self, batch):
    ## Step 3-B: wrap the encode and transformation function
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        label_list, text_list, lengths = [], [], []
        for _label, _text in batch:
            label_list.append(self.label_pipeline(_label))
            # print(_label, _text)
            processed_text = torch.tensor(self.text_pipeline(_text), dtype=torch.int64)
            text_list.append(processed_text)
            lengths.append(processed_text.size(0))
        label_list = torch.tensor(label_list, dtype=torch.float32)
        lengths = torch.tensor(lengths)
        # padding is appended to the end of token vector of each sentence
        padded_text_list = nn.utils.rnn.pad_sequence(text_list, batch_first=True)
        # print(label_list, text_list)
        return padded_text_list.to(device), label_list.to(device), lengths.to(device)

    # slice sentence by word.
    def tokenizer_single(self, text:str):
        return list(text)
    
    def tokenizer_kmer(self, text:str):
        step, res = 6, []
        for i in range(0, len(text)-step+1):
            res.append(text[i:i+step])
        return res
