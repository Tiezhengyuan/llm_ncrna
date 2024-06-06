import torch
import torch.nn as nn

class RNN(nn.Module):
    
    def __init__(self, input_size, embed_dim, rnn_hidden_size, fc_hidden_size):
        super().__init__()
        self.embedding = nn.Embedding(input_size, embed_dim, padding_idx=0) 
        # layer: long-short term memory
        self.rnn = nn.LSTM(embed_dim, rnn_hidden_size, num_layers=2, batch_first=True)
        # fully-connected layer
        self.fc1 = nn.Linear(rnn_hidden_size, fc_hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(fc_hidden_size, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, text, lengths):
        out = self.embedding(text)
        # Packs a Tensor containing padded sequences of variable length.
        lengths = lengths.cpu().numpy()
        out = nn.utils.rnn.pack_padded_sequence(out, lengths, batch_first=True, enforce_sorted=False)
        # 
        out, (hidden, cell) = self.rnn(out)
        out = hidden[-1, :, :]
        out = self.fc1(out)
        out = self.relu(out)
        out = self.fc2(out)
        out = self.sigmoid(out)
        return out

   
