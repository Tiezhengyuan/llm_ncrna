import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score, matthews_corrcoef

from my_dataset import MyDataset

class MyProcess:
    loss_fn = nn.BCELoss()
    batch_size = 32

    def __init__(self, model, collate_fn, num_epochs=None):
        self.model = model
        self.optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        self.collate_fn = collate_fn
        self.num_epochs = 5 if num_epochs is None else num_epochs
    
    def run(self, train_dataset, valid_dataset):
        info = pd.DataFrame(np.zeros((self.num_epochs, 5), dtype=float), columns=\
            ['acc_train', 'loss_train', 'acc_valid', 'loss_valid', 'rec_valid'])

        ## Step 4: batching the datasets. shuffle data for each epoch
        train_dl = DataLoader(train_dataset, batch_size=self.batch_size, \
            shuffle=True, collate_fn=self.collate_fn)
        valid_dl = DataLoader(valid_dataset, batch_size=self.batch_size, \
            shuffle=False, collate_fn=self.collate_fn)

        for epoch in range(self.num_epochs):
            acc_train, loss_train = self.train(train_dl)
            acc_valid, loss_valid, rec_valid = self.evaluate(valid_dl)
            print(f'Epoch {epoch} acc: {acc_train:.4f} val_acc: {acc_valid:.4f}, val_recall: {rec_valid:.4f}')
            info.loc[epoch] = [acc_train, loss_train, acc_valid, loss_valid, rec_valid]
        return self.model, info

    def train(self, dataloader):
        self.model.train()
        total_acc, total_loss = 0, 0
        for text_batch, label_batch, lengths in dataloader:
            self.optimizer.zero_grad()
            pred = self.model(text_batch, lengths)[:, 0]
            loss = self.loss_fn(pred, label_batch)
            loss.backward()
            self.optimizer.step()
            total_acc += ((pred>=0.5).float() == label_batch).float().sum().item()
            total_loss += loss.item()*label_batch.size(0)
        return total_acc/len(dataloader.dataset), total_loss/len(dataloader.dataset)

    def evaluate(self, dataloader):
        self.model.eval()
        total_acc, total_loss = 0, 0
        np_target, np_pred = np.array([]), np.array([])
        with torch.no_grad():
            for text_batch, label_batch, lengths in dataloader:
                pred = self.model(text_batch, lengths)[:, 0]
                loss = self.loss_fn(pred, label_batch)
                total_acc += ((pred>=0.5).float() == label_batch).float().sum().item()
                total_loss += loss.item()*label_batch.size(0)
                
                # 
                _np_target = label_batch.view(-1).cpu().numpy()
                _np_pred = pred.view(-1).cpu().numpy()
                _np_pred = (_np_pred >= .5).astype(int)
                np_target = np.concatenate([np_target, _np_target])
                np_pred = np.concatenate([np_pred, _np_pred])
            recall = recall_score(np_target, np_pred, average = 'macro', zero_division=np.nan)
        return total_acc/len(dataloader.dataset), total_loss/len(dataloader.dataset), recall

    def test(self, datasets:dict):
        info = pd.DataFrame(np.zeros((len(datasets), 3), dtype=float), \
            columns=['acc', 'loss', 'recall'], index=list(datasets))

        for key, dataset in datasets.items():
            # dataset = MyDataset(dataset[0], dataset[1])
            dl = DataLoader(dataset, batch_size=self.batch_size,
                shuffle=False, collate_fn=self.collate_fn)
            acc, loss, rec = self.evaluate(dl)
            print(f'dataset: {key}, accuracy: {acc:.4f}, loss: {loss:.4f}, recall: {rec:.4f}')
            info.loc[key] = [acc, loss, rec]
        return info


    def predict(self, texts, labels):
        total_acc, total_pred = 0, np.array([])

        dataset = MyDataset(texts, labels)
        dataloader = DataLoader(dataset, batch_size=self.batch_size,
                shuffle=False, collate_fn=self.collate_fn)
        self.model.eval()
        with torch.no_grad():
            for text_batch, label_batch, lengths in dataloader:
                pred = self.model(text_batch, lengths)[:, 0]
                _pred = pred.view(-1).cpu().numpy()
                total_pred = np.concatenate([total_pred, _pred])
                total_acc += ((pred>=.5).float() == label_batch).float().sum().item()
        df = pd.DataFrame({'text':texts, 'predict':total_pred})
        return total_acc/len(dataloader.dataset), df
