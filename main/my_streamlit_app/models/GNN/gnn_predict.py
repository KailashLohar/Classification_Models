import os
import torch
import numpy as np
import pandas as pd
import torch.nn.functional as F 
from torch_geometric.loader import DataLoader  
from torch.nn import Linear, Dropout, BatchNorm1d
from torch_geometric.data import Data, InMemoryDataset
from torch_geometric.nn import GCNConv, SAGEConv, GraphConv
from torch_geometric.nn import global_mean_pool, global_max_pool, global_add_pool
from ogb.utils import smiles2graph

torch.manual_seed(42)

class MolecularGraphNeuralNetwork(torch.nn.Module):
    def __init__(self):
        super(MolecularGraphNeuralNetwork, self).__init__()
        embedding_size = 128  
        self.initial_conv = GCNConv(9, embedding_size)
        self.conv1 = GCNConv(embedding_size, embedding_size)
        self.conv2 = GCNConv(embedding_size, embedding_size)
        self.conv3 = GCNConv(embedding_size, embedding_size)
        self.out = torch.nn.Linear(embedding_size, 1)
        self.bn1 = torch.nn.BatchNorm1d(embedding_size)
        self.bn2 = torch.nn.BatchNorm1d(embedding_size)
        self.bn3 = torch.nn.BatchNorm1d(embedding_size)
        self.dropout = torch.nn.Dropout(0.2)

    def forward(self, x, edge_index, batch_index):
        x = self.initial_conv(x, edge_index)
        x = F.leaky_relu(x, negative_slope=0.01)
        x = self.bn1(x)
        x = self.dropout(x)
        x = self.conv1(x, edge_index)
        x = F.leaky_relu(x, negative_slope=0.01)
        x = self.bn2(x)
        x = self.dropout(x)
        x = self.conv2(x, edge_index)
        x = F.leaky_relu(x, negative_slope=0.01)
        x = self.bn3(x)
        x = self.dropout(x)
        x = self.conv3(x, edge_index)
        x = F.leaky_relu(x, negative_slope=0.01)
        x = global_mean_pool(x, batch_index)
        x = self.out(x)
        return x

class CustomMoleculeNetDataset_predict(InMemoryDataset):
    def __init__(self, data_list):
        super(CustomMoleculeNetDataset_predict, self).__init__(".", transform=None, pre_transform=None)
        self.data_list = data_list
        self.data, self.slices = self.collate(data_list)

    @staticmethod
    def create_data_list(df):
        data_list = []
        for _, row in df.iterrows():
            graph = smiles2graph(row['SMILES'])
            data = Data(
                x=torch.tensor(graph['node_feat']),
                edge_index=torch.tensor(graph['edge_index']),
                edge_attr=torch.tensor(graph['edge_feat'])
            )
            data.smiles = row['SMILES']
            data_list.append(data)
        return data_list

def predict_gnn(df):
    NUM_FOLDS = 5
    num_graphs_per_batch = 16
    test_data = CustomMoleculeNetDataset_predict.create_data_list(df)
    test_loader = DataLoader(test_data, batch_size=num_graphs_per_batch)

    models = []
    for fold in range(NUM_FOLDS):
        model = MolecularGraphNeuralNetwork()
        model_checkpoint_path = os.path.join('models', 'GNN', 'model_files', f'model_fold_{fold+1}.pth')
        checkpoint = torch.load(model_checkpoint_path, map_location=torch.device('cpu')) 

        if 'module.' in list(checkpoint.keys())[0]:
            checkpoint = {k.replace('module.', ''): v for k, v in checkpoint.items()}

        model.load_state_dict(checkpoint)  
        model.eval()
        models.append(model)

    predictions = []
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    for batch in test_loader:
        batch = batch.to(device)
        batch_predictions = []
        for model in models:
            model = model.to(device)  
            with torch.no_grad():
                pred = model(batch.x.float().to(device), batch.edge_index.to(device), batch.batch.to(device))
                batch_predictions.append(pred.cpu().numpy())

        batch_predictions = np.concatenate(batch_predictions, axis=1)
        mean_predictions = batch_predictions.mean(axis=1)
        mean_predictions = (mean_predictions > 0.5).astype(int)
        predictions.extend(mean_predictions)

    test_results = pd.DataFrame({'SMILES': df['SMILES'], 'GNN_Prediction': predictions})
    return test_results
