import os
import torch
import numpy as np
import pandas as pd
import torch.nn as nn
import torch.nn.functional as F 
from torch_geometric.loader import DataLoader  
from torch.nn import Linear, Dropout, BatchNorm1d
from torch_geometric.data import Data, InMemoryDataset
from torch_geometric.nn import GCNConv, GINConv, SAGEConv, GATConv, TopKPooling, global_mean_pool, global_max_pool
from torch_geometric.nn import global_mean_pool, global_max_pool, global_add_pool
from ogb.utils import smiles2graph

torch.manual_seed(42)

embedding_size = 256
dropout_rate = 0.05
leaky_relu_slope = 0.01

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
                edge_index=torch.tensor(graph['edge_index'], dtype=torch.long),
                edge_attr=torch.tensor(graph['edge_feat'], dtype=torch.float)
            )
            data.smiles = row['SMILES']
            data_list.append(data)
        return data_list

# df = pd.read_csv('training_data_DPRA.csv')
# data_list = CustomMoleculeNetDataset_predict.create_data_list(df)
# dataset = CustomMoleculeNetDataset_predict(data_list)
# edge_feature_dim = dataset[0].edge_attr.size(1)

class MolecularGraphNeuralNetwork(nn.Module):
    def __init__(self, num_features, embedding_size, dropout_rate, leaky_relu_slope, edge_feature_dim):
        super(MolecularGraphNeuralNetwork, self).__init__()
        self.initial_conv = GCNConv(num_features, embedding_size)  # Use num_features from argument
        self.edge_embedding = nn.Linear(edge_feature_dim, embedding_size)  
        self.conv1 = GINConv(nn.Sequential(nn.Linear(embedding_size, embedding_size),
                                           nn.ReLU(),
                                           nn.Linear(embedding_size, embedding_size)))
        self.conv2 = SAGEConv(embedding_size, embedding_size)
        self.conv3 = GATConv(embedding_size, embedding_size, heads=4, concat=False)
        self.graph_embedding = nn.Linear(embedding_size + 1, embedding_size)  
        self.out = nn.Linear(embedding_size, 1)
        self.bn1 = nn.BatchNorm1d(embedding_size)
        self.bn2 = nn.BatchNorm1d(embedding_size)
        self.bn3 = nn.BatchNorm1d(embedding_size)
        self.dropout = nn.Dropout(dropout_rate)
        self.leaky_relu_slope = leaky_relu_slope

    def forward(self, x, edge_index, batch_index):
        x = self.initial_conv(x, edge_index)
        x = F.leaky_relu(x, negative_slope=self.leaky_relu_slope)
        x = self.bn1(x)
        x = self.dropout(x)
        x = self.conv1(x, edge_index)
        x = F.leaky_relu(x, negative_slope=self.leaky_relu_slope)
        x = self.bn2(x)
        x = self.dropout(x)
        x = self.conv2(x, edge_index)
        x = F.leaky_relu(x, negative_slope=self.leaky_relu_slope)
        x = self.bn3(x)
        x = self.dropout(x)
        x = self.conv3(x, edge_index)
        x = F.leaky_relu(x, negative_slope=self.leaky_relu_slope)
        x = global_mean_pool(x, batch_index)
        x = self.out(x)
        return x

        
def predict(df, model_choice):
    NUM_FOLDS = 5
    num_graphs_per_batch = 16
    test_data = CustomMoleculeNetDataset_predict.create_data_list(df)
    test_loader = DataLoader(test_data, batch_size=num_graphs_per_batch)

    model_path_map = {"in vitro (H-CLAT)": "skin_hCLAT",
                      "in vitro (KeratinoSens)": "skin_KeratinoSens",
                      "in vivo (LLNA)": "skin_LLNA",
                      "in chemico (DPRA)": "skin_DPRA",
                      "human": "skin_Human",}

    model_path = model_path_map.get(model_choice)
    if not model_path:
        raise ValueError(f"Invalid model choice: {model_choice}")

    models = []
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Initialize dataset and retrieve features
    data_list = CustomMoleculeNetDataset_predict.create_data_list(df)
    dataset = CustomMoleculeNetDataset_predict(data_list)
    num_features = dataset[0].x.size(1)  # Get number of node features
    edge_feature_dim = dataset[0].edge_attr.size(1)  # Get number of edge features

    for fold in range(NUM_FOLDS):
        model = MolecularGraphNeuralNetwork(
            num_features=num_features,
            embedding_size=embedding_size,
            dropout_rate=dropout_rate,
            leaky_relu_slope=leaky_relu_slope,
            edge_feature_dim=edge_feature_dim
        ).to(device)

        model_checkpoint_path = os.path.join(model_path, f'model_fold_{fold+1}.pth')
        checkpoint = torch.load(model_checkpoint_path, map_location=device)

        if 'module.' in list(checkpoint.keys())[0]:
            checkpoint = {k.replace('module.', ''): v for k, v in checkpoint.items()}

        model.load_state_dict(checkpoint)
        model.eval()
        models.append(model)

    predictions = []
    confidences = []

    for batch in test_loader:
        batch = batch.to(device)
        batch_predictions = []
        for model in models:
            model = model.to(device)
            with torch.no_grad():
                pred = model(batch.x.float().to(device), batch.edge_index.to(device), batch.batch.to(device))
                pred = torch.sigmoid(pred)
                batch_predictions.append(pred.cpu().numpy())

        batch_predictions = np.concatenate(batch_predictions, axis=1)
        mean_predictions = batch_predictions.mean(axis=1)
        prob_class_1 = mean_predictions
        prob_class_0 = 1 - prob_class_1
        softmax_predictions = np.stack([prob_class_0, prob_class_1], axis=1)
        max_probs = np.max(softmax_predictions, axis=1)
        binary_predictions = np.argmax(softmax_predictions, axis=1)
        predictions.extend(binary_predictions)
        confidences.extend(max_probs * 100)

    test_results = pd.DataFrame({'SMILES': df['SMILES'],
                                 'Prediction': predictions,
                                 'Confidence (%)': np.round(confidences, 2)})

    return test_results
