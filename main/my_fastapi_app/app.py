import torch
import mlflow
import numpy as np
import pandas as pd

from ogb.utils import smiles2graph
from torch.utils.data import DataLoader
from fastapi.responses import HTMLResponse
from torch_geometric.data import Data, Batch
from fastapi.templating import Jinja2Templates
from torch_geometric.data import InMemoryDataset
from fastapi import FastAPI, HTTPException, Request, Form

app = FastAPI()

model_uri = '/mlflow/mlruns/900649188862028013/33a1cfef715648b4909e2abdff31b24f/artifacts/best_model_fold_4'
model = mlflow.pytorch.load_model(model_uri)
model.eval()
model = model.to("cpu")


templates = Jinja2Templates(directory="templates")

class CustomMoleculeNetDatasetPredict(InMemoryDataset):
    def __init__(self, data_list):
        super().__init__(".", transform=None, pre_transform=None)
        self.data_list = data_list
        self.data, self.slices = self.collate(data_list)

    @staticmethod
    def create_data_list(df):
        data_list = []
        for _, row in df.iterrows():
            graph = smiles2graph(row['SMILES'])
            data = Data(
                x=torch.tensor(graph['node_feat'], dtype=torch.float32),
                edge_index=torch.tensor(graph['edge_index'], dtype=torch.long),
                edge_attr=torch.tensor(graph['edge_feat'], dtype=torch.float32)
            )
            data.smiles = row['SMILES']
            data_list.append(data)
        return data_list

@app.get("/", response_class=HTMLResponse)
def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.post("/predict", response_class=HTMLResponse)
def predict(request: Request, smiles: str = Form(...)):
    try:
        smiles_data = pd.DataFrame([{'SMILES': smiles}])
        test_data = CustomMoleculeNetDatasetPredict.create_data_list(smiles_data)

        def collate_fn(data_list):
            return Batch.from_data_list(data_list)

        test_loader = DataLoader(test_data, batch_size=1, collate_fn=collate_fn)

        optimal_threshold = 0.5
        predictions = []
        for batch in test_loader:
            batch = batch.to("cpu")
            with torch.no_grad():
                pred = model(batch.x, batch.edge_index, batch.batch)
                pred_prob = torch.sigmoid(pred).cpu().numpy()
                batch_predictions = (pred_prob > optimal_threshold).astype(int).flatten()
                predictions.extend(batch_predictions)

        return templates.TemplateResponse("result.html", {
            "request": request,
            "smiles": smiles,
            "prediction": int(predictions[0])
        })

    except Exception as e:
        return templates.TemplateResponse("error.html", {
            "request": request,
            "error_message": str(e)
        })
