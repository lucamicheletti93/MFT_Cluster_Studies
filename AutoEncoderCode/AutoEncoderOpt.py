import shap
import torch
from torchvision import datasets
from torchvision import transforms
import matplotlib.pyplot as plt
import ROOT
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data
from torch.utils.data import DataLoader, TensorDataset
import uproot  # This is to read the root files
import optuna
from optuna.visualization import plot_optimization_history


root_file = ROOT.TFile("RootFiles/new/Input.root")
tree = root_file.Get("input")

branches = [ "PdgID", "AvgClSize", "Eta", "ClSizeLayer0", "ClSizeLayer1", "ClSizeLayer2", "ClSizeLayer3", "ClSizeLayer4", "ClSizeLayer5", "ClSizeLayer6", "ClSizeLayer7", "PattIDLayer0", "PattIDLayer1", "PattIDLayer2", "PattIDLayer3", "PattIDLayer4", "PattIDLayer5", "PattIDLayer6", "PattIDLayer7"]  

# Preprocessing
data_dict = {}
for branch in branches:
    data_dict[branch] = []
    for entry in tree:
        data_dict[branch] = np.array(getattr(tree, branch))

df = pd.DataFrame.from_dict(data_dict)


df = df.drop(df[(df["PdgID"] != 111) & (df["PdgID"]!=211) & (df["PdgID"]!=130) & (df["PdgID"]!=310) & (df["PdgID"]!=311) & (df["PdgID"]!=321) & (df["PdgID"]!=2212) & (df["PdgID"]!=13) & (df["PdgID"]!=1000020030) ].index)

def sortPdgID(PdgID):
    if (PdgID==111) | (PdgID==211) | (PdgID==130) | (PdgID==310) | (PdgID==311) | (PdgID==321) | (PdgID==2212) | (PdgID==13):
        PdgID=1
    elif PdgID==1000020030:
        PdgID=2
    return PdgID

df['PdgID'] = df['PdgID'].apply(sortPdgID)


# Preparation for Training
x_train, x_test, y_train, y_test = train_test_split(df.values, df.values[:,0:1], test_size=0.2, random_state=111)

scaler = MinMaxScaler()
data_scaled = scaler.fit(x_train)
train_data_scaled = data_scaled.transform(x_train)
test_data_scaled = data_scaled.transform(x_test)



normal_train_data = pd.DataFrame(train_data_scaled).add_prefix('c').query('c0 == 0').values[:,1:]
anomaly_train_data = pd.DataFrame(train_data_scaled).add_prefix('c').query('c0 > 0').values[:, 1:]
normal_test_data = pd.DataFrame(test_data_scaled).add_prefix('c').query('c0 == 0').values[:,1:]
anomaly_test_data = pd.DataFrame(test_data_scaled).add_prefix('c').query('c0 > 0').values[:, 1:]

normal_val_data = normal_train_data[-int(0.2*len(normal_train_data)):]
anomaly_val_data = anomaly_train_data[-int(0.2*len(anomaly_train_data)):]

data_tensor = torch.FloatTensor(normal_train_data).unsqueeze(1) 
data_test_tensor =  torch.FloatTensor(normal_test_data).unsqueeze(1) 
dim = len(normal_train_data[0])

validation_percentage = 0.2

num_validation = int(validation_percentage * len(data_tensor))

# Split the data into training and validation sets
train_data = data_tensor[:-num_validation]
val_data = data_tensor[-num_validation:]



class AE(torch.nn.Module):
    def __init__(self, input_dim, layers, nodes_per_layer, activation):
        super().__init__()
# Define the encoder layers
        encoder_layers = []
        prev_dim = input_dim
        for i in range(layers):
            encoder_layers.append(nn.Linear(prev_dim, nodes_per_layer[i]))
            encoder_layers.append(activation)
            prev_dim = nodes_per_layer[i]

        self.encoder = nn.Sequential(*encoder_layers)

        # Define the decoder layers
        decoder_layers = []
        for i in range(layers-2, -1, -1):
            decoder_layers.append(nn.Linear(prev_dim, nodes_per_layer[i]))
            decoder_layers.append(activation)
            prev_dim = nodes_per_layer[i]

        decoder_layers.append(nn.Linear(prev_dim, input_dim))
        self.decoder = nn.Sequential(*decoder_layers)


    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

def objective(trial):

    n_layers = trial.suggest_int('n_layers', 3, 7)
    nodes_per_layer = [trial.suggest_int(f'nodes_l{i}', 2, 20) for i in range(n_layers)]
    lr = trial.suggest_float("lr", 1e-5, 1e-1, log=True)

    model = AE(dim, n_layers, nodes_per_layer, nn.ReLU())


    optimizer = torch.optim.Adam(model.parameters(),
            lr = lr,
            weight_decay = 1e-8)
    batch_size = 100

    criterion = nn.L1Loss()
    train_loader = data.DataLoader(train_data, batch_size=batch_size, shuffle=True)
    val_loader = data.DataLoader(val_data, batch_size=batch_size)




    device = 'cuda' if torch.cuda.is_available() else 'cpu'


    # Evaluate the model on the entire dataset to get the reconstruction errors
    model.eval()
    loss_dist_anomal = []
    loss_dist = []
    for i in range(len(anomaly_val_data)):
        data_val = torch.from_numpy(anomaly_val_data[i]).float()
        sample_val = model(data_val.to(device))
        loss = criterion(data_val.to(device), sample_val)
        loss_dist_anomal.append(loss.item())


    for i in range(len(normal_val_data)):
        data_val = torch.from_numpy(normal_val_data[i]).float()
        sample_val = model(data_val.to(device))
        loss = criterion(data_val.to(device), sample_val)
        loss_dist.append(loss.item())


    loss_dist = np.array(loss_dist)
    loss_dist_anomal = np.array(loss_dist_anomal)



    return loss_dist_anomal.mean()-loss_dist.mean() 


study = optuna.create_study(direction="maximize", study_name="example-study", storage="sqlite:///example.db")
study.optimize(objective, n_trials=100)

trial = study.best_trial
print("Best trial:")
print("  Value: ", trial.value)
print("  Params: ")
for key, value in trial.params.items():
    print("    {}: {}".format(key, value))

fig = plot_optimization_history(study)

fig.write_image("Plots/new/AutoEncoder/optimization_history_plot.png")

fig = optuna.visualization.plot_parallel_coordinate(study, params=["n_layers", "lr", "nodes_l0"])


fig.write_image("Plots/new/AutoEncoder/parallel_coordinate.png")

fig = optuna.visualization.plot_param_importances(study)

fig.write_image("Plots/new/AutoEncoder/hyperparameter_importance.png")





