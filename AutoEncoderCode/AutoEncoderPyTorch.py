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
from torch.utils.data import DataLoader, TensorDataset, random_split
import uproot  # This is to read the root files


root_file = ROOT.TFile("RootFiles/new/Input.root")
tree = root_file.Get("input")

branches = [ "PdgID", "AvgClSize", "Eta", "ClSizeLayer0", "ClSizeLayer1", "ClSizeLayer2", "ClSizeLayer3", "ClSizeLayer4", "ClSizeLayer5", "ClSizeLayer6", "ClSizeLayer7",  "PattIDLayer0", "PattIDLayer1", "PattIDLayer2", "PattIDLayer3", "PattIDLayer4", "PattIDLayer5", "PattIDLayer6", "PattIDLayer7"]  

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

data_tensor = torch.FloatTensor(normal_train_data).unsqueeze(1) 
anomaly_tensor = torch.FloatTensor(anomaly_train_data[:100]).unsqueeze(1)
data_test_tensor =  torch.FloatTensor(normal_test_data).unsqueeze(1) 
dim = len(normal_train_data[0])

# validation split
train_data, val_data = random_split(data_tensor, [round(0.8*len(data_tensor)), round(0.2*len(data_tensor))])
#val_data+=anomaly_tensor
layers = 4
nodes_per_layers = [4,9,2,8]
lr = 0.0003


class AE(torch.nn.Module):
    def __init__(self, input_dim, layers, nodes_per_layer):
        super().__init__()
# Define the encoder layers
        encoder_layers = []
        prev_dim = input_dim
        for i in range(layers):
            encoder_layers.append(nn.Linear(prev_dim, nodes_per_layer[i]))
            encoder_layers.append(nn.ReLU())
            prev_dim = nodes_per_layer[i]

        self.encoder = nn.Sequential(*encoder_layers)

        # Define the decoder layers
        decoder_layers = []
        for i in range(layers-2, -1, -1):
            decoder_layers.append(nn.Linear(prev_dim, nodes_per_layer[i]))
            decoder_layers.append(nn.ReLU())
            prev_dim = nodes_per_layer[i]

        decoder_layers.append(nn.Linear(prev_dim, input_dim))
        self.decoder = nn.Sequential(*decoder_layers)


    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

model = AE(dim, layers, nodes_per_layers)


optimizer = torch.optim.Adam(model.parameters(),
        lr = lr,
        weight_decay = 1e-8)
batch_size = 100
criterion = nn.L1Loss()
train_loader = data.DataLoader(train_data, batch_size=batch_size, shuffle=True)
val_loader = data.DataLoader(val_data, batch_size=batch_size, shuffle=True)

train_losses = []
val_losses = []



num_epochs = 50
for epoch in range(num_epochs):
    model.train()
    running_loss = 0.0
    for data_batch in train_loader:
        optimizer.zero_grad()
        reconstructions = model(data_batch)
        loss = criterion(reconstructions, data_batch)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()
    epoch_train_loss = running_loss / len(train_loader)
    train_losses.append(epoch_train_loss)
    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {epoch_train_loss:.4f}")

    model.eval()  # Set the model in evaluation mode
    running_val_loss = 0.0

    with torch.no_grad():
        for val_batch in val_loader:
            val_reconstructions = model(val_batch)
            val_loss = criterion(val_reconstructions, val_batch)
            running_val_loss += val_loss.item()

    # Compute the average validation loss for this epoch
    epoch_val_loss = running_val_loss / len(val_loader)
    val_losses.append(epoch_val_loss)
    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {epoch_val_loss:.4f}")

device = 'cuda' if torch.cuda.is_available() else 'cpu'

# Evaluate the model on the entire dataset to get the reconstruction errors
model.eval()
loss_dist_anomal = []
loss_dist = []
for i in range(len(anomaly_test_data)):
    data_test = torch.from_numpy(anomaly_test_data[i]).float()
    sample_test = model(data_test.to(device))
    loss = criterion(data_test.to(device), sample_test)
    loss_dist_anomal.append(loss.item())


for i in range(len(normal_test_data)):
    data_test = torch.from_numpy(normal_test_data[i]).float()
    sample_test = model(data_test.to(device))
    loss = criterion(data_test.to(device), sample_test)
    loss_dist.append(loss.item())


loss_dist = np.array(loss_dist)
loss_dist_anomal = np.array(loss_dist_anomal)

azureblue = (85./255, 170./255, 216./255)
darkyellow=(0.800000011920929, 0.800000011920929, 0.20000000298023224)

fig, ax = plt.subplots()
plt.plot(train_losses[1:], label='Training Loss', color=azureblue)
plt.plot(val_losses[:-1], label='Validation Loss', color='red')
plt.xlabel('Epoch')
plt.ylabel('Loss')
ax.xaxis.set_label_coords(0.85,-0.1)
plt.legend(frameon=False)
plt.text(10,0.20, "Run 3 enriched MC", fontsize=14)
plt.text(10,0.18, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig('Plots/new/AutoEncoder/Loss.png')
plt.close()


normal_train_data = normal_train_data[:100]
data_shap = torch.tensor(normal_train_data, dtype=torch.float32)
explainer = shap.DeepExplainer(model, data_shap)

# Calculate SHAP values for the test set
shap_values = explainer.shap_values(data_shap)

feature_names = [branches[i] for i in range(1,len(normal_train_data[0])+1)]
shap.summary_plot(shap_values[0], data_shap, plot_type="bar", show=False, feature_names=feature_names)

plt.savefig('Plots/new/AutoEncoder/shap_summary_plot.png')
plt.close()




def CalculateEffPurAE(threshold, loss_dist, loss_dist_anomal):

    normal_eff = (loss_dist < threshold).sum()/len(loss_dist)
    anomal_eff = (loss_dist_anomal > threshold).sum()/len(loss_dist_anomal)
    below_pur = (loss_dist < threshold).sum()/((loss_dist < threshold).sum()+(loss_dist_anomal < threshold).sum())
    above_pur = (loss_dist_anomal > threshold).sum()/((loss_dist > threshold).sum()+(loss_dist_anomal > threshold).sum())

    anomal_sig = (loss_dist_anomal > threshold).sum()/np.sqrt((loss_dist > threshold).sum()+(loss_dist_anomal > threshold).sum())

    return (anomal_eff, above_pur, normal_eff, below_pur, anomal_sig)


def CalculateEffPur(AvgClThreshold):


    normal_df = df[(df["PdgID"]==1)]
    anomal_df = df[(df["PdgID"]==2)]
    below_df = df[(df["AvgClSize"]<AvgClThreshold)]
    above_df = df[(df["AvgClSize"]>AvgClThreshold)]

    normal_eff = len(normal_df[normal_df["AvgClSize"] < AvgClThreshold])/len(normal_df)
    anomal_eff = len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold])/len(anomal_df)
    below_pur = len(below_df[below_df["PdgID"] == 1])/len(below_df)
    above_pur = len(above_df[above_df["PdgID"] == 2])/len(above_df)
    anomal_sig = len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold])/np.sqrt(len(normal_df[normal_df["AvgClSize"] > AvgClThreshold])+len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold]))
    return (anomal_eff, above_pur, normal_eff, below_pur, anomal_sig)



thresholdsAE = np.linspace(0.1,0.3,20)
anomal_effAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[0] for x in thresholdsAE] 
above_purAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[1] for x in thresholdsAE]
normal_effAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[2] for x in thresholdsAE]
below_purAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[3] for x in thresholdsAE]
anomal_sigAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[4] for x in thresholdsAE]


thresholds = np.linspace(4.5,7.5,20)
anomal_eff = [CalculateEffPur(x)[0] for x in thresholds]
above_pur = [CalculateEffPur(x)[1] for x in thresholds]
normal_eff = [CalculateEffPur(x)[2] for x in thresholds]
below_pur = [CalculateEffPur(x)[3] for x in thresholds]
anomal_sig = [CalculateEffPur(x)[4] for x in thresholds]



KerasPur = np.array([0.42857142857142855, 0.5096525096525096, 0.5919282511210763, 0.6839378238341969, 0.75, 0.8300653594771242, 0.8857142857142857, 0.9172932330827067, 0.9375, 0.9583333333333334, 0.9557522123893806, 0.9532710280373832, 0.9702970297029703, 0.9787234042553191, 0.9777777777777777, 0.9759036144578314, 0.9743589743589743, 0.9722222222222222, 0.9836065573770492, 0.9830508474576272])
KerasEff = np.array([0.9777777777777777, 0.9777777777777777, 0.9777777777777777, 0.9777777777777777, 0.9555555555555556, 0.9407407407407408, 0.9185185185185185, 0.9037037037037037, 0.8888888888888888, 0.8518518518518519, 0.8, 0.7555555555555555, 0.725925925925926, 0.6814814814814815, 0.6518518518518519, 0.6, 0.562962962962963, 0.5185185185185185, 0.4444444444444444, 0.42962962962962964] )

fig, ax = plt.subplots()
plt.scatter(anomal_effAE, above_purAE, label='AutoEncoder', color=azureblue)
plt.scatter(anomal_eff, above_pur, label='<Cluster Size>', color=darkyellow)
plt.xlabel('Efficiency')
plt.ylabel('Purity')
plt.legend(loc='lower left', frameon = False)
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(0.15,0.8, "Run 3 enriched MC", fontsize=14)
plt.text(0.15,0.74, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/new/AutoEncoder/AnomalEffPur.png")


fig, ax = plt.subplots()
plt.scatter(anomal_effAE, above_purAE, label='PyTorch', color=azureblue)
plt.scatter(KerasEff, KerasPur, label='Keras', color=darkyellow)
plt.xlabel('Efficiency')
plt.ylabel('Purity')
plt.legend(loc='lower left', frameon = False)
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(0.15,0.8, "Run 3 enriched MC", fontsize=14)
plt.text(0.15,0.74, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/new/AutoEncoder/KerasPyTorch.png")


fig, ax = plt.subplots()
plt.scatter(thresholds, anomal_sig, label='Cut on Avg Cluster Size')
plt.xlabel('Cut')
plt.ylabel('Anomaly Significance')
plt.legend(loc='lower left')
plt.savefig("Plots/new/AutoEncoder/AnomalSig.png")

fig, ax = plt.subplots()
plt.scatter(thresholdsAE, anomal_sigAE, label='AutoEncoder')
plt.xlabel('Cut')
plt.ylabel('Anomaly Significance')
plt.legend(loc='lower left')
plt.savefig("Plots/new/AutoEncoder/AnomalSigAE.png")


fig, ax = plt.subplots()
plt.scatter(normal_effAE, below_purAE, label='AutoEncoder')
plt.scatter(normal_eff, below_pur, label='Cut on Avg Cl Size')
plt.xlabel('Normal Efficiency')
plt.ylabel('Normal Purity')
plt.legend(loc='lower left')
plt.savefig("Plots/new/AutoEncoder/NormalEffPur.png")

fig, ax = plt.subplots()
plt.hist(loss_dist, bins=50, label="MIP", color=azureblue, alpha =0.65)
plt.hist(loss_dist_anomal, bins=50, label="Nuclei", color='red', alpha =0.65)
plt.yscale('log')
plt.xlabel('AE score')
plt.legend(loc='upper right', frameon=False)
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(0.15,700, "Run 3 enriched MC", fontsize=14)
plt.text(0.15,400, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/new/AutoEncoder/LossPyTorch.png")


anomal_sigAE = [float(i)/max(anomal_sigAE) for i in anomal_sigAE]
fig, ax = plt.subplots()
plt.scatter(thresholdsAE, anomal_effAE, label='Efficiency', color='red')
plt.scatter(thresholdsAE, above_purAE, label='Purity', color=darkyellow)
ax.set_xlabel('AE Score')
ax.xaxis.set_label_coords(0.85,-0.1)
plt.legend(loc='lower center', frameon=False)
plt.ylim([0.4,1.1])
plt.text(0.1,1.05, "Run 3 enriched MC", fontsize=14)
plt.text(0.1,1.01, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/new/AutoEncoder/AnomalPerformanceAE.png")

anomal_sig = [float(i)/max(anomal_sig) for i in anomal_sig]
fig, ax = plt.subplots()
plt.scatter(thresholds, anomal_eff, label='Efficiency', color='red')
plt.scatter(thresholds, above_pur, label='Purity', color=darkyellow)
ax.set_xlabel('<Cluster Size>')
ax.xaxis.set_label_coords(0.85,-0.1)
plt.legend(loc='lower right', frameon=False)
plt.ylim([0.4,1.1])
plt.text(4.5,1.05, "Run 3 enriched MC", fontsize=14)
plt.text(4.5,1.01, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/new/AutoEncoder/AnomalPerformance.png")
