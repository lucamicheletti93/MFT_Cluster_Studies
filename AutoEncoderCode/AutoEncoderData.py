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


def getData(inputfile):
    root_file = ROOT.TFile("RootFiles/new/"+inputfile)
    tree = root_file.Get("input")

    branches = [ "AvgClSize", "Eta", "ClSizeLayer0", "ClSizeLayer1", "ClSizeLayer2", "ClSizeLayer3", "ClSizeLayer4", "ClSizeLayer5", "ClSizeLayer6", "ClSizeLayer7",  "PattIDLayer0", "PattIDLayer1", "PattIDLayer2", "PattIDLayer3", "PattIDLayer4", "PattIDLayer5", "PattIDLayer6", "PattIDLayer7"]  

    # Preprocessing
    data_dict = {}
    for branch in branches:
        data_dict[branch] = []
        for entry in tree:
            data_dict[branch] = np.array(getattr(tree, branch))


    df = pd.DataFrame.from_dict(data_dict)
    return df

df = getData("DataInput.root")


original_array = df.values

df_no_selection = getData("DataInputNoSelection.root")

no_selection = df_no_selection.values

df_MFT = getData("DataInputMFT.root")
print(len(df_MFT[df_MFT['AvgClSize']>7]))

MFT = df_MFT.values 


root_file = ROOT.TFile("RootFiles/new/Input.root")
tree = root_file.Get("input")

branches = [ "PdgID",  "AvgClSize", "Eta", "ClSizeLayer0", "ClSizeLayer1", "ClSizeLayer2", "ClSizeLayer3", "ClSizeLayer4", "ClSizeLayer5", "ClSizeLayer6", "ClSizeLayer7",  "PattIDLayer0", "PattIDLayer1", "PattIDLayer2", "PattIDLayer3", "PattIDLayer4", "PattIDLayer5", "PattIDLayer6", "PattIDLayer7"]  

# Preprocessing
data_dict = {}
for branch in branches:
    data_dict[branch] = []
    for entry in tree:
        data_dict[branch] = np.array(getattr(tree, branch))


MCdf = pd.DataFrame.from_dict(data_dict)
normal_df = MCdf[(MCdf["PdgID"] == 111) | (MCdf["PdgID"]==211) | (MCdf["PdgID"]==130) | (MCdf["PdgID"]==310) | (MCdf["PdgID"]==311) | (MCdf["PdgID"]==321) | (MCdf["PdgID"]==2212) | (MCdf["PdgID"]==13)]
anomal_df = MCdf[(MCdf["PdgID"] == 1000020030) ]

normal_test_data = normal_df.values[:,1:]
anomaly_test_data = anomal_df.values[:,1:]




# Get the number of elements in the array
num_elements = len(original_array)

# Shuffle the indices of the array
shuffled_indices = np.random.permutation(num_elements)

# Calculate the split index for 80/20 ratio
split_index = int(0.8 * num_elements)

# Split the shuffled indices into train and test indices
train_indices = shuffled_indices[:split_index]
test_indices = shuffled_indices[split_index:]

# Extract train and test arrays using the indices
train_array = original_array[train_indices]
test_array = original_array[test_indices]

data_tensor = torch.FloatTensor(train_array).unsqueeze(1) 
#data_test_tensor =  torch.FloatTensor(normal).unsqueeze(1) 
dim = len(train_array[0])

validation_percentage = 0.2

num_validation = int(validation_percentage * len(data_tensor))

num_elements = len(data_tensor)

# Shuffle the indices of the array
shuffled_indices = np.random.permutation(num_elements)

# Calculate the split index for 80/20 ratio
split_index = int(0.8 * num_elements)

# Split the shuffled indices into train and test indices
train_indices = shuffled_indices[:split_index]
val_indices = shuffled_indices[split_index:]


# Split the data into training and validation sets
train_data = data_tensor[train_indices]
val_data = data_tensor[val_indices]

layers = 5
nodes_per_layers = [5,3,3,7,6,16]
lr = 0.03


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
val_loader = data.DataLoader(val_data, batch_size=batch_size)

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
    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}")

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

device = 'cuda' if torch.cuda.is_available() else 'cpu'

# Evaluate the model on the entire dataset to get the reconstruction errors
model.eval()

def testAE(test_array):
    loss_dist=[]
    for i in range(len(test_array)):
        data_test = torch.from_numpy(test_array[i]).float()
        sample_test = model(data_test.to(device))
        loss = criterion(data_test.to(device), sample_test)
        loss_dist.append(loss.item())
    loss_dist = np.array(loss_dist)
    return loss_dist


loss_dist_anomal = testAE(anomaly_test_data)
loss_dist = testAE(normal_test_data)
loss_data_dist = testAE(test_array)
loss_no_selection_dist= testAE(no_selection)
loss_MFT = testAE(MFT)

azureblue = (85./255, 170./255, 216./255)
darkyellow=(0.800000011920929, 0.800000011920929, 0.20000000298023224)


plt.figure(figsize=(8, 6))
plt.plot(train_losses, label='Training Loss')
plt.plot(val_losses, label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.title('Training and Validation Loss over Epochs')
plt.savefig('Plots/data/AutoEncoder/Loss.png')
plt.close()


#normal_train_data = normal_train_data[:100]
#data_shap = torch.tensor(normal_train_data, dtype=torch.float32)
#explainer = shap.DeepExplainer(model, data_shap)

# Calculate SHAP values for the test set
#shap_values = explainer.shap_values(data_shap)

#feature_names = [branches[i] for i in range(2,len(normal_train_data[0])+1)]
#shap.summary_plot(shap_values[0], data_shap, plot_type="bar", show=False, feature_names=feature_names)
#print('shap_values', shap_values)

#plt.savefig('Plots/Chi2/AutoEncoder/shap_summary_plot.png')
#plt.close()




def CalculateEffPurAE(threshold, loss_dist, loss_dist_anomal):

    normal_eff = (loss_dist < threshold).sum()/len(loss_dist)
    anomal_eff = (loss_dist_anomal > threshold).sum()/len(loss_dist_anomal)
    below_pur = (loss_dist < threshold).sum()/((loss_dist < threshold).sum()+(loss_dist_anomal < threshold).sum())
    above_pur = (loss_dist_anomal > threshold).sum()/((loss_dist > threshold).sum()+(loss_dist_anomal > threshold).sum())

    anomal_sig = (loss_dist_anomal > threshold).sum()/np.sqrt((loss_dist > threshold).sum()+(loss_dist_anomal > threshold).sum())

    return (anomal_eff, above_pur, normal_eff, below_pur, anomal_sig)


def CalculateEffPur(AvgClThreshold):


    below_df = MCdf[(MCdf["AvgClSize"]<AvgClThreshold)]
    above_df = MCdf[(MCdf["AvgClSize"]>AvgClThreshold)]

    normal_eff = len(normal_df[normal_df["AvgClSize"] < AvgClThreshold])/len(normal_df)
    anomal_eff = len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold])/len(anomal_df)
    below_pur = len(normal_df[normal_df["AvgClSize"] < AvgClThreshold])/(len(normal_df[normal_df["AvgClSize"] < AvgClThreshold])+len(anomal_df[anomal_df["AvgClSize"] < AvgClThreshold]))
    above_pur = len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold])/(len(normal_df[normal_df["AvgClSize"] > AvgClThreshold])+len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold]))
    anomal_sig = len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold])/np.sqrt(len(normal_df[normal_df["AvgClSize"] > AvgClThreshold])+len(anomal_df[anomal_df["AvgClSize"] > AvgClThreshold]))
    return (anomal_eff, above_pur, normal_eff, below_pur, anomal_sig)



thresholdsAE = np.linspace(5,12,20)
anomal_effAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[0] for x in thresholdsAE] 
above_purAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[1] for x in thresholdsAE]
normal_effAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[2] for x in thresholdsAE]
below_purAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[3] for x in thresholdsAE]
anomal_sigAE = [CalculateEffPurAE(x, loss_dist, loss_dist_anomal)[4] for x in thresholdsAE]
for i in zip(above_purAE, anomal_effAE, anomal_sigAE, thresholdsAE):
    print(i)


thresholds = np.linspace(2.5,7,20)
anomal_eff = [CalculateEffPur(x)[0] for x in thresholds]
above_pur = [CalculateEffPur(x)[1] for x in thresholds]
normal_eff = [CalculateEffPur(x)[2] for x in thresholds]
below_pur = [CalculateEffPur(x)[3] for x in thresholds]
anomal_sig = [CalculateEffPur(x)[4] for x in thresholds]




fig, ax = plt.subplots()
plt.scatter(anomal_effAE, above_purAE, label='AutoEncoder', color=azureblue)
plt.scatter(anomal_eff, above_pur, label='<Cluster Size>', color=darkyellow)
plt.xlabel('Efficiency')
plt.ylabel('Purity')
plt.legend(loc='lower left', frameon = False)
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(0.45,0.8, "Run 3 enriched MC/Data", fontsize=14)
plt.text(0.45,0.74, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/data/AutoEncoder/AnomalEffPur.png")

fig, ax = plt.subplots()
plt.scatter(thresholds, anomal_sig, label='Cut on Avg Cluster Size')
plt.xlabel('Cut')
plt.ylabel('Anomaly Significance')
plt.legend(loc='lower left')
plt.savefig("Plots/data/AutoEncoder/AnomalSig.png")

fig, ax = plt.subplots()
plt.scatter(thresholdsAE, anomal_sigAE, label='AutoEncoder')
plt.xlabel('Cut')
plt.ylabel('Anomaly Significance')
plt.legend(loc='lower left')
plt.savefig("Plots/data/AutoEncoder/AnomalSigAE.png")



fig, ax = plt.subplots()
plt.hist(loss_dist, bins=50, label="enriched MC MIP", color=azureblue, alpha =0.65)
plt.hist(loss_dist_anomal, bins=50, label="enriched MC Nuclei", color='red', alpha =0.65)
plt.hist(loss_data_dist, bins=50, label="Data", color=darkyellow, alpha=0.65)
plt.yscale('log')
ax.set_ylim([1,11000])
plt.xlabel('AE score')
plt.legend(loc='upper right', frameon=False)
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(1.5,7000, "Run 3 enriched MC / Data", fontsize=14)
plt.text(1.5,4000, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/data/AutoEncoder/LossPyTorch.png")


fig, ax = plt.subplots()
plt.hist(loss_no_selection_dist, bins=50, label="Global")
ax.set_yscale('log')
plt.xlabel('Loss')
plt.legend(loc='upper right')
plt.savefig("Plots/data/AutoEncoder/LossDataGlobal.png")

print("total tracks", len(loss_MFT))
index = np.argmax(loss_MFT)
print(MFT[index])
fig, ax = plt.subplots()
plt.hist(loss_MFT, bins=100, label="MFT", color=azureblue)
ax.set_yscale('log')
#ax.set_xscale('log')
plt.xlabel('AE score')
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(1000000,90000, "Run 3 ", fontsize=14)
plt.text(1000000,30000, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.legend(loc='upper right', frameon=False)
plt.savefig("Plots/data/AutoEncoder/LossDataMFT.png")

loss_MFT = loss_MFT[(loss_MFT < 50)]
fig, ax = plt.subplots()
plt.hist(loss_MFT, bins=100, label="MFT", color=azureblue)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_xlim([0,50])
plt.xlabel('AE score')
ax.xaxis.set_label_coords(0.85,-0.1)
plt.text(5,70000, "Run 3 ", fontsize=14)
plt.text(5,40000, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.legend(loc='upper right', frameon=False)
plt.savefig("Plots/data/AutoEncoder/LossDataMFTzoomed.png")


anomal_sigAE = [float(i)/max(anomal_sigAE) for i in anomal_sigAE]
fig, ax = plt.subplots()
plt.scatter(thresholdsAE, anomal_effAE, label='Efficiency', color='red')
plt.scatter(thresholdsAE, above_purAE, label='Purity', color=darkyellow)
ax.set_xlabel('AE Score')
ax.xaxis.set_label_coords(0.85,-0.1)
plt.legend(loc='lower center', frameon=False)
plt.ylim([0.4,1.1])
plt.text(5,1.05, "Run 3 enriched MC/ Data", fontsize=14)
plt.text(5,1.01, r"pp, $\sqrt{s}$=13 TeV", fontsize=14)
plt.savefig("Plots/data/AutoEncoder/AnomalPerformanceAE.png")

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
plt.savefig("Plots/data/AutoEncoder/AnomalPerformance.png")
