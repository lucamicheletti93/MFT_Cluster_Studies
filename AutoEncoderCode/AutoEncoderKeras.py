import ROOT
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf

from tensorflow.keras import layers, models
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model
from sklearn.preprocessing import MinMaxScaler, StandardScaler

root_file = ROOT.TFile("RootFiles/new/Input.root")
tree = root_file.Get("input")

branches = [ "PdgID", "Eta", "AvgClSize", "ClSizeLayer0", "ClSizeLayer1", "ClSizeLayer2", "ClSizeLayer3", "ClSizeLayer4", "ClSizeLayer5", "ClSizeLayer6", "ClSizeLayer7", "ClSizeLayer8", "ClSizeLayer9", "PattIDLayer0", "PattIDLayer1", "PattIDLayer2", "PattIDLayer3", "PattIDLayer4", "PattIDLayer5", "PattIDLayer6", "PattIDLayer7", "PattIDLayer8", "PattIDLayer9"]  

# Preprocessing
data = {}
for branch in branches:
    data[branch] = []
    for entry in tree:
        data[branch] = np.array(getattr(tree, branch))

df = pd.DataFrame.from_dict(data)

df = df.drop(df[(df["PdgID"] != 111) & (df["PdgID"]!=211) & (df["PdgID"]!=130) & (df["PdgID"]!=310) & (df["PdgID"]!=311) & (df["PdgID"]!=321) & (df["PdgID"]!=2212) & (df["PdgID"]!=13) & (df["PdgID"]!=1000020030) ].index)

def sortPdgID(PdgID):
    if (PdgID==111) | (PdgID==211) | (PdgID==130) | (PdgID==310) | (PdgID==311) | (PdgID==321) | (PdgID==2212) | (PdgID==13):
        PdgID=1
    elif PdgID==1000020030:
        PdgID=2
    return PdgID

df['PdgID'] = df['PdgID'].apply(sortPdgID)



x_train, x_test, y_train, y_test = train_test_split(df.values, df.values[:,0:1], test_size=0.2, random_state=111)

scaler = MinMaxScaler()
data_scaled = scaler.fit(x_train)
train_data_scaled = data_scaled.transform(x_train)
test_data_scaled = data_scaled.transform(x_test)

normal_train_data = pd.DataFrame(train_data_scaled).add_prefix('c').query('c0 == 0').values[:,1:]
anomaly_train_data = pd.DataFrame(train_data_scaled).add_prefix('c').query('c0 > 0').values[:, 1:]
normal_test_data = pd.DataFrame(test_data_scaled).add_prefix('c').query('c0 == 0').values[:,1:]
anomaly_test_data = pd.DataFrame(test_data_scaled).add_prefix('c').query('c0 > 0').values[:, 1:]



# Building the model

input_dim = len(normal_train_data[0])  

class AutoEncoder(Model):
    def __init__(self):
        super(AutoEncoder, self).__init__()
        self.encoder = tf.keras.Sequential([
        tf.keras.layers.Dense(4, activation="relu"),
        tf.keras.layers.Dense(9, activation="relu"),
        tf.keras.layers.Dense(2, activation="relu"),
        tf.keras.layers.Dense(8, activation="relu")
        ])
        self.decoder = tf.keras.Sequential([
        tf.keras.layers.Dense(2, activation="relu"),
        tf.keras.layers.Dense(9, activation="relu"),
        tf.keras.layers.Dense(4, activation="relu"),
        tf.keras.layers.Dense(input_dim, activation="sigmoid")
        ])
    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

model = AutoEncoder()
#early_stopping = tf.keras.callbacks.EarlyStopping(monitor="val_loss", patience=2, mode="min")
model.compile(optimizer='adam', loss="mae")
history = model.fit(normal_train_data, normal_train_data, epochs=50, batch_size=100,
        validation_split=0.2,
        shuffle=True,
        #   callbacks=[early_stopping]
        )

# Evaluation

reconstructed_data = model.predict(normal_test_data)
reconstruction_errors = np.mean(np.square(normal_test_data - reconstructed_data), axis=1)
test_loss = tf.keras.losses.mae(reconstructed_data, normal_test_data)



reconstructed_anomaly_data = model.predict(anomaly_test_data)
reconstruction_anomaly_errors = np.mean(np.square(anomaly_test_data - reconstructed_anomaly_data), axis=1)
test_loss_anomal = tf.keras.losses.mae(reconstructed_anomaly_data, anomaly_test_data)


test_loss=np.array(test_loss)
test_loss_anomal=np.array(test_loss_anomal)

def CalculateEffPurAE(threshold, loss_dist, loss_dist_anomal):

    normal_eff = (loss_dist < threshold).sum()/len(loss_dist)
    anomal_eff = (loss_dist_anomal > threshold).sum()/len(loss_dist_anomal)
    below_pur = (loss_dist < threshold).sum()/((loss_dist < threshold).sum()+(loss_dist_anomal < threshold).sum())
    above_pur = (loss_dist_anomal > threshold).sum()/((loss_dist > threshold).sum()+(loss_dist_anomal > threshold).sum())

    anomal_sig = (loss_dist_anomal > threshold).sum()/np.sqrt((loss_dist > threshold).sum()+(loss_dist_anomal > threshold).sum())

    return (anomal_eff, above_pur, normal_eff, below_pur, anomal_sig)

thresholdsAE = np.linspace(0.1,0.2,20)
anomal_effAE = [CalculateEffPurAE(x, test_loss, test_loss_anomal)[0] for x in thresholdsAE] 
above_purAE = [CalculateEffPurAE(x, test_loss, test_loss_anomal)[1] for x in thresholdsAE]
normal_effAE = [CalculateEffPurAE(x, test_loss, test_loss_anomal)[2] for x in thresholdsAE]
below_purAE = [CalculateEffPurAE(x, test_loss, test_loss_anomal)[3] for x in thresholdsAE]
anomal_sigAE = [CalculateEffPurAE(x, test_loss, test_loss_anomal)[4] for x in thresholdsAE]

for i in zip(above_purAE, anomal_effAE, anomal_sigAE, thresholdsAE):
    print(i)

print(above_purAE)
print(anomal_effAE)
# Plotting

plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.savefig("Plots/loss.png")
plt.clf()


fig, ax = plt.subplots()
plt.hist(test_loss, bins=50, label='normal')
plt.hist(test_loss_anomal, bins=50, label='anomaly')
plt.legend(loc='upper right')
plt.xlabel("Loss")
plt.title("Keras")
plt.savefig("Plots/LossComparison.png")

