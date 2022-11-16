import pandas as pd
import numpy as np
import uproot
from alive_progress import alive_bar
import torch
import tensorflow as tf
from tensorflow.keras.models import Model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from hipe4ml import plot_utils
import matplotlib.pyplot as plt
import sys
from ROOT import TCanvas, TFile, TTree, TTreeReader, TTreeReaderArray, TH1F, TTreeReaderValue, TH2F, TLatex, TLegend, TF1
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack

DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Configurations
fInName = "../output/MFTAssessment_test.root"
treeName = "treeTrackClusterSize"
cols_to_keep = ['mIsMIP'] # mClsizeRecLayer taken by default
test_frac = 0.5
seed_split = 42
OutPutDir = "./"
anomaly_def = 'mIsMIP == False'
doRootToParquet = True # once you have the parquet file, set this to False

# Data preparation
print('\033[1m' + 'Preparing data...' + '\033[0m')

# Read the ROOT file and convert it to a pandas dataframe
if doRootToParquet:
    df_origin = uproot.open(fInName)[treeName].arrays(library='pd')
    for layer in range(0, 10):
        cols_to_keep.append(f'mClsizeRecLayer{layer}')

        # normalize the cluster size (helps the AE)
        df_origin[f'mClsizeRecLayer{layer}'] = df_origin[f'mClsizeRecLayer{layer}']/df_origin[f'mClsizeRecLayer{layer}'].max()
        df_origin[f'mClsizeRecLayer{layer}'] = np.where(df_origin[f'mClsizeRecLayer{layer}'] < 0,
                                                        float("nan"),
                                                        df_origin[f'mClsizeRecLayer{layer}'])

    df = df_origin[cols_to_keep] # drop columns not needed
    df = df.dropna() # drop rows with nan
    df['mClsizeMean'] = df.iloc[:, 1:11].mean(axis=1) # add the mean cluster size
    cols_to_keep = ['mIsMIP', 'mClsizeMean'] # mClsizeRecLayer and mPattIdRecLayer taken by default
    df = df[cols_to_keep]
    print(df)
    input('Press Enter to continue...')
    #df = df.sample(frac=0.5) # drop data to speed up the process (for testing)
else:
    df = pd.read_parquet(f'{OutPutDir}parquet/df_MC.parquet.gzip')


# Split the data into training and testing sets
df_nucl = df.query(anomaly_def)
df_mip = df[~df.index.isin(df_nucl.index)] # Select all the rows that are not in df_nucl
df_nucl_train = df_nucl.sample(frac=1-test_frac, random_state=seed_split)
df_mip_train = df_mip.sample(frac=1-test_frac, random_state=seed_split)
df_nucl_test = df_nucl.drop(df_nucl_train.index)
df_mip_test = df_mip.drop(df_mip_train.index)


# Plot training variables for MIP and nuclei
plot_utils.plot_distr([df_nucl, df_mip], cols_to_keep, 100, ['Nuclei', 'MIP'], figsize=(24, 14),
                      alpha=0.4, log=True, grid=False, density=True)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig(f'{OutPutDir}/TrainVarsDistr.pdf')
plt.close('all')

df_nucl_train = df_nucl_train.drop(['mIsMIP'], axis=1)
df_mip_train = df_mip_train.drop(['mIsMIP'], axis=1)
df_nucl_test = df_nucl_test.drop(['mIsMIP'], axis=1)
df_mip_test = df_mip_test.drop(['mIsMIP'], axis=1)
df_nucl_train = df_nucl_train.to_numpy()
df_mip_train = df_mip_train.to_numpy()
df_nucl_test = df_nucl_test.to_numpy()
df_mip_test = df_mip_test.to_numpy()

print('\r\033[1m' + 'Data preparation complete.' + '\033[0m')


# Model preparationfrom sklearn.metrics import accuracy_score
print('\033[1m' + 'Training the model...' + '\033[0m')
class AutoEncoder(Model):
    def __init__(self, input_dim):
        super(AutoEncoder, self).__init__()
        self.encoder = tf.keras.Sequential([
                      tf.keras.layers.Dense(input_dim*6, activation='relu'),
                      tf.keras.layers.Dense(input_dim*5, activation="relu"),
                      tf.keras.layers.Dense(input_dim*4, activation="relu"),
                      tf.keras.layers.Dense(input_dim*3, activation="relu"),
                      tf.keras.layers.Dense(input_dim*2, activation="relu"),
                      tf.keras.layers.Dense(input_dim, activation="relu")
                  ])
        self.decoder = tf.keras.Sequential([
                       tf.keras.layers.Dense(input_dim, activation="relu"),
                       tf.keras.layers.Dense(input_dim*2, activation="relu"),
                       tf.keras.layers.Dense(input_dim*3, activation="relu"),
                       tf.keras.layers.Dense(input_dim*4, activation="relu"),
                       tf.keras.layers.Dense(input_dim*5, activation="relu"),
                       tf.keras.layers.Dense(input_dim*6, activation='relu'),
                       tf.keras.layers.Dense(input_dim, activation='sigmoid')
                  ])
    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

model = AutoEncoder(df_nucl_train.shape[1])
model.compile(optimizer='adam', loss="mae")
history = model.fit(df_mip_train, df_mip_train, epochs=450, batch_size=50,
                    validation_data=(df_nucl_test, df_nucl_test),
                    verbose=True, shuffle=True)
print('\r\033[1m' + 'Model training complete.' + '\033[0m')


# Apply the model to the data
print('\033[1m' + 'Applying the model to the data...' + '\033[0m')
# mip
encoder_out = model.encoder(df_mip_test).numpy()
decoder_out = model.decoder(encoder_out).numpy()
# nucl
encoder_out_a = model.encoder(df_nucl_test).numpy()
decoder_out_a = model.decoder(encoder_out_a).numpy()

reconstruction = model.predict(df_mip_test)
train_loss = tf.keras.losses.mae(reconstruction, df_mip_test)

threshold = np.mean(train_loss) + 2*np.std(train_loss)
reconstruction_a = model.predict(df_nucl_test)
train_loss_a = tf.keras.losses.mae(reconstruction_a, df_nucl_test)

df_decoder_out = pd.DataFrame(decoder_out)
df_decoder_out_a = pd.DataFrame(decoder_out_a)
print('\r\033[1m' + 'Decoding complete:' + '\033[0m')
print(df_decoder_out.head(5))
print(df_decoder_out_a.head(5))
input('Press Enter to continue...')

print('\r\033[1m' + 'Model application complete.' + '\033[0m')


# Plot the results
print('\033[1m' + 'Plotting the results...' + '\033[0m')

fig = plt.figure()
plt.plot(history.history['loss'], label='train')
plt.plot(history.history['val_loss'], label='test')
plt.legend()
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('Loss vs Epoch')
plt.savefig(f'{OutPutDir}/LossVsEpoch.pdf')

fig = plt.figure()
plt1 = fig.add_subplot(221)
plt2 = fig.add_subplot(222)
plt3 = fig.add_subplot(223)
plt4 = fig.add_subplot(224)

plt1.plot(df_mip_test[0], 'b')
plt1.plot(decoder_out[0], 'r')
plt1.set_title("Model performance on Normal data")

plt2.plot(df_nucl_test[0], 'b')
plt2.plot(decoder_out_a[0], 'r')
plt2.set_title("Model performance on Anomaly data")

plt3.hist(train_loss_a, bins=50)
plt3.set_title("loss on anomaly test data")

plt4.hist(train_loss, bins=50, label='normal')
plt4.hist(train_loss_a, bins=50, label='anomaly')
plt4.axvline(threshold, color='r', linewidth=3, linestyle='dashed', label='{:0.3f}'.format(threshold))
plt4.legend(loc='upper right')
plt4.set_title("Normal and Anomaly Loss")
fig.savefig("Results.pdf")

fig = plt.figure()
train_loss = pd.DataFrame(train_loss)
train_loss_a = pd.DataFrame(train_loss_a)
plt.hist(train_loss, bins=20, label='normal')
plt.hist(train_loss_a, bins=20, label='anomaly')
#plt.axvline(threshold, color='r', linewidth=3, linestyle='dashed', label='{:0.3f}'.format(threshold))
plt.legend(loc='upper right')
plt.title("Normal and Anomaly Loss")
fig.savefig("ResultsMFT.pdf")

print('\r\033[1m' + 'Plotting complete.' + '\033[0m')

input("Press Enter to exit")
sys.exit()
