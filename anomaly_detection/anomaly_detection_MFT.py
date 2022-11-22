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
from sklearn.metrics import confusion_matrix, precision_recall_curve, roc_curve, auc
from sklearn.preprocessing import StandardScaler
from hipe4ml import plot_utils
import matplotlib.pyplot as plt
import sys
from ROOT import TCanvas, TFile, TTree, TTreeReader, TTreeReaderArray, TH1F, TTreeReaderValue, TH2F, TLatex, TLegend, TF1
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack

DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Configurations
fInName = "../output/MFTAssessment_test.root"
treeName = "treeTrackClusterSize"
cols_to_keep = ['mIsMIP', 'mClsizeMean', 'mPattIdMean'] # mClsizeRecLayer taken by default, mClsizeMean is the mean of the cluster sizes
test_frac = 0.2
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
        cols_to_keep.append(f'mPattIdRecLayer{layer}')

        # normalize the cluster size (helps the AE)
        #df_origin[f'mClsizeRecLayer{layer}'] = df_origin[f'mClsizeRecLayer{layer}']/df_origin[f'mClsizeRecLayer{layer}'].max()
        if dropNan:
            df_origin[f'mClsizeRecLayer{layer}'] = np.where(df_origin[f'mClsizeRecLayer{layer}'] < 0,
                                                            float("nan"),
                                                            df_origin[f'mClsizeRecLayer{layer}'])

    if 'mClsizeMean' in cols_to_keep:
        df_origin['mClsizeMean'] = [(df_origin[f'mClsizeRecLayer0'][i] + df_origin[f'mClsizeRecLayer1'][i] + df_origin[f'mClsizeRecLayer2'][i] + df_origin[f'mClsizeRecLayer3'][i] + df_origin[f'mClsizeRecLayer4'][i] + df_origin[f'mClsizeRecLayer5'][i] + df_origin[f'mClsizeRecLayer6'][i] + df_origin[f'mClsizeRecLayer7'][i] + df_origin[f'mClsizeRecLayer8'][i] + df_origin[f'mClsizeRecLayer9'][i])/ 10 for i in range(len(df_origin[f'mClsizeRecLayer0']))]
        df_origin['mPattIdMean'] = [(df_origin[f'mPattIdRecLayer0'][i] + df_origin[f'mPattIdRecLayer1'][i] + df_origin[f'mPattIdRecLayer2'][i] + df_origin[f'mPattIdRecLayer3'][i] + df_origin[f'mPattIdRecLayer4'][i] + df_origin[f'mPattIdRecLayer5'][i] + df_origin[f'mPattIdRecLayer6'][i] + df_origin[f'mPattIdRecLayer7'][i] + df_origin[f'mPattIdRecLayer8'][i] + df_origin[f'mPattIdRecLayer9'][i])/ 10 for i in range(len(df_origin[f'mPattIdRecLayer0']))]
    else:
        print('WARNING: mClsizeMean not in cols_to_keep, the mean cluster size will not be computed')
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
df_nucl = df.query(anomaly_def) # select only the nucl
df_mip = df[~df.index.isin(df_nucl.index)] # select only the MIPs
plot_utils.plot_distr([df_nucl, df_mip], cols_to_keep, 100, ['Nuclei', 'MIP'], figsize=(24, 14),
                      alpha=0.4, log=True, grid=False, density=True)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig(f'{OutPutDir}/TrainVarsDistr.pdf')
plt.close('all')

df_mip_train = df_mip.sample(frac=1-test_frac, random_state=seed_split) # select the MIPs for training
df_mip_test = df_mip.drop(df_mip_train.index) # select the MIPs for testing

df_nucl_test = df_nucl # select the nucl for testing

df_test = pd.concat([df_nucl_test, df_mip_test]) # concatenate the MIPs and nucl for testing
df_test_labels = df_test['mIsMIP']  # define the labels for testing
df_test.drop(columns=['mIsMIP'], inplace=True) # drop the labels from the testing set

df_mip_train = df_mip_train.drop(columns=['mIsMIP']) # drop the label
df_mip_test = df_mip_test.drop(columns=['mIsMIP']) # drop the label
df_nucl_test = df_nucl_test.drop(columns=['mIsMIP']) # drop the label

df_mip_train = df_mip_train.to_numpy() # convert the training set to numpy
df_mip_test = df_mip_test.to_numpy() # convert the testing set to numpy
df_nucl_test = df_nucl_test.to_numpy() # convert the testing set to numpy
df_test = df_test.to_numpy() # convert the testing set to numpy

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
                       tf.keras.layers.Dense(input_dim, activation='relu')
                  ])
    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

model = AutoEncoder(df_mip_train.shape[1])
model.compile(optimizer='adam', loss="mae")
history = model.fit(df_mip_train, df_mip_train, epochs=20, batch_size=100,
                    validation_data=(df_mip_test, df_mip_test),
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

reconstruction = model.predict(df_mip_test) # apply the model to the MIPs --> low loss expected
train_loss = tf.keras.losses.mae(reconstruction, df_mip_test)

threshold = np.mean(train_loss) + 2*np.std(train_loss)
reconstruction_a = model.predict(df_nucl_test) # apply the model to the nucl --> high loss expected
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

# Plot the loss
fig = plt.figure()
plt.hist(reconstruction, bins=100, density=True, alpha=0.5, label='MIP')
plt.hist(reconstruction_a, bins=100, density=True, alpha=0.5, label='Nuclei')
plt.legend(loc='upper right')
plt.xlabel('Mean Cl size')
plt.ylabel('Pred')
fig.savefig(f'{OutPutDir}/Pred.pdf')
# Plot the Precision-Recall curve
fig = plt.figure()
valid_x_predictions = model.predict(df_test)
mse = np.mean(np.power(df_test - valid_x_predictions, 2), axis=1)
error_df = pd.DataFrame({'Reconstruction_error': mse,
                         'True_class': df_test_labels})
precision_rt, recall_rt, threshold_rt = precision_recall_curve(error_df.True_class,
                                                               error_df.Reconstruction_error)
plt.plot(threshold_rt, precision_rt[1:], label="Precision",linewidth=5)
plt.plot(threshold_rt, recall_rt[1:], label="Recall",linewidth=5)
plt.title('Precision and recall for different threshold values')
plt.xlabel('Threshold')
plt.ylabel('Precision/Recall')
plt.legend()
fig.savefig(f'{OutPutDir}/PrecisionRecall.pdf')
# Plot the ROC curve
fig = plt.figure()
false_pos_rate, true_pos_rate, thresholds = roc_curve(error_df.True_class,
                                                      error_df.Reconstruction_error)
roc_auc = auc(false_pos_rate, true_pos_rate)
print('AUC: ', roc_auc)
plt.plot( false_pos_rate, true_pos_rate, linewidth=5, label='AUC = %0.3f'% roc_auc)
plt.plot([0,1],[0,1], linewidth=5)
plt.xlim([-0.01, 1])
plt.ylim([0, 1.01])
plt.legend(loc='lower right')
plt.title('Receiver operating characteristic curve (ROC)')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
fig.savefig(f'{OutPutDir}/ROC.pdf')
# Plot the loss hitory
fig = plt.figure()
plt.plot(history.history['loss'], label='train')
plt.plot(history.history['val_loss'], label='test')
plt.legend()
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.title('Loss vs Epoch')
plt.savefig(f'{OutPutDir}/LossVsEpoch.pdf')

fig = plt.figure()
train_loss = pd.DataFrame(train_loss)
train_loss_a = pd.DataFrame(train_loss_a)
plt.hist(train_loss, label='normal', log=True, alpha=0.5, bins=range(0, 10, 100))
plt.hist(train_loss_a, label='anomaly', log=True, alpha=0.5, bins=range(0, 10, 100))
plt.axvline(threshold, color='r', linewidth=3, linestyle='dashed', label='{:0.3f}'.format(threshold))
plt.legend(loc='upper right')
plt.title("Normal and Anomaly Loss")
fig.savefig("ResultsMFT.pdf")

hLoss = TH1F("hLoss", "hLoss", 100, 0, 1)
hLoss_a = TH1F("hLoss_a", "hLoss_a", 100, 0, 1)
for i in range(len(train_loss)):
    hLoss.Fill(train_loss.iloc[i][0])
for i in range(len(train_loss_a)):
    hLoss_a.Fill(train_loss_a.iloc[i][0])
canv  = TCanvas("canv", "canv", 800, 600)
leg = TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(hLoss, "MIP", "l")
leg.AddEntry(hLoss_a, "Nuclei", "l")
hLoss.SetLineColor(kRed)
hLoss_a.SetLineColor(kBlue)
hLoss.Draw()
hLoss_a.Draw("same")
canv.SaveAs("ResultsMFT_root.pdf")

print('\r\033[1m' + 'Plotting complete.' + '\033[0m')

input("Press Enter to exit")
sys.exit()
