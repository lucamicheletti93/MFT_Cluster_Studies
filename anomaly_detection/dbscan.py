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
from sklearn import metrics
from sklearn.metrics import confusion_matrix, precision_recall_curve, roc_curve, auc
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from hipe4ml import plot_utils
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib
import sys
from ROOT import TCanvas, TFile, TTree, TTreeReader, TTreeReaderArray, TH1F, TTreeReaderValue, TH2F, TLatex, TLegend, TF1
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack

DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Configurations
fInName = "../output/MFTAssessmentMC.root"
treeName = "treeTrackClusterSize"
cols_to_keep = ['mIsMIP', 'mClsizeMean', 'mPattIdMean'] # mClsizeRecLayer taken by default, mClsizeMean is the mean of the cluster sizes
test_frac = 0.2
seed_split = 42
dropNan = True # if True, drop the rows with NaN values
doNorm = False # if True, normalize the data with MinMaxScaler
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
    if dropNan:
        print('Dropping NaN values...')
        df = df.dropna() # drop rows with NaN values
    if doNorm:
        df = pd.DataFrame(MinMaxScaler().fit_transform(df), columns=df.columns)
    print(df)
else:
    df = pd.read_parquet(f'{OutPutDir}parquet/df_MC.parquet.gzip')

labels_true = df['mIsMIP'].values
labels_true = np.where(labels_true == True, 0, -1)
print('labels_true', labels_true)
df = df.drop(columns=['mIsMIP']) # drop the column with the labels

# n_neighbors = 5 as kneighbors function returns distance of point to itself (i.e. first column will be zeros)
nbrs = NearestNeighbors(n_neighbors=42).fit(df)
# Find the k-neighbors of a point
neigh_dist, neigh_ind = nbrs.kneighbors(df)
# sort the neighbor distances (lengths to points) in ascending order
# axis = 0 represents sort along first axis i.e. sort along row
sort_neigh_dist = np.sort(neigh_dist, axis=0)

k_dist = sort_neigh_dist[:, 40] # 44 is the 45th element (index 44) as first element is 0
fig = plt.figure()
plt.plot(k_dist)
plt.axhline(y=10, linewidth=1, linestyle='dashed', color='k')
plt.ylabel("k-NN distance")
plt.xlabel("Sorted observations (4th NN)")
fig.savefig("dbs_map.pdf")

db = DBSCAN(eps=9, min_samples=20).fit(df)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
print(labels)

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0) # -1 is the noise
colors = [color for name, color in mcolors.TABLEAU_COLORS.items()]  # more colors: 'mcolors.XKCD_COLORS'
n_noise_ = list(labels).count(-1)

df_out = df
df_out['mIsMIP'] = labels_true
df_out['mIsMIP_pred'] = labels
print(f'df_out: {df_out}')

# compare predicted labels with true labels distribution
df_nucl_true = df_out[df_out['mIsMIP'] == -1]
df_mip_true = df_out[df_out['mIsMIP'] == 0]
df_nucl_pred = df_out[df_out['mIsMIP_pred'] == -1]
df_mip_pred = df_out[df_out['mIsMIP_pred'] == 0]

plot_utils.plot_distr([df_nucl_true, df_nucl_pred], cols_to_keep, 100, ['True', 'Pred'], figsize=(24, 14),
                      alpha=0.4, log=True, grid=False, density=True)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig(f'{OutPutDir}/dbscan_true_pred_nucl.pdf')
plt.close('all')

plot_utils.plot_distr([df_mip_true, df_mip_pred], cols_to_keep, 100, ['True', 'Pred'], figsize=(24, 14),
                      alpha=0.4, log=True, grid=False, density=True)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig(f'{OutPutDir}/dbscan_true_pred_mip.pdf')
plt.close('all')


fig = plt.figure()
plt.scatter(df['mClsizeRecLayer0'], df['mClsizeRecLayer9'], c=labels, cmap=matplotlib.colors.ListedColormap(colors), s=10)
fig.savefig("dbs_corr.pdf")

fig = plt.figure()
delta = np.abs(labels_true - labels)
plt.hist(delta, bins=5, range=(-2.5, 2.5))
fig.savefig("dbs_hist.pdf")

fig = plt.figure()
plt.hist(labels_true, bins=5, range=(-2.5, 2.5))
plt.hist(labels, bins=5, range=(-2.5, 2.5), alpha=0.5)
fig.savefig("dbs_labels.pdf")

print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of noise points: %d' % n_noise_)
print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
print("Adjusted Rand Index: %0.3f"
      % metrics.adjusted_rand_score(labels_true, labels))
print("Adjusted Mutual Information: %0.3f"
      % metrics.adjusted_mutual_info_score(labels_true, labels))
print("Signal Efficiency: %0.3f"% metrics.accuracy_score(labels_true, labels)) # signal efficiency = TP/(TP+FN)

input('Press Enter to exit')
sys.exit()