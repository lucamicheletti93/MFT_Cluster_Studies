import pandas as pd
import numpy as np
import uproot
from alive_progress import alive_bar
import torch
import seaborn as sn
import tensorflow as tf
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from hipe4ml import plot_utils
import matplotlib.pyplot as plt
import sys
from ROOT import TCanvas, TFile, TTree, TTreeReader, TTreeReaderArray, TH1F, TTreeReaderValue, TH2F, TLatex, TLegend, TF1
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack, kFullCircle, kOpenCircle, kFullSquare, kOpenSquare
sys.path.append('..')
from utils.StyleFormatter import SetObjectStyle
from utils.MLmodels import AutoEncoder, MLplots, plot_confusion_matrix
from utils.plot_library import plot_hist_from_df, plot_hist2d_from_df


#____________________________________________
# Configurations
fInName = "../output/MFTAssessmentMC.root"
treeName = "treeTrackClusterSize"
cols_to_keep = ['mClsizeRecLayer0', 'mClsizeRecLayer1', 'mClsizeRecLayer2',
                'mClsizeRecLayer3', 'mClsizeRecLayer4', 'mClsizeRecLayer5',
                'mClsizeRecLayer6', 'mClsizeRecLayer7', 'mClsizeRecLayer8',
                'mClsizeRecLayer9', 'mPattIdRecLayer0', 'mPattIdRecLayer1',
                'mPattIdRecLayer2', 'mPattIdRecLayer3', 'mPattIdRecLayer4',
                'mPattIdRecLayer5', 'mPattIdRecLayer6', 'mPattIdRecLayer7',
                'mPattIdRecLayer8', 'mPattIdRecLayer9',
                'mClsizeMean', 'mPattIdMean', 'mIsMIP'] # columns to keep from the input tree
seed_split = 42
dropNan = True   # if True, drop the rows with NaN values
doNorm = True # if True, normalize the data with MinMaxScaler; automatically set to True if model == 'AE'
OutPutDir = "./"
anomaly_def = 'mIsMIP == False'
test_frac = 0.5
doRootToParquet = True # dump Tree to .parquet: once the .parquet file exists, set this to False to speed up data preparation
mlmodel = 'AE' # 'DBSCAN' or 'AE'

#____________________________________________
# Data preparation
if doRootToParquet:
    print('\033[1m' + 'Preparing data...' + '\033[0m')
    df_origin = uproot.open(fInName)[treeName].arrays(library='pd')

    # Quality cuts
    print('Considering only tracks with at least 4 hits...')
    n_tracks = len(df_origin)
    for irow, (cl0, cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8, cl9) in enumerate(zip(df_origin['mClsizeRecLayer0'], df_origin['mClsizeRecLayer1'], df_origin['mClsizeRecLayer2'], df_origin['mClsizeRecLayer3'], df_origin['mClsizeRecLayer4'], df_origin['mClsizeRecLayer5'], df_origin['mClsizeRecLayer6'], df_origin['mClsizeRecLayer7'], df_origin['mClsizeRecLayer8'], df_origin['mClsizeRecLayer9'])):
        cls = [cl0, cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8, cl9]
        while -999 in cls:
            cls.remove(-999)
        if len(cls) < 5:
            df_origin.drop(irow)
    print(f'Total number of tracks: {n_tracks}')
    print(f'Number of tracks with at least 4 hits: {len(df_origin)}')
    print(f'% of tracks with at least 4 hits: {len(df_origin)/n_tracks*100:.2f}%')

    # Additional df info
    df_origin['mClsizeMean'] = df_origin['mMeanClsizePerTrackRec'] # dummy, to be fixed
    df_origin['mPattIdMean'] = [(df_origin[f'mPattIdRecLayer0'][i] + df_origin[f'mPattIdRecLayer1'][i] + df_origin[f'mPattIdRecLayer2'][i] + df_origin[f'mPattIdRecLayer3'][i] + df_origin[f'mPattIdRecLayer4'][i] + df_origin[f'mPattIdRecLayer5'][i] + df_origin[f'mPattIdRecLayer6'][i] + df_origin[f'mPattIdRecLayer7'][i] + df_origin[f'mPattIdRecLayer8'][i] + df_origin[f'mPattIdRecLayer9'][i])/ 10 for i in range(len(df_origin[f'mPattIdRecLayer0']))]
    if doNorm or mlmodel == 'AE': # normalize the data, especially for the AE
        df_origin = pd.DataFrame(MinMaxScaler().fit_transform(df_origin), columns=df_origin.columns)
    hit_mean_number = [] # mean number of hits per track
    for i, (cl0, cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8, cl9) in enumerate(zip(df_origin['mClsizeRecLayer0'], df_origin['mClsizeRecLayer1'], df_origin['mClsizeRecLayer2'], df_origin['mClsizeRecLayer3'], df_origin['mClsizeRecLayer4'], df_origin['mClsizeRecLayer5'], df_origin['mClsizeRecLayer6'], df_origin['mClsizeRecLayer7'], df_origin['mClsizeRecLayer8'], df_origin['mClsizeRecLayer9'])):
        cl_list = [cl0, cl1, cl2, cl3, cl4, cl5, cl6, cl7, cl8, cl9]
        hit_mean_number.append(len([cl for cl in cl_list if cl > 0])/10)
    print('Mean number of hits per track: ', np.mean(hit_mean_number))
    df_origin['mMeanHitNumberPerTrackRec'] = hit_mean_number

    # dropping columns
    for col in df_origin.columns:
        if col not in cols_to_keep:
            df_origin.drop(col, axis=1, inplace=True)

    # dropping NaN values
    if dropNan:
        df = df_origin
        print('Dropping NaN values...')
        for layer in range(0, 10): # loop over the layers
            if dropNan:
                df[f'mClsizeRecLayer{layer}'] = np.where(df[f'mClsizeRecLayer{layer}'] < 0,
                                                            float("nan"),
                                                            df[f'mClsizeRecLayer{layer}'])
        df_origin = df.dropna() # drop rows with NaN values
        print('% of dropped candidates: ', round((1 - len(df)/len(df_origin))*100), 3)
    df_origin.to_parquet(f'{OutPutDir}parquet/df_MC.parquet.gzip', compression='gzip') # saving .parquet to speed up future loading process
else:
    print('\033[1mSkipping data preparation, loading .parquet file...' + '\033[0m')
    df_origin = pd.read_parquet(f'{OutPutDir}parquet/df_MC.parquet.gzip')
print('\033[1m' + 'Data preparation: complete' + '\033[0m')

df = df_origin
labels_true = df['mIsMIP'].values
labels_true = np.where(labels_true is True, 0, -1) # set the MIPs to 0 and Nucl. to -1
print('labels_true', labels_true)
df_mip_true = df.query('mIsMIP == True')
df_nucl_true = df.query('mIsMIP == False')
df = df.drop(columns=['mIsMIP']) # drop the column with the labels to avoid cheating
print('\nFinal DF')
print(df)
print('\nMIP DF')
print(df_mip_true)
print('\nNucl DF')
print(df_nucl_true)

#____________________________________________
# Data visualization
cols_to_keep = df.keys()
for col in cols_to_keep:
    if 'Gen' in col:
        cols_to_keep.drop(col)
outfile = TFile(f'{OutPutDir}/{mlmodel}_mftpid.root', 'recreate')
outfile.mkdir('DF_feature distr')
outfile.cd('DF_feature distr')

plot_utils.plot_distr([df_mip_true, df_nucl_true], cols_to_keep, 100, ['MIP', 'Nucl.'], figsize=(24, 14),
                      alpha=0.4, log=True, grid=False, density=True)
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
plt.savefig(f'{OutPutDir}/DF_feature distr.pdf')
plt.close('all')
for key in cols_to_keep:
    canvas = TCanvas(f"canvas{key}", "canvas", 800, 600)
    h_mip = plot_hist_from_df(df_mip_true, f'{key}', 0, 1, 100, title=f'{key}_mip')
    h_nucl = plot_hist_from_df(df_nucl_true, f'{key}', 0, 1, 100, title=f'{key}_nucl')
    SetObjectStyle(h_mip, linecolor=kAzure+4, markercolor=kAzure+4)
    SetObjectStyle(h_nucl, linecolor=kRed+1, markercolor=kRed+1)
    h_mip.DrawNormalized('hist same')
    h_nucl.DrawNormalized('hist same')
    canvas.Write()
    h_mip.Write()
    h_nucl.Write()
for key1 in cols_to_keep:
    for key2 in cols_to_keep:
        if key1 >= key2:
            continue
        h_mip = plot_hist2d_from_df(df_mip_true, f'{key1}', f'{key2}', 0, 1, 100, 0, 1, 100, title=f'{key1}_{key2}_mip')
        h_nucl = plot_hist2d_from_df(df_nucl_true, f'{key1}', f'{key2}', 0, 1, 100, 0, 1, 100, title=f'{key1}_{key2}_nucl')
        h_mip.Write()
        h_nucl.Write()
input('Press enter to continue...')

outfile.cd()

#____________________________________________
# ML model application
print('\033[1m' + 'ML model application...' + '\033[0m')
print (f'ML modele selected: {mlmodel}')
outfile.mkdir(f'{mlmodel}')
outfile.cd(f'{mlmodel}')
if mlmodel == 'DBSCAN':
    ######################################### DBSCAN #########################################
    # n_neighbors = 5 as kneighbors function returns distance of point to itself (i.e. first column will be zeros)
    print('\033[1m' + 'DBSCAN...' + '\033[0m')
    nbrs = NearestNeighbors(n_neighbors=42).fit(df)
    # Find the k-neighbors of a point
    neigh_dist, neigh_ind = nbrs.kneighbors(df)
    # sort the neighbor distances (lengths to points) in ascending order
    # axis = 0 represents sort along first axis i.e. sort along row
    sort_neigh_dist = np.sort(neigh_dist, axis=0)

    #____________________________________________
    # DBSCAN hyperparameters determination (https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html)
    # nearest neighbor distance plot for epsilon determination
    # (https://stackoverflow.com/questions/15050389/estimating-choosing-optimal-hyperparameters-for-dbscan/15063143#15063143)
    k_dist = sort_neigh_dist[:, 40] # 44 is the 45th element (index 44) as first element is 0
    fig = plt.figure()
    plt.plot(k_dist)
    plt.axhline(y=10, linewidth=1, linestyle='dashed', color='k')
    plt.ylabel("k-NN distance")
    plt.xlabel("Sorted observations (4th NN)")
    fig.savefig("dbs_map.pdf")
    # min_samples = 2*dim + 1 (https://stackoverflow.com/questions/15050389/estimating-choosing-optimal-hyperparameters-for-dbscan/15063143#15063143)
    min_samples = 2*len(df.columns) + 1

    db = DBSCAN(eps=10, min_samples=min_samples).fit(df) # fit the model
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    df_out = df
    df_out['mIsMIP'] = labels_true
    df_out['mIsMIP_pred'] = labels
    print(f'df_out: {df_out}')
    df_out.to_parquet(f'{OutPutDir}parquet/df_MCout.parquet.gzip', compression='gzip') # saving output df in .parquet
    
    # compare predicted labels with true labels distribution
    df_nucl_true = df_out.query('mIsMIP == -1', inplace=False)
    df_mip_true = df_out.query('mIsMIP == 0', inplace=False)
    df_nucl_pred = df_out.query('mIsMIP_pred == -1', inplace=False)
    df_mip_pred = df_out.query('mIsMIP_pred == 0', inplace=False)

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

    #____________________________________________
    # Precision and recall determination
    tp = len(df_out[(df_out['mIsMIP'] == -1) & (df_out['mIsMIP_pred'] == -1)]) # true positive (NUCL is positive)
    fn = len(df_out[(df_out['mIsMIP'] == -1) & (df_out['mIsMIP_pred'] == 0)]) # false negative
    tn = len(df_out[(df_out['mIsMIP'] == 0) & (df_out['mIsMIP_pred'] == 0)]) # true negative
    fp = len(df_out[(df_out['mIsMIP'] == 0) & (df_out['mIsMIP_pred'] == -1)]) # false positive
    print(f'tp: {tp}, fn: {fn}, tn: {tn}, fp: {fp}')
    # Precision = TP/(TP+FP)
    precision = tp/(tp+fp)
    # Recall = TP/(TP+FN)
    recall = tp/(tp+fn)

    print(f'\033[1mPrecision: {precision}\033[0m')
    print(f'\033[1mRecall: {recall}\033[0m')
    fig = plt.figure()
    plt.scatter(precision, recall, marker='o', color='r', s=100)
    plt.xlabel("Precision")
    plt.ylabel("Recall")
    plt.legend()
    fig.savefig("dbs_precision_recall.pdf")

    # confusion matrix
    fig = plt.figure()
    tn /= len(df_out[(df_out['mIsMIP'] == 0)])
    fn /= len(df_out[(df_out['mIsMIP'] == -1)])
    tp /= len(df_out[(df_out['mIsMIP'] == -1)])
    fp /= len(df_out[(df_out['mIsMIP'] == 0)])
    df_cm = pd.DataFrame([[tn, fp], [fn, tp]], index=['isMIP', 'isNucl'], columns=['isMipPred', 'isNuclPred'])
    plt.title("Confusion matrix")
    sn.heatmap(df_cm, annot=True, fmt='g', cmap='Blues')
    fig.savefig("dbs_confusion_matrix.pdf")

    fig = plt.figure()
    delta = np.abs(labels_true - labels)
    plt.hist(delta, bins=5, range=(-2.5, 2.5))
    fig.savefig("dbs_hist.pdf")

    fig = plt.figure()
    plt.hist(labels_true, bins=5, range=(-2.5, 2.5))
    plt.hist(labels, bins=5, range=(-2.5, 2.5), alpha=0.5)
    fig.savefig("dbs_labels.pdf")

    hNuclMeanClSize = TH1F("hNuclMeanClSize", "hNuclMeanClSize", 100, 0, 10)
    hNuclPredMeanClSize = TH1F("hNuclPredMeanClSize", "hNuclPredMeanClSize", 100, 0, 10)
    hMIPMeanClSize = TH1F("hMIPMeanClSize", "hMIPMeanClSize", 100, 0, 10)
    hMIPPredMeanClSize = TH1F("hMIPPredMeanClSize", "hMIPPredMeanClSize", 100, 0, 10)
    SetObjectStyle(hNuclMeanClSize, markerstyle=kFullCircle, markersize=1, markercolor=kRed+1,
                   fillstyle=3345, fillcolor=kRed+1, linewidth=1, linecolor=kRed+1, fillalpha=0.3)
    SetObjectStyle(hNuclPredMeanClSize, markerstyle=kOpenCircle, markersize=1, markercolor=kRed+1,
                   fillstyle=3354, fillcolor=kRed+1, linewidth=1, linecolor=kRed+1, fillalpha=0.3)
    SetObjectStyle(hMIPMeanClSize, markerstyle=kFullSquare, markersize=1, markercolor=kAzure+4,
                   fillstyle=3345, fillcolor=kAzure+4, linewidth=1, linecolor=kAzure+4, fillalpha=0.3)
    SetObjectStyle(hMIPPredMeanClSize, markerstyle=kOpenSquare, markersize=1, markercolor=kAzure+4,
                   fillstyle=3354, fillcolor=kAzure+4, linewidth=1, linecolor=kAzure+4, fillalpha=0.3)

    for i, (meanCl, isMip, isMipPred) in enumerate(zip(df_out['mClsizeMean'], df_out['mIsMIP'], df_out['mIsMIP_pred'])):
        if isMip == -1:
            hNuclMeanClSize.Fill(meanCl)
        if isMipPred == -1:
            hNuclPredMeanClSize.Fill(meanCl)
        if isMip == 0:
            hMIPMeanClSize.Fill(meanCl)
        if isMipPred == 0:
            hMIPPredMeanClSize.Fill(meanCl)
    canvas = TCanvas("canvas", "canvas", 800, 600)
    canvas.SetLogy()
    hNuclMeanClSize.SetStats(0)
    hNuclMeanClSize.Draw("hist e")
    hMIPMeanClSize.Draw("hist e same")
    hNuclPredMeanClSize.Draw("hist e same")
    hMIPPredMeanClSize.Draw("hist e same")
    canvas.SaveAs("dbs_meanClSize.pdf")

if mlmodel == 'AE':
    ############################ Autoencoder ############################
    # see: https://www.tensorflow.org/tutorials/generative/autoencoder#third_example_anomaly_detection
    df_mip_train = df_mip_true.sample(frac=1-test_frac, random_state=seed_split) # select the MIPs for training
    df_mip_test = df_mip_true.drop(df_mip_train.index) # select the MIPs for testing

    df_nucl_test = df_nucl_true # select the nucl for testing
    df_test = pd.concat([df_nucl_test, df_mip_test]) # concatenate the MIPs and nucl for testing
    df_out = df_test.copy() # copy the testing set for output
    df_test_labels = df_test['mIsMIP'] # get the labels for the testing set
    print(df_out['mIsMIP'])
    

    # drop the label
    df_mip_train.drop(columns=['mIsMIP'], inplace=True)
    df_mip_test.drop(columns=['mIsMIP'], inplace=True)
    df_nucl_test.drop(columns=['mIsMIP'], inplace=True)
    df_test.drop(columns=['mIsMIP'], inplace=True)

    df_mip_train = df_mip_train.to_numpy() # convert the training set to numpy
    df_mip_test = df_mip_test.to_numpy() # convert the testing set to numpy
    df_nucl_test = df_nucl_test.to_numpy() # convert the testing set to numpy
    df_test = df_test.to_numpy() # convert the testing set to numpy

    # Model preparation
    print('\033[1m' + 'Training the model...' + '\033[0m')
    model = AutoEncoder(df_mip_train.shape[1])
    model.compile(optimizer='adam', loss="mae")
    history = model.fit(df_mip_train, df_mip_train, epochs=50, batch_size=100,
                    validation_data=(df_mip_test, df_mip_test),
                    verbose=True, shuffle=True)
    MLplots._plot_loss_vs_epoch(history, OutPutDir) # plot the loss
    print('\r\033[1m' + 'Model training complete.' + '\033[0m')

    # Apply the model to the data
    print('\033[1m' + 'Applying the model to the data...' + '\033[0m')
    # mip
    encoder_out = model.encoder(df_mip_test).numpy()
    decoder_out = model.decoder(encoder_out).numpy()
    # nucl
    encoder_out_a = model.encoder(df_nucl_test).numpy()
    decoder_out_a = model.decoder(encoder_out_a).numpy()
    print('\r\033[1m' + 'Model application complete.' + '\033[0m')
    # loss & threshold
    ismip_pred = model.predict(df_mip_train)
    train_loss = tf.keras.losses.mae(ismip_pred, df_mip_train)
    test_loss = tf.keras.losses.mae(model.predict(df_test), df_test)
    df_out['mLoss'] = test_loss
    hist_train_loss = TH1F('hist_train_loss', 'hist_train_loss', 100, 0, 1)
    hist_test_loss = TH1F('hist_test_loss', 'hist_test_loss', 100, 0, 1)
    for i in range(len(train_loss)):
        hist_train_loss.Fill(train_loss[i])
    for i in range(len(test_loss)):
        hist_test_loss.Fill(test_loss[i])
    hist_train_loss.Write()
    hist_test_loss.Write()
    threshold = np.mean(train_loss) + 1*np.std(train_loss) # 1 std above the mean
    print(f'threshold: {threshold}')

    # predict
    preds = model.get_label_pred(df_test, threshold)
    preds = np.array(preds)
    df_out['mIsMIP_pred'] = preds
    print('\r\033[1m' + 'Model prediction complete.' + '\033[0m')
    print(f'preds: {len(preds)}')
    print(f'df_test_labels: {len(df_test_labels)}')
    # plot
    print(df_out)
    
    tp = len(df_out[(df_out['mIsMIP'] == False) & (df_out['mIsMIP_pred'] == False)])
    tn = len(df_out[(df_out['mIsMIP'] == True) & (df_out['mIsMIP_pred'] == True)])
    fp = len(df_out[(df_out['mIsMIP'] == True) & (df_out['mIsMIP_pred'] == False)])
    fn = len(df_out[(df_out['mIsMIP'] == False) & (df_out['mIsMIP_pred'] == True)])
    tn /= len(df_out[(df_out['mIsMIP'] == True)])
    fn /= len(df_out[(df_out['mIsMIP'] == False)])
    tp /= len(df_out[(df_out['mIsMIP'] == False)])
    fp /= len(df_out[(df_out['mIsMIP'] == True)])
    print(f'tp: {tp}, fn: {fn}, tn: {tn}, fp: {fp}')

    plot_confusion_matrix(tp, fp, fn, tn, OutPutDir,
                          labels=['MIP', 'Nucl'],
                          title='Confusion matrix AE')
    AutoEncoder.get_stats(df_test_labels, preds)

    hNuclMeanClSize = TH1F("hNuclMeanClSize", "hNuclMeanClSize", 100, 0, 10)
    hNuclPredMeanClSize = TH1F("hNuclPredMeanClSize", "hNuclPredMeanClSize", 100, 0, 10)
    hMIPMeanClSize = TH1F("hMIPMeanClSize", "hMIPMeanClSize", 100, 0, 10)
    hMIPPredMeanClSize = TH1F("hMIPPredMeanClSize", "hMIPPredMeanClSize", 100, 0, 10)
    SetObjectStyle(hNuclMeanClSize, markerstyle=kFullCircle, markersize=1, markercolor=kRed+1,
                   fillstyle=3345, fillcolor=kRed+1, linewidth=1, linecolor=kRed+1, fillalpha=0.3)
    SetObjectStyle(hNuclPredMeanClSize, markerstyle=kOpenCircle, markersize=1, markercolor=kRed+1,
                   fillstyle=3354, fillcolor=kRed+1, linewidth=1, linecolor=kRed+1, fillalpha=0.3)
    SetObjectStyle(hMIPMeanClSize, markerstyle=kFullSquare, markersize=1, markercolor=kAzure+4,
                   fillstyle=3345, fillcolor=kAzure+4, linewidth=1, linecolor=kAzure+4, fillalpha=0.3)
    SetObjectStyle(hMIPPredMeanClSize, markerstyle=kOpenSquare, markersize=1, markercolor=kAzure+4,
                   fillstyle=3354, fillcolor=kAzure+4, linewidth=1, linecolor=kAzure+4, fillalpha=0.3)

    for i, (meanCl, isMip, isMipPred) in enumerate(zip(df_out['mClsizeMean'], df_out['mIsMIP'], df_out['mIsMIP_pred'])):
        if isMipPred:
            hMIPPredMeanClSize.Fill(meanCl)
        else:
            hNuclPredMeanClSize.Fill(meanCl)
        if isMip == 0:
            hNuclMeanClSize.Fill(meanCl)
        else:
            hMIPMeanClSize.Fill(meanCl)
    canvas = TCanvas("canvas", "canvas", 800, 600)
    canvas.SetLogy()
    hNuclMeanClSize.SetStats(0)
    hNuclMeanClSize.Draw("hist e")
    hMIPMeanClSize.Draw("hist e same")
    hNuclPredMeanClSize.Draw("hist e same")
    hMIPPredMeanClSize.Draw("hist e same")
    canvas.Write()
    hNuclMeanClSize.Write()
    hNuclMeanClSize.Write()
    hMIPMeanClSize.Write()
    hNuclPredMeanClSize.Write()
    hMIPPredMeanClSize.Write()

    outfile.Close()
    df_out.to_parquet(f'{OutPutDir}/parquet/{mlmodel}_df_out.parquet.gzip', compression='gzip')
    


input('Press Enter to exit')
sys.exit()