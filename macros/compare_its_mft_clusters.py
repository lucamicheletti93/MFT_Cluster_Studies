# python script to compare the clusters of ITS and MFT
# run: python3 compare_its_mft_clusters.py
import pandas as pd
import uproot as ur
import numpy as np
import sys
from ROOT import TFile, gStyle, TH1, TH2F, TH1F, TH1D, TCanvas, TLegend, TGraph, TGraphErrors, kRed, kBlue, kCyan, kAzure, kSpring, kGreen, kOrange, kGray, kBlack, TGaxis,gPad, TLatex, kFullCircle, kOpenCross, kFullSquare, kFullCross # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetObjectStyle # pylint: disable=import-error,no-name-in-module


#_______________________________________________________________________________
# load ITS clusters
df_its = pd.read_parquet('/data/shared/ITS/ML/particles_pid_520143_itstpc.parquet') # taking V0s
df_its.query('p > 1.4', inplace=True) # taking only high momentum particles (p > 1.4 GeV/c)
df_its['delta_p'] = df_its['pITS'] - df_its['pTPC']
df_its.query('20 < rofBC < 500 and tpcITSchi2 < 5 and nClusTPC > 100 and -0.2 < delta_p < 0.2', inplace=True) # taking only good tracks'
for i in range(7):
        df_its[f'ClSizeL{i}'] = np.where(df_its[f'ClSizeL{i}'] < 0, 0., df_its[f'ClSizeL{i}'])
df_its['mean_cl_size'] = (df_its['ClSizeL0'] + df_its['ClSizeL1'] + df_its['ClSizeL2'] + df_its['ClSizeL3'] + df_its['ClSizeL4'] + df_its['ClSizeL5'] + df_its['ClSizeL6']) / 10
df_its['mean_patt_id'] = (df_its['PattIDL0'] + df_its['PattIDL1'] + df_its['PattIDL2'] + df_its['PattIDL3'] + df_its['PattIDL4'] + df_its['PattIDL5'] + df_its['PattIDL6']) / 10

histo_its_clsize = TH1F('histo_its_clsize', 'ITS cluster size', 50, 0, 10)
histo_its_pattid = TH1F('histo_its_pattid', 'ITS pattern ID', 50, 0, 10)
SetObjectStyle(histo_its_clsize, marker=kFullCircle, color=kRed+1,
               fill_color=kRed+1, fill_style=0, fillalpha=0.2)
SetObjectStyle(histo_its_pattid, marker=kFullCircle, color=kRed+1,
               fill_color=kRed+1, fill_style=0, fillalpha=0.2)
for i, (mean_cl_size, mean_patt_id) in enumerate(zip(df_its['mean_cl_size'], df_its['mean_patt_id'])):
        histo_its_clsize.Fill(mean_cl_size)
        histo_its_pattid.Fill(mean_patt_id)
print(df_its.head())

#_______________________________________________________________________________
# load MFT clusters
df_mft = ur.open('../output/MFTAssessmentMC.root')['treeTrackClusterSize'].arrays(library='pd')
df_mft.query('mIsMIP == True', inplace=True) # taking only MIPs (mIsMIP == True)
for i in range(10):
        df_mft[f'mClsizeRecLayer{i}'] = np.where(df_mft[f'mClsizeRecLayer{i}'] == -999, 0, df_mft[f'mClsizeRecLayer{i}'])
df_mft['mean_cl_size'] = (df_mft['mClsizeRecLayer0'] + df_mft['mClsizeRecLayer1'] + df_mft['mClsizeRecLayer2'] + df_mft['mClsizeRecLayer3'] + df_mft['mClsizeRecLayer4'] + df_mft['mClsizeRecLayer5'] + df_mft['mClsizeRecLayer6'] + df_mft['mClsizeRecLayer7'] + df_mft['mClsizeRecLayer8'] + df_mft['mClsizeRecLayer9']) / 10
df_mft['mean_patt_id'] = (df_mft['mPattIdRecLayer0'] + df_mft['mPattIdRecLayer1'] + df_mft['mPattIdRecLayer2'] + df_mft['mPattIdRecLayer3'] + df_mft['mPattIdRecLayer4'] + df_mft['mPattIdRecLayer5'] + df_mft['mPattIdRecLayer6'] + df_mft['mPattIdRecLayer7'] + df_mft['mPattIdRecLayer8'] + df_mft['mPattIdRecLayer9']) / 10

histo_mft_clsize = TH1F('histo_mft_clsize', 'MFT cluster size', 50, 0, 10)
histo_mft_pattid = TH1F('histo_mft_pattid', 'MFT pattern ID', 50, 0, 10)
SetObjectStyle(histo_mft_clsize, marker=kFullCircle, color=kAzure+4,
                fill_color=kAzure+1, fill_style=0, fillalpha=0.2)
SetObjectStyle(histo_mft_pattid, marker=kFullCircle, color=kAzure+1,
                fill_color=kAzure+1, fill_style=0, fillalpha=0.2)
for i, (mean_cl_size, mean_patt_id) in enumerate(zip(df_mft['mean_cl_size'], df_mft['mean_patt_id'])):
        histo_mft_clsize.Fill(mean_cl_size)
        histo_mft_pattid.Fill(mean_patt_id)
print(df_mft.head())

#_______________________________________________________________________________
# compare clusters
canvas = TCanvas('canvas', 'canvas', 800, 600)
leg = TLegend(0.6, 0.6, 0.8, 0.8)
leg.SetBorderSize(0)
leg.AddEntry(histo_its_clsize, 'ITS', 'lep')
leg.AddEntry(histo_mft_clsize, 'MFT', 'lep')
#canvas.Divide(2, 1)
canvas.cd(1).SetLogy()
histo_its_clsize.SetStats(0)
histo_its_clsize.DrawNormalized('hist e')
histo_mft_clsize.DrawNormalized('hist e same')
leg.Draw()
#canvas.cd(2)
#histo_its_pattid.Draw('hist e')
#histo_mft_pattid.Draw('hist e same')

canvas.SaveAs('clusters.png')