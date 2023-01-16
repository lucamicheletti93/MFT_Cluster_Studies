# python script to plot the efficiency and purity of the anomaly detection
# run: python plot_eff_purity.py

import pandas as pd
import uproot as ur
import numpy as np
import sys
from ROOT import TFile, gStyle, TH1, TH2F, TH1F, TH1D, TCanvas, TLegend, TGraph, TGraphErrors, kRed, kBlue, kCyan, kAzure, kSpring, kGreen, kOrange, kGray, kBlack, TGaxis,gPad, TLatex, kFullCircle, kOpenCross, kFullSquare, kFullCross # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetObjectStyle # pylint: disable=import-error,no-name-in-module
from utils.plot_library import plot_hist_from_df, plot_hist2d_from_df # pylint: disable=import-error,no-name-in-module

df = pd.read_parquet('parquet/AE_df_out.parquet.gzip')
print(df['mIsMIP_pred'], df['mIsMIP'])
thr = 0.02106 # threshold for the MIP prediction
df_nucl = df[df['mIsMIP'] == 0]
df_mip = df[df['mIsMIP'] == 1]

hist_nucl = plot_hist_from_df(df_nucl, 'mClsizeMean', 0, 1, 100, 'hist_nucl')
hist_mip = plot_hist_from_df(df_mip, 'mClsizeMean', 0, 1, 100, 'hist_mip')
hist_loss = plot_hist_from_df(df, 'mLoss', 0, 0.1, 100, 'hist_loss')
canv = TCanvas('canv', 'canv', 800, 600)
hist_nucl.SetLineColor(kRed)
hist_nucl.Draw()
hist_mip.SetLineColor(kBlue)
hist_mip.Draw('same')
hist_loss.SetLineColor(kBlack)
hist_loss.Draw('same')
input()

eff = []
purity = []
thr_high = thr
thr_low = thr
max_loss = max(df['mLoss'])
thrs = np.linspace(0, max_loss, 100)
print(max_loss)
for thr_high in thrs:
    print(f'Threshold: {thr_high}')
    isMIp_pred_thr = []
    for i, (mLoss) in enumerate(zip(df['mLoss'])):
        isMIp_pred_thr.append(1 if (mLoss[0] < thr_high) else 0)
    df['mIsMIP_pred_thr'] = isMIp_pred_thr
    tp = len(df[(df['mIsMIP_pred_thr'] == 0) & (df['mIsMIP'] == 0)])/len(df[df['mIsMIP'] == 0]) # true positive: predicted nucl, is nucl
    fp = len(df[(df['mIsMIP_pred_thr'] == 0) & (df['mIsMIP'] == 1)])/len(df[df['mIsMIP'] == 1]) # false positive: predicted nucl, is mip
    fn = len(df[(df['mIsMIP_pred_thr'] == 1) & (df['mIsMIP'] == 0)])/len(df[df['mIsMIP'] == 0]) # false negative: predicted mip, is nucl
    tn = len(df[(df['mIsMIP_pred_thr'] == 1) & (df['mIsMIP'] == 1)])/len(df[df['mIsMIP'] == 1]) # true negative: predicted mip, is mip

    eff.append(tp/(tp+fn))
    purity.append(tp/(tp+fp))

hist_eff_purity = TGraph(len(eff), np.array(eff), np.array(purity))
hist_eff_purity.SetTitle('Efficiency vs Purity')
hist_eff_purity.GetXaxis().SetTitle('Efficiency')
hist_eff_purity.GetYaxis().SetTitle('Purity')
SetObjectStyle(hist_eff_purity, markercolor=kBlack, markersize=0.5, markerstyle=kFullCircle)
canv = TCanvas('canv', 'canv', 800, 600)
hist_eff_purity.Draw('AP same')
input('Press enter to exit...')
