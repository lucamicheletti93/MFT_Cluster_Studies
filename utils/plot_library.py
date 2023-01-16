import matplotlib.pyplot as plt
import array as arr
import numpy as np
import os
import sys
import argparse
import ROOT
from os import path
from ROOT import TCanvas, TLatex, TF1, TFile, TPaveText, TMath, TH1F, TH2F, TString, TLegend, TRatioPlot, TGaxis
from ROOT import gROOT, gBenchmark, gPad, gStyle, kTRUE, kFALSE

def plot_hist_from_df(df, columns, min, max, bins, title='Histogram'):
    histo = TH1F(title, title, bins, min, max)
    for i in df[columns]:
        histo.Fill(i)
    return histo

def plot_hist2d_from_df(df, x, y, minx, maxx, binsx, miny, maxy,  binsy, title='Histogram'):
    histo = TH2F(title, title, binsx, minx, maxx, binsy, miny, maxy)
    for i, (ix, iy) in enumerate(zip(df[x], df[y])):
        histo.Fill(ix, iy)
    return histo