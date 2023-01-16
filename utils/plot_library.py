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

def SetLatex(latex):
    latex.SetTextSize(0.035)
    latex.SetNDC()
    latex.SetTextFont(42)

def SetLegend(legend):
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetFillStyle(1)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)

def SetHistStyle(hist):
    hist.SetTitle("")
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitle("Entries")
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)

def SetHistStyle2(hist, title, xaxis, yaxis, color, linewidth, fillstyle, trasparency):
    hist.SetTitle(title)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitle(xaxis)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitle(yaxis)
    hist.SetLineColor(color)
    hist.SetLineWidth(linewidth)
    #hist.SetFillStyle(fillstyle)
    hist.SetFillColorAlpha(color, trasparency)

def SetHistStyle3(hist, title, xaxis, yaxis, color, linewidth, markerstyle, markersize):
    hist.SetTitle(title)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetTitle(xaxis)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitle(yaxis)
    hist.SetMarkerStyle(markerstyle)
    hist.SetMarkerSize(markersize)
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.SetLineWidth(linewidth)

def LoadStyle():
    gStyle.SetOptStat(0)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadBottomMargin(0.15)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetEndErrorSize(0.0)
    gStyle.SetTitleSize(0.05,"X")
    gStyle.SetTitleSize(0.045,"Y")
    gStyle.SetLabelSize(0.045,"X")
    gStyle.SetLabelSize(0.045,"Y")
    gStyle.SetTitleOffset(1.2,"X")
    gStyle.SetTitleOffset(1.35,"Y")