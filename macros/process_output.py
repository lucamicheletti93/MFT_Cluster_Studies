from ast import parse
from itertools import count
from time import process_time_ns
from tkinter import Canvas
from turtle import color, right
import matplotlib.pyplot as plt
import array as arr
import numpy as np
import pandas as pd
from pandas.plotting import scatter_matrix
import os
import sys
import argparse
import ROOT
from os import path
from ROOT import TCanvas, TFile, TTree, TTreeReader, TH1F, TH2F, TLatex, TLegend, TF1
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack
from ctypes import cdll
sys.path.append('../utils')
from plot_library import LoadStyle, SetLegend, SetHistStyle, SetLatex, SetHistStyle2

###
def read_sim(fInName):
    '''
    Function to read the reduced tables and produce the output for ML training
    '''
    gStyle.SetOptStat(0)
    LoadStyle()
    ROOT.TGaxis.SetMaxDigits(3)

    treeNames = ["treeTrackClusterSize"]

    # init histograms
    nMFTLayers = 10
    histClsizeRecPerLayer = []
    histClsizeRecPerLayerPiKP = []
    histClsizeRecPerLayerHe3 = []
    for i in range(0, nMFTLayers):
        histClsizeRecPerLayer.append(ROOT.TH1I("histClsizeRecLayer{}".format(i), "MFT Layer {}".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayer[i], "MFT Layer {}".format(i), "cluster size", "Entries", ROOT.kBlack, 1, 3001, 0)

        histClsizeRecPerLayerPiKP.append(ROOT.TH1I("histClsizeRecLayerPiKP{}".format(i), "MFT Layer {} - #pi, K, p".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerPiKP[i], "MFT Layer {} - PiKP".format(i), "cluster size", "Entries", ROOT.kAzure+2, 1, 3001, 0.6)

        histClsizeRecPerLayerHe3.append(ROOT.TH1I("histClsizeRecLayerHe3{}".format(i), "MFT Layer {} - He3".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerHe3[i], "MFT Layer {} - He3".format(i), "cluster size", "Entries", ROOT.kRed+1, 1, 3001, 0.6)

    histClsizeRec = ROOT.TH1I("histClsizeRec", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRec, "", "cluster size", "Entries", ROOT.kBlack, 2, 3001, 0)

    histClsizeRecPiKP = ROOT.TH1I("histClsizeRecPiKP", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRecPiKP, "", "cluster size", "Entries", ROOT.kAzure+2, 2, 3001, 0.6)

    histClsizeRecHe3 = ROOT.TH1I("histClsizeRecHe3", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRecHe3, "", "cluster size", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histMeanClsizePerTrackRec = ROOT.TH1F("histMeanClsizePerTrackRec", "", 40, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRec, "", "<cluster size>", "Entries", ROOT.kBlack, 2, 3001, 0)

    histMeanClsizePerTrackRecPiKP = ROOT.TH1F("histMeanClsizePerTrackRecPiKP", "", 40, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecPiKP, "", "<cluster size>", "Entries", ROOT.kAzure+2, 2, 3005, 0.6)

    histMeanClsizePerTrackRecHe3 = ROOT.TH1F("histMeanClsizePerTrackRecHe3", "", 40, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecHe3, "", "<cluster size>", "Entries", ROOT.kRed+1, 2, 3001, 0.6)
    
    for treeName in treeNames:
        fIn = TFile.Open(fInName)
        treeReaderInput = TTreeReader(treeName, fIn)
        mChargeRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeRec")
        mClsizeRec = ROOT.TTreeReaderArray(ROOT.int)(treeReaderInput, "mClsizeRec")
        mLayerRec = ROOT.TTreeReaderArray(ROOT.int)(treeReaderInput, "mLayerRec")
        mMeanClsizePerTrackRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mMeanClsizePerTrackRec")
        mPdgCodeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPdgCodeGen")

        print("Processing {}...".format(fIn.GetName()))
        while treeReaderInput.Next():
            histMeanClsizePerTrackRec.Fill(mMeanClsizePerTrackRec.Get()[0])
            if abs(mPdgCodeGen.Get()[0]) == 1000020030:
                histMeanClsizePerTrackRecHe3.Fill(mMeanClsizePerTrackRec.Get()[0])
            if abs(mPdgCodeGen.Get()[0]) == 211 or abs(mPdgCodeGen.Get()[0]) == 321 or abs(mPdgCodeGen.Get()[0]) == 2212:
                histMeanClsizePerTrackRecPiKP.Fill(mMeanClsizePerTrackRec.Get()[0])
            counter = 0
            for nCl in mClsizeRec:
                histClsizeRec.Fill(nCl)
                layer = mLayerRec[counter]
                histClsizeRecPerLayer[layer].Fill(nCl)
                if abs(mPdgCodeGen.Get()[0]) == 1000020030:
                    histClsizeRecHe3.Fill(nCl)
                    histClsizeRecPerLayerHe3[layer].Fill(nCl)
                if abs(mPdgCodeGen.Get()[0]) == 211 or abs(mPdgCodeGen.Get()[0]) == 321 or abs(mPdgCodeGen.Get()[0]) == 2212:
                    histClsizeRecPiKP.Fill(nCl)
                    histClsizeRecPerLayerPiKP[layer].Fill(nCl)
                counter = counter + 1
    
    legend = ROOT.TLegend(0.70, 0.40, 0.89, 0.73, " ", "brNDC")
    SetLegend(legend)
    legend.AddEntry(histClsizeRec, "All", "F")
    legend.AddEntry(histClsizeRecHe3, "^{3}He / ^{3}#bar{He}", "F")
    legend.AddEntry(histClsizeRecPiKP, "#pi, K, p", "F")

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.050)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)
    
    # Cluster size for all reconstructed tracks
    canvasClsizeRec = ROOT.TCanvas("canvasClsizeRec", "canvasClsizeRec", 800, 600)
    gPad.SetLogy(1)
    histClsizeRec.Draw("HE")
    histClsizeRecPiKP.Draw("HEsame")
    histClsizeRecHe3.Draw("HEsame")
    legend.Draw("same")
    latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitle.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    latexTitle.DrawLatex(0.6, 0.73, "Run3 simulation")
    canvasClsizeRec.Update()
    canvasClsizeRec.SaveAs("MFT_clsize_per_track.pdf")

    # Cluster size for all reconstructed tracks per layer
    canvasClsizeRecPerLayer = ROOT.TCanvas("canvasClsizeRecPerLayer", "canvasClsizeRecPerLayer", 3000, 1200)
    canvasClsizeRecPerLayer.Divide(5, 2)
    for i in range(0, 10):
        canvasClsizeRecPerLayer.cd(i+1)
        gPad.SetLogx(1)
        gPad.SetLogy(1)
        histClsizeRecPerLayer[i].Draw("HE")
        histClsizeRecPerLayerPiKP[i].Draw("HEsame")
        histClsizeRecPerLayerHe3[i].Draw("HEsame")
        if i == 0:
            legend.Draw("same")
        latexTitle.DrawLatex(0.5, 0.85, "<#pi, K, p> = %3.2f" % (histClsizeRecPerLayerPiKP[i].GetMean()))
        latexTitle.DrawLatex(0.5, 0.79, "<^{3}He> = %3.2f" % (histClsizeRecPerLayerHe3[i].GetMean()))
    canvasClsizeRecPerLayer.Update()
    canvasClsizeRecPerLayer.SaveAs("MFT_clsize_mean_per_track_per_layer.pdf")

    # Mean cluster size for all reconstructed tracks
    canvasMeanClsizePerTrackRec = ROOT.TCanvas("canvasMeanClsizePerTrackRec", "canvasMeanClsizePerTrackRec", 800, 600)
    histMeanClsizePerTrackRec.Draw("HE")
    histMeanClsizePerTrackRecPiKP.Draw("HEsame")
    histMeanClsizePerTrackRecHe3.Draw("HEsame")
    legend.Draw("same")
    latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitle.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    latexTitle.DrawLatex(0.6, 0.73, "Run3 simulation")
    canvasMeanClsizePerTrackRec.Update()
    canvasMeanClsizePerTrackRec.SaveAs("MFT_clsize_mean_per_track.pdf")

    #input()

    canvasClsizeRec.Close()
    
    exit()

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read_sim", help="Read the simulation and produce output histograms", action="store_true")
    args = parser.parse_args()

    if args.read_sim:
        read_sim("../output/MFTAssessmentMC.root")

main()