from ast import parse
from itertools import count
from time import process_time_ns
from tkinter import Canvas
from turtle import color, right
import matplotlib.pyplot as plt
from alive_progress import alive_bar
import time
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
def read(fInName, fOutName):
    '''
    Function to read the reduced tables and produce the output for ML training
    '''
    gStyle.SetOptStat(0)
    LoadStyle()
    ROOT.TGaxis.SetMaxDigits(3)

    #treeNames = ["treeTrackClusterSizeMC"]
    treeNames = ["treeTrackClusterSizeData"]
    isMC = False

    # init histograms
    nMFTLayers = 10
    histClsizeRecPerLayer = []
    histClsizeRecPerLayerPiKP = []
    histClsizeRecPerLayerHe3 = []
    histPattIdRecPerLayer = []
    for i in range(0, nMFTLayers):
        histClsizeRecPerLayer.append(ROOT.TH1I("histClsizeRecLayer{}".format(i), "MFT Layer {}".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayer[i], "MFT Layer {}".format(i), "cluster size", "Entries", ROOT.kBlack, 1, 3001, 0)

        histClsizeRecPerLayerPiKP.append(ROOT.TH1I("histClsizeRecLayerPiKP{}".format(i), "MFT Layer {} - #pi, K, p".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerPiKP[i], "MFT Layer {} - PiKP".format(i), "cluster size", "Entries", ROOT.kAzure+2, 1, 3001, 0.6)

        histClsizeRecPerLayerHe3.append(ROOT.TH1I("histClsizeRecLayerHe3{}".format(i), "MFT Layer {} - He3".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerHe3[i], "MFT Layer {} - He3".format(i), "cluster size", "Entries", ROOT.kRed+1, 1, 3001, 0.6)

        histPattIdRecPerLayer.append(ROOT.TH1I("histPattIdRecPerLayer{}".format(i), "MFT Layer {}".format(i), 5000, -0.5, 4999.5))
        SetHistStyle2(histPattIdRecPerLayer[i], "MFT Layer {}".format(i), "pattern ID", "Entries", ROOT.kBlack, 1, 3001, 0)

    histClsizeRec = ROOT.TH1I("histClsizeRec", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRec, "", "cluster size", "Entries", ROOT.kBlack, 2, 3001, 0)

    histClsizeRecPiKP = ROOT.TH1I("histClsizeRecPiKP", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRecPiKP, "", "cluster size", "Entries", ROOT.kAzure+2, 2, 3001, 0.6)

    histClsizeRecHe3 = ROOT.TH1I("histClsizeRecHe3", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRecHe3, "", "cluster size", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histPattIdRec = ROOT.TH1I("histPattIdRec", "", 100, -0.5, 99.5)
    SetHistStyle2(histPattIdRec, "", "pattern ID", "Entries", ROOT.kBlack, 2, 3001, 0)

    histMeanClsizePerTrackRec = ROOT.TH1F("histMeanClsizePerTrackRec", "", 80, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRec, "", "<cluster size>", "Entries", ROOT.kBlack, 2, 3001, 0)

    histMeanClsizePerTrackRecPiKP = ROOT.TH1F("histMeanClsizePerTrackRecPiKP", "", 80, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecPiKP, "", "<cluster size>", "Entries", ROOT.kAzure+2, 2, 3005, 0.6)

    histMeanClsizePerTrackRecPi = ROOT.TH1F("histMeanClsizePerTrackRecPi", "", 80, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecPi, "", "<cluster size>", "Entries", ROOT.kRed+2, 2, 3005, 0.3)

    histMeanClsizePerTrackRecK = ROOT.TH1F("histMeanClsizePerTrackRecK", "", 80, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecK, "", "<cluster size>", "Entries", ROOT.kGreen+2, 2, 3005, 0.3)
    
    histMeanClsizePerTrackRecP = ROOT.TH1F("histMeanClsizePerTrackRecP", "", 80, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecP, "", "<cluster size>", "Entries", ROOT.kBlue+2, 2, 3005, 0.3)

    histMeanClsizePerTrackRecHe3 = ROOT.TH1F("histMeanClsizePerTrackRecHe3", "", 80, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecHe3, "", "<cluster size>", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    effChargeMatch = ROOT.TEfficiency("effChargeMatch", "Charge Match;p_t [GeV];#epsilon", 20, 0, 10)
    
    for treeName in treeNames:
        fIn = TFile.Open(fInName)
        treeReaderInput = TTreeReader(treeName, fIn)
        mChargeRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeRec")
        #mClsizeRec = ROOT.TTreeReaderArray(ROOT.int)(treeReaderInput, "mClsizeRec")
        mMeanClsizePerTrackRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mMeanClsizePerTrackRec")
        if isMC:
            mPtGen = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPtGen")
            mChargeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeGen")
            mPdgCodeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPdgCodeGen")
        # Clsize per layer
        mClsizeRecLayer0 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer0")
        mClsizeRecLayer1 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer1")
        mClsizeRecLayer2 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer2")
        mClsizeRecLayer3 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer3")
        mClsizeRecLayer4 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer4")
        mClsizeRecLayer5 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer5")
        mClsizeRecLayer6 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer6")
        mClsizeRecLayer7 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer7")
        mClsizeRecLayer8 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer8")
        mClsizeRecLayer9 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRecLayer9")
        # Pattern ID per layer
        mPattIdRecLayer0 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer0")
        mPattIdRecLayer1 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer1")
        mPattIdRecLayer2 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer2")
        mPattIdRecLayer3 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer3")
        mPattIdRecLayer4 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer4")
        mPattIdRecLayer5 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer5")
        mPattIdRecLayer6 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer6")
        mPattIdRecLayer7 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer7")
        mPattIdRecLayer8 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer8")
        mPattIdRecLayer9 = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRecLayer9")

        nEv = treeReaderInput.GetEntries()
        evCounter = 0
        print("n. Events = ", nEv)
        print("Processing {}...".format(fIn.GetName()))
        with alive_bar(nEv, force_tty=True) as bar:
            while treeReaderInput.Next():
                #if evCounter > 1000:
                    #break
                evCounter = evCounter + 1
                time.sleep(0.0000001)
                bar()

                histMeanClsizePerTrackRec.Fill(mMeanClsizePerTrackRec.Get()[0])
                if isMC:
                    deltaCharge = mChargeRec.Get()[0] - mChargeGen.Get()[0]
                    effChargeMatch.Fill(not deltaCharge, mPtGen.Get()[0])
                    if abs(mPdgCodeGen.Get()[0]) == 1000020030:
                        histMeanClsizePerTrackRecHe3.Fill(mMeanClsizePerTrackRec.Get()[0])
                    if abs(mPdgCodeGen.Get()[0]) == 211 or abs(mPdgCodeGen.Get()[0]) == 321 or abs(mPdgCodeGen.Get()[0]) == 2212:
                        histMeanClsizePerTrackRecPiKP.Fill(mMeanClsizePerTrackRec.Get()[0])
                    if abs(mPdgCodeGen.Get()[0]) == 211:
                        histMeanClsizePerTrackRecPi.Fill(mMeanClsizePerTrackRec.Get()[0])
                    if abs(mPdgCodeGen.Get()[0]) == 321:
                        histMeanClsizePerTrackRecK.Fill(mMeanClsizePerTrackRec.Get()[0])
                    if abs(mPdgCodeGen.Get()[0]) == 2212:
                        histMeanClsizePerTrackRecP.Fill(mMeanClsizePerTrackRec.Get()[0])

                mClsizeRec = []
                mClsizeRec.append(mClsizeRecLayer0.Get()[0])
                mClsizeRec.append(mClsizeRecLayer1.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer2.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer3.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer4.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer5.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer6.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer7.Get()[0])  
                mClsizeRec.append(mClsizeRecLayer8.Get()[0])
                mClsizeRec.append(mClsizeRecLayer9.Get()[0])

                mPattIdRec = []
                mPattIdRec.append(mPattIdRecLayer0.Get()[0])
                mPattIdRec.append(mPattIdRecLayer1.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer2.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer3.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer4.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer5.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer6.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer7.Get()[0])  
                mPattIdRec.append(mPattIdRecLayer8.Get()[0])
                mPattIdRec.append(mPattIdRecLayer9.Get()[0]) 

                for layer in range(0, 10):
                    nPixel = mClsizeRec[layer]
                    pattId = mPattIdRec[layer]
                    if nPixel <= 0:
                        continue
                    histClsizeRec.Fill(nPixel)
                    histClsizeRecPerLayer[layer].Fill(nPixel)

                    if pattId <= 0:
                        continue
                    histPattIdRec.Fill(pattId)
                    histPattIdRecPerLayer[layer].Fill(pattId)

                    if isMC:
                        if abs(mPdgCodeGen.Get()[0]) == 1000020030:
                            histClsizeRecHe3.Fill(nPixel)
                            histClsizeRecPerLayerHe3[layer].Fill(nPixel)
                        if abs(mPdgCodeGen.Get()[0]) == 211 or abs(mPdgCodeGen.Get()[0]) == 321 or abs(mPdgCodeGen.Get()[0]) == 2212:
                            histClsizeRecPiKP.Fill(nPixel)
                            histClsizeRecPerLayerPiKP[layer].Fill(nPixel)
    
    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.050)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    fOut = TFile.Open(fOutName, "RECREATE")
    fOut.cd()
    for i in range(0, 10):
        histClsizeRecPerLayer[i].Write()
        histPattIdRecPerLayer[i].Write()
    histClsizeRec.Write()
    histPattIdRec.Write()
    histMeanClsizePerTrackRec.Write()

    if isMC:
        legend = ROOT.TLegend(0.70, 0.40, 0.89, 0.73, " ", "brNDC")
        SetLegend(legend)
        legend.AddEntry(histClsizeRec, "All", "F")
        legend.AddEntry(histClsizeRecHe3, "^{3}He / ^{3}#bar{He}", "F")
        legend.AddEntry(histClsizeRecPiKP, "#pi, K, p", "F")

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
        canvasClsizeRec.SaveAs("figures/MFT_clsize_per_track.pdf")

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
        canvasClsizeRecPerLayer.SaveAs("figures/MFT_clsize_mean_per_track_per_layer.pdf")

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
        canvasMeanClsizePerTrackRec.SaveAs("figures/MFT_clsize_mean_per_track.pdf")

        # Normalize histograms
        histMeanClsizePerTrackRecPi.Scale(1. / histMeanClsizePerTrackRecPi.Integral())
        histMeanClsizePerTrackRecK.Scale(1. / histMeanClsizePerTrackRecK.Integral())
        histMeanClsizePerTrackRecP.Scale(1. / histMeanClsizePerTrackRecP.Integral())

        legendLF = ROOT.TLegend(0.70, 0.40, 0.89, 0.73, " ", "brNDC")
        SetLegend(legendLF)
        legendLF.AddEntry(histMeanClsizePerTrackRecPi, "#pi", "PE")
        legendLF.AddEntry(histMeanClsizePerTrackRecK, "K", "PE")
        legendLF.AddEntry(histMeanClsizePerTrackRecP, "p", "PE")

        canvasMeanClsizePerTrackRecLF = ROOT.TCanvas("canvasMeanClsizePerTrackRecLF", "canvasMeanClsizePerTrackRecLF", 800, 600)
        histMeanClsizePerTrackRecPi.Draw("EP")
        histMeanClsizePerTrackRecK.Draw("EPsame")
        histMeanClsizePerTrackRecP.Draw("EPsame")
        legendLF.Draw("same")
        latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
        latexTitle.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
        latexTitle.DrawLatex(0.6, 0.73, "Run3 simulation")
        canvasMeanClsizePerTrackRecLF.Update()
        canvasMeanClsizePerTrackRecLF.SaveAs("figures/MFT_clsize_mean_per_track_LF.pdf")

        # Charge efficiency reconstruction
        histTotalChargeMatch = effChargeMatch.GetCopyTotalHisto()
        histPassedChargeMatch = effChargeMatch.GetCopyPassedHisto()
        histEffChargeMatch = ROOT.TH1D("histEffChargeMatch", "", 20, 0, 10)
        histEffChargeMatch.Divide(histPassedChargeMatch, histTotalChargeMatch, 1, 1, "B")
        histEffChargeMatch.SetMarkerColor(ROOT.kRed+1)
        histEffChargeMatch.SetLineColor(ROOT.kRed+1)
        histEffChargeMatch.SetLineWidth(3)
        histEffChargeMatch.GetXaxis().SetRangeUser(0, 5)
        histEffChargeMatch.GetXaxis().SetTitle("#it{p}_{T}^{Gen} (GeV/#it{c})")
        histEffChargeMatch.GetYaxis().SetRangeUser(0.5, 1.2)
        histEffChargeMatch.GetYaxis().SetTitle("Efficiency")

        lineUnity = ROOT.TLine(0, 1, 5, 1)
        lineUnity.SetLineWidth(2)
        lineUnity.SetLineColor(ROOT.kGray+1)
        lineUnity.SetLineStyle(ROOT.kDashed)

        canvasChargeMatch = ROOT.TCanvas("canvasChargeMatch", "", 800, 600)
        histEffChargeMatch.Draw("E")
        lineUnity.Draw("same")
        latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
        latexTitle.DrawLatex(0.6, 0.80, "Run3 simulation")
        canvasChargeMatch.Update()
        canvasChargeMatch.SaveAs("figures/MFT_charge_match_eff.pdf")

    

    # Cluster size for all reconstructed tracks per layer
    canvasClsizeRecPerLayerData = ROOT.TCanvas("canvasClsizeRecPerLayerData", "canvasClsizeRecPerLayerData", 3000, 1200)
    canvasClsizeRecPerLayerData.Divide(5, 2)
    for i in range(0, 10):
        canvasClsizeRecPerLayerData.cd(i+1)
        gPad.SetLogx(1)
        gPad.SetLogy(1)
        histClsizeRecPerLayer[i].Draw("HE")
    canvasClsizeRecPerLayerData.Update()
    #canvasClsizeRecPerLayerData.SaveAs("figures/MFT_clsize_mean_per_track_per_layer_data.pdf")

    # Mean cluster size for all reconstructed tracks
    canvasMeanClsizePerTrackRecData = ROOT.TCanvas("canvasMeanClsizePerTrackRecData", "canvasMeanClsizePerTrackRecData", 800, 600)
    histMeanClsizePerTrackRec.Draw("HE")
    latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitle.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    latexTitle.DrawLatex(0.6, 0.73, "Run3 simulation")
    canvasMeanClsizePerTrackRecData.Update()
    #canvasMeanClsizePerTrackRecData.SaveAs("figures/MFT_clsize_mean_per_track_data.pdf")

    exit()

def compare_data_mc():
    gStyle.SetOptStat(0)
    LoadStyle()
    ROOT.TGaxis.SetMaxDigits(3)

    fInData = TFile.Open("../output/MFTAssessmentData.root")
    histClsizeData = fInData.Get("mMFTClusterSize")
    histClsizeData.SetMarkerStyle(20)
    histClsizeData.SetMarkerColor(ROOT.kBlack)
    histClsizeData.SetLineColor(ROOT.kBlack)
    histClsizeData.Scale(1. / histClsizeData.Integral())

    histMeanClsizePerTrackData = fInData.Get("mMFTTrackMeanClusterSize")
    histMeanClsizePerTrackData.SetMarkerStyle(20)
    histMeanClsizePerTrackData.SetMarkerColor(ROOT.kBlack)
    histMeanClsizePerTrackData.SetLineColor(ROOT.kBlack)
    histMeanClsizePerTrackData.Rebin(10)
    histMeanClsizePerTrackData.Scale(1. / histClsizeData.Integral())

    fInMC = TFile.Open("../output/MFTAssessmentMC.root")
    histClsizeMC = fInMC.Get("mMFTClusterSize")
    SetHistStyle2(histClsizeMC, "", "cluster size", "Entries", ROOT.kRed+1, 2, 3005, 0.5)
    histClsizeMC.Scale(1. / histClsizeMC.Integral())

    legend = ROOT.TLegend(0.70, 0.65, 0.89, 0.85, " ", "brNDC")
    SetLegend(legend)
    legend.AddEntry(histClsizeData, "Data", "EP")
    legend.AddEntry(histClsizeMC, "MC", "F")

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.050)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    canvasClsizeComp = ROOT.TCanvas("canvasClsizeComp", "canvasClsizeComp", 800, 600)
    gPad.SetLogy(1)
    histClsizeMC.Draw("H")
    histClsizeData.Draw("EPsame")
    legend.Draw("same")
    latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitle.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    canvasClsizeComp.Update()

    canvasClsizePerTrackComp = ROOT.TCanvas("canvasClsizePerTrackComp", "canvasClsizePerTrackComp", 800, 600)
    gPad.SetLogy(1)
    histMeanClsizePerTrackData.Draw("EPsame")
    legend.Draw("same")
    latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitle.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    canvasClsizePerTrackComp.Update()

    input()
    canvasClsizeComp.SaveAs("figures/Data_MC_clsize_comparison.pdf")

    exit()

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read", help="Read the input and produce output histograms", action="store_true")
    parser.add_argument("--compare_data_mc", help="Compare Data and MC distributions", action="store_true")
    args = parser.parse_args()

    if args.read:
        #read("../output/MFTAssessment_test.root")
        #read("../output/MFTAssessmentMC.root")
        #read("../output/MFTAssessmentData_LHC22n.root", "../output/MFTAssessmentData_LHC22n_skimmed.root")
        read("../output/MFTAssessmentData_LHC22n_526195_0420.root", "../output/MFTAssessmentData_LHC22n_526195_0420_skimmed.root")
    if args.compare_data_mc:
        compare_data_mc()

main()