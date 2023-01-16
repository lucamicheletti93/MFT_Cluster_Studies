from ast import parse
from itertools import count
from time import process_time_ns
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
from ROOT import TCanvas, TFile, TTree, TTreeReader, TH1F, TH2F, TLatex, TLegend, TF1, TLine
from ROOT import gStyle, gPad, kRed, kBlue, kGreen, kOrange, kAzure, kBlack
from ctypes import cdll
sys.path.append('../utils')
from plot_library import LoadStyle, SetLegend, SetHistStyle, SetLatex, SetHistStyle2, SetHistStyle3

###
def read(fInName, fOutName):
    '''
    Function to read the reduced tables and produce QC plots
    '''
    gStyle.SetOptStat(0)
    LoadStyle()
    ROOT.TGaxis.SetMaxDigits(3)

    treeNames = ["treeTrackClusterSize"]
    #treeNames = ["treeTrackClusterSizeData"]
    isMC = True

    # init histograms
    nMFTLayers = 10
    histClsizeRecPerLayer = []
    histClsizeRecPerLayerMIP = []
    histClsizeRecPerLayerPiKP = []
    histClsizeRecPerLayerHe3 = []
    histPattIdRecPerLayer = []
    histPattIdRecPerLayerMIP = []
    histPattIdRecPerLayerPiKP = []
    for i in range(0, nMFTLayers):
        histClsizeRecPerLayer.append(ROOT.TH1F("histClsizeRecLayer{}".format(i), "MFT Layer {}".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayer[i], "MFT Layer {}".format(i), "cluster size", "Entries", ROOT.kBlack, 1, 3001, 0)

        histClsizeRecPerLayerMIP.append(ROOT.TH1F("histClsizeRecLayer{}_MIP".format(i), "MFT Layer {}".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerMIP[i], "MFT Layer {}".format(i), "cluster size", "Entries", ROOT.kBlack, 1, 3001, 0)

        histClsizeRecPerLayerPiKP.append(ROOT.TH1F("histClsizeRecLayer{}_PiKP".format(i), "MFT Layer {} - #pi, K, p".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerPiKP[i], "MFT Layer {} - PiKP".format(i), "cluster size", "Entries", ROOT.kAzure+2, 1, 3001, 0.6)

        histClsizeRecPerLayerHe3.append(ROOT.TH1F("histClsizeRecLayer{}_He3".format(i), "MFT Layer {} - He3".format(i), 1000, -0.5, 999.5))
        SetHistStyle2(histClsizeRecPerLayerHe3[i], "MFT Layer {} - He3".format(i), "cluster size", "Entries", ROOT.kRed+1, 1, 3001, 0.6)

        histPattIdRecPerLayer.append(ROOT.TH1F("histPattIdRecLayer{}".format(i), "MFT Layer {}".format(i), 10000, -0.5, 9999.5))
        SetHistStyle2(histPattIdRecPerLayer[i], "MFT Layer {}".format(i), "pattern ID", "Entries", ROOT.kBlack, 1, 3001, 0)

        histPattIdRecPerLayerMIP.append(ROOT.TH1F("histPattIdRecLayer{}_MIP".format(i), "MFT Layer {}".format(i), 10000, -0.5, 9999.5))
        SetHistStyle2(histPattIdRecPerLayerMIP[i], "MFT Layer {}".format(i), "pattern ID", "Entries", ROOT.kBlack, 1, 3001, 0)

        histPattIdRecPerLayerPiKP.append(ROOT.TH1F("histPattIdRecLayer{}_PiKP".format(i), "MFT Layer {} - #pi, K, p".format(i), 10000, -0.5, 9999.5))
        SetHistStyle2(histPattIdRecPerLayerPiKP[i], "MFT Layer {} - PiKP".format(i), "pattern ID", "Entries", ROOT.kBlack, 1, 3001, 0)

    histClsizeRec = ROOT.TH1F("histClsizeRec", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRec, "", "cluster size", "Entries", ROOT.kBlack, 2, 3001, 0)

    histClsizeRecPiKP = ROOT.TH1F("histClsizeRecPiKP", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRecPiKP, "", "cluster size", "Entries", ROOT.kAzure+2, 2, 3001, 0.6)

    histClsizeRecHe3 = ROOT.TH1F("histClsizeRecHe3", "", 100, -0.5, 99.5)
    SetHistStyle2(histClsizeRecHe3, "", "cluster size", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histPattIdRec = ROOT.TH1F("histPattIdRec", "", 1000, -0.5, 999.5)
    SetHistStyle2(histPattIdRec, "", "pattern ID", "Entries", ROOT.kBlack, 2, 3001, 0)

    histMeanClsizePerTrackRec = ROOT.TH1F("histMeanClsizePerTrackRec", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRec, "", "<cluster size>", "Entries", ROOT.kBlack, 2, 3001, 0)

    histEtaRec = ROOT.TH1F("histEtaRec", "", 200, -10, 10)
    SetHistStyle2(histEtaRec, "", "#eta", "Entries", ROOT.kBlack, 2, 3001, 0)

    histPhiRec = ROOT.TH1F("histPhiRec", "", 100, 0, 2*ROOT.TMath.Pi())
    SetHistStyle2(histPhiRec, "", "#phi", "Entries", ROOT.kBlack, 2, 3001, 0)

    histMeanClsizePerTrackRecPiKP = ROOT.TH1F("histMeanClsizePerTrackRecPiKP", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecPiKP, "", "<cluster size>", "Entries", ROOT.kAzure+2, 2, 3005, 0.6)

    histMeanClsizePerTrackRecPi = ROOT.TH1F("histMeanClsizePerTrackRecPi", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecPi, "", "<cluster size>", "Entries", ROOT.kRed+2, 2, 3005, 0.3)

    histMeanClsizePerTrackRecK = ROOT.TH1F("histMeanClsizePerTrackRecK", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecK, "", "<cluster size>", "Entries", ROOT.kGreen+2, 2, 3005, 0.3)
    
    histMeanClsizePerTrackRecP = ROOT.TH1F("histMeanClsizePerTrackRecP", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecP, "", "<cluster size>", "Entries", ROOT.kBlue+2, 2, 3005, 0.3)

    histMeanClsizePerTrackRecHe3 = ROOT.TH1F("histMeanClsizePerTrackRecHe3", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecHe3, "", "<cluster size>", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histMeanClsizePerTrackRecHe3 = ROOT.TH1F("histMeanClsizePerTrackRecHe3", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecHe3, "", "<cluster size>", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histMeanClsizePerTrackRecMIP = ROOT.TH1F("histMeanClsizePerTrackRecMIP", "", 800, 0., 8.)
    SetHistStyle2(histMeanClsizePerTrackRecMIP, "", "<cluster size>", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histEtaRecMIP = ROOT.TH1F("histEtaRecMIP", "", 200, -10, 10)
    SetHistStyle2(histEtaRecMIP, "", "#eta", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histPhiRecMIP = ROOT.TH1F("histPhiRecMIP", "", 100, 0, 2*ROOT.TMath.Pi())
    SetHistStyle2(histPhiRecMIP, "", "#phi", "Entries", ROOT.kRed+1, 2, 3001, 0.6)

    histEtaVsMeanClsizePerTrackRecMIP = ROOT.TH2F("histEtaVsMeanClsizePerTrackRecMIP", "", 800, 0, 8, 200, -10, 10)

    effChargeMatch = ROOT.TEfficiency("effChargeMatch", "Charge Match;p_t [GeV];#epsilon", 20, 0, 10)
    
    for treeName in treeNames:
        fIn = TFile.Open(fInName)
        treeReaderInput = TTreeReader(treeName, fIn)
        mChargeRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeRec")
        mEtaRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mEtaRec")
        mPhiRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPhiRec")
        #mClsizeRec = ROOT.TTreeReaderArray(ROOT.int)(treeReaderInput, "mClsizeRec")
        mMeanClsizePerTrackRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mMeanClsizePerTrackRec")
        if isMC:
            mPtGen = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPtGen")
            mChargeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeGen")
            mPdgCodeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPdgCodeGen")
            mIsMIP = ROOT.TTreeReaderValue(ROOT.bool)(treeReaderInput, "mIsMIP")
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
                histEtaRec.Fill(mEtaRec.Get()[0])
                histPhiRec.Fill(mPhiRec.Get()[0])
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
                    if mIsMIP.Get()[0]:
                        histMeanClsizePerTrackRecMIP.Fill(mMeanClsizePerTrackRec.Get()[0])
                        histEtaRecMIP.Fill(mEtaRec.Get()[0])
                        histPhiRecMIP.Fill(mPhiRec.Get()[0])
                        histEtaVsMeanClsizePerTrackRecMIP.Fill(mMeanClsizePerTrackRec.Get()[0], mEtaRec.Get()[0])

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
                    if mIsMIP.Get()[0]:
                        histClsizeRecPerLayerMIP[layer].Fill(nPixel)

                    if pattId <= 0:
                        continue
                    histPattIdRec.Fill(pattId)
                    histPattIdRecPerLayer[layer].Fill(pattId)
                    if mIsMIP.Get()[0]:
                        histPattIdRecPerLayerMIP[layer].Fill(pattId)

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
        if isMC:
            histClsizeRecPerLayerMIP[i].Write()
            histPattIdRecPerLayerMIP[i].Write()
    histClsizeRec.Write()
    histPattIdRec.Write()
    histMeanClsizePerTrackRec.Write()
    histEtaRec.Write()
    histPhiRec.Write()
    if isMC:
        histMeanClsizePerTrackRecMIP.Write()
        histEtaRecMIP.Write()
        histPhiRecMIP.Write()
        histEtaVsMeanClsizePerTrackRecMIP.Write()

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

def read_fast():
    '''
    Function to read the MFTAssesment output for data
    '''
    os.system("root -l process_output_fast.C++")

def compare_data_mc():
    LoadStyle()


    fInData = TFile.Open("../output/MFTAssessmentData_LHC22r_529341_1000_skimmed.root")
    histClsizePerLayerData = []
    histPattIdPerLayerData = []
    histMeanClsizePerTrackRecData = []
    histEtaRecData = []
    histPhiRecData = []
    integralData = []

    histMeanClsizePerTrackRecData = fInData.Get("histMeanClsizePerTrackRec")
    histMeanClsizePerTrackRecData.Rebin(40)
    histMeanClsizePerTrackRecData.Scale(1. / histMeanClsizePerTrackRecData.Integral())
    SetHistStyle3(histMeanClsizePerTrackRecData, "", "<cluster size>", "Entries", ROOT.kBlack, 2, 24, 0.5)
    histEtaRecData = fInData.Get("histEtaRec")
    histEtaRecData.Scale(1. / histEtaRecData.Integral())
    SetHistStyle3(histEtaRecData, "", "#eta", "Entries", ROOT.kBlack, 2, 24, 0.5)
    histPhiRecData = fInData.Get("histPhiRec")
    histPhiRecData.Scale(1. / histPhiRecData.Integral())
    SetHistStyle3(histPhiRecData, "", "#phi", "Entries", ROOT.kBlack, 2, 24, 0.5)
    for iLayer in range(0, 10):
        histClsizePerLayerData.append(fInData.Get("histClsizeRecLayer{}".format(iLayer)))
        integralData.append(histClsizePerLayerData[iLayer].Integral())
        SetHistStyle3(histClsizePerLayerData[iLayer], "MFT Layer {}".format(iLayer), "cluster size", "Entries", ROOT.kBlack, 2, 24, 0.5)
        histPattIdPerLayerData.append(fInData.Get("histPattIdRecLayer{}".format(iLayer)))
        SetHistStyle3(histPattIdPerLayerData[iLayer], "MFT Layer {}".format(iLayer), "pattern ID", "Entries", ROOT.kBlack, 2, 24, 0.5)


    fInMC = TFile.Open("../output/MFTAssessmentMC_skimmed.root")
    histClsizePerLayerMC = []
    histPattIdPerLayerMC = []
    histMeanClsizePerTrackRecMC = []
    histEtaRecMC = []
    histPhiRecMC = []

    histMeanClsizePerTrackRecMC = fInMC.Get("histMeanClsizePerTrackRecMIP")
    histMeanClsizePerTrackRecMC.Rebin(40)
    histMeanClsizePerTrackRecMC.Scale(1. / histMeanClsizePerTrackRecMC.Integral())
    SetHistStyle2(histMeanClsizePerTrackRecMC, "", "<cluster size>", "Entries", ROOT.kRed+1, 1, 3001, 0)
    histEtaRecMC = fInMC.Get("histEtaRec")
    SetHistStyle2(histEtaRecMC, "", "#eta", "Entries", ROOT.kRed+1, 1, 3001, 0)
    histEtaRecMC.Scale(1. / histEtaRecMC.Integral())
    histPhiRecMC = fInMC.Get("histPhiRec")
    histPhiRecMC.Scale(1. / histPhiRecMC.Integral())
    SetHistStyle2(histPhiRecMC, "", "#phi", "Entries", ROOT.kRed+1, 1, 3001, 0)
    for iLayer in range(0, 10):
        histClsizePerLayerMC.append(fInMC.Get("histClsizeRecLayer{}_MIP".format(iLayer)))
        histClsizePerLayerMC[iLayer].Scale(integralData[iLayer] / histClsizePerLayerMC[iLayer].Integral())
        SetHistStyle2(histClsizePerLayerMC[iLayer], "MFT Layer {}".format(iLayer), "cluster size", "Entries", ROOT.kRed+1, 1, 3001, 0)
        histPattIdPerLayerMC.append(fInMC.Get("histPattIdRecLayer{}_MIP".format(iLayer)))
        histPattIdPerLayerMC[iLayer].Scale(integralData[iLayer] / histPattIdPerLayerMC[iLayer].Integral())
        SetHistStyle2(histPattIdPerLayerMC[iLayer], "MFT Layer {}".format(iLayer), "cluster size", "Entries", ROOT.kRed+1, 1, 3001, 0)

    histRatioClsizePerLayer = []
    histRatioPattIdPerLayer = []
    for iLayer in range(0, 10):
        histRatioClsizePerLayer.append(histClsizePerLayerData[iLayer].Clone("histRatioClsizeRecLayer{}".format(iLayer)))
        histRatioClsizePerLayer[iLayer].Divide(histClsizePerLayerMC[iLayer])
        histRatioPattIdPerLayer.append(histPattIdPerLayerData[iLayer].Clone("histRatioPattIdRecLayer{}".format(iLayer)))
        histRatioPattIdPerLayer[iLayer].Divide(histPattIdPerLayerMC[iLayer])

    canvasMeanClsizePerTrackRec = ROOT.TCanvas("canvasMeanClsizePerTrackRec", "canvasMeanClsizePerTrackRec", 800, 600)
    gPad.SetLogy(1)
    histMeanClsizePerTrackRecData.Draw("EP")
    histMeanClsizePerTrackRecMC.Draw("Hsame")
    canvasMeanClsizePerTrackRec.Update()

    canvasEtaRec = ROOT.TCanvas("canvasEtaRec", "canvasEtaRec", 800, 600)
    gPad.SetLogy(1)
    histEtaRecData.Draw("EP")
    histEtaRecMC.Draw("Hsame")
    canvasEtaRec.Update()

    canvasPhiRec = ROOT.TCanvas("canvasPhiRec", "canvasPhiRec", 800, 600)
    gPad.SetLogy(1)
    histPhiRecData.Draw("EP")
    histPhiRecMC.Draw("Hsame")
    canvasPhiRec.Update()


    canvasClsizeRecPerLayer = ROOT.TCanvas("canvasClsizeRecPerLayer", "canvasClsizeRecPerLayer", 3000, 1200)
    canvasClsizeRecPerLayer.Divide(5, 2)
    for iLayer in range(0, 10):
        canvasClsizeRecPerLayer.cd(iLayer+1)
        gPad.SetLogx(1)
        gPad.SetLogy(1)
        histClsizePerLayerData[iLayer].Draw("EP")
        histClsizePerLayerMC[iLayer].Draw("Hsame")
        canvasClsizeRecPerLayer.Update()

    canvasPattIdRecPerLayer = ROOT.TCanvas("canvasPattIdRecPerLayer", "canvasPattIdRecPerLayer", 3000, 1200)
    canvasPattIdRecPerLayer.Divide(5, 2)
    for iLayer in range(0, 10):
        canvasPattIdRecPerLayer.cd(iLayer+1)
        gPad.SetLogx(1)
        gPad.SetLogy(1)
        histPattIdPerLayerData[iLayer].Draw("EP")
        histPattIdPerLayerMC[iLayer].Draw("Hsame")
        canvasPattIdRecPerLayer.Update()



    lineUnity = TLine(0, 1., 1000, 1.)
    lineUnity.SetLineStyle(ROOT.kDashed)
    lineUnity.SetLineWidth(2)
    lineUnity.SetLineColor(ROOT.kGray+1)

    canvasRatioClsizeRecPerLayer = ROOT.TCanvas("canvasRatioClsizeRecPerLayer", "canvasRatioClsizeRecPerLayer", 3000, 1200)
    canvasRatioClsizeRecPerLayer.Divide(5, 2)
    for iLayer in range(0, 10):
        canvasRatioClsizeRecPerLayer.cd(iLayer+1)
        gPad.SetLogx(1)
        gPad.SetLogy(1)
        histRatioClsizePerLayer[iLayer].Draw("EP")
        lineUnity.Draw("same")
        canvasRatioClsizeRecPerLayer.Update()

    canvasRatioPattIdRecPerLayer = ROOT.TCanvas("canvasRatioPattIdRecPerLayer", "canvasRatioPattIdRecPerLayer", 3000, 1200)
    canvasRatioPattIdRecPerLayer.Divide(5, 2)
    for iLayer in range(0, 10):
        canvasRatioPattIdRecPerLayer.cd(iLayer+1)
        gPad.SetLogx(1)
        gPad.SetLogy(1)
        histRatioPattIdPerLayer[iLayer].Draw("EP")
        lineUnity.Draw("same")
        canvasRatioPattIdRecPerLayer.Update()

    input()







    

    exit()

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--read", help="Read the input and produce output histograms", action="store_true")
    parser.add_argument("--read_fast", help="Read fast the input and produce output histograms", action="store_true")
    parser.add_argument("--compare_data_mc", help="Compare Data and MC distributions", action="store_true")
    args = parser.parse_args()

    if args.read:
        read("../output/MFTAssessmentMC.root", "../output/MFTAssessmentMC_skimmed.root")
    if args.read_fast:
        read_fast()
    if args.compare_data_mc:
        compare_data_mc()

main()