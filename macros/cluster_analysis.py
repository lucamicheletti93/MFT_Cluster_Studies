import os
import sys
import argparse
import yaml
import random
import ROOT
from ROOT import gROOT, gBenchmark, gPad, gStyle, kTRUE, kFALSE
sys.path.append('../utils')
from plot_library import LoadStyle, SetLegend, SetHistStyle


def analysis():
    '''
    function for producing QC plots
    '''
    LoadStyle()

    ROOT.TGaxis.SetMaxDigits(3)

    fInData = ROOT.TFile.Open("../output/MFTAssessmentData.root", "READ")
    mMFTTrackMeanClusterSize = fInData.Get("mMFTTrackMeanClusterSize")
    mMFTTrackMeanClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrackMeanClusterSize.SetMarkerColor(ROOT.kBlack)
    mMFTTrackMeanClusterSize.SetMarkerStyle(20)
    mMFTTrackMeanClusterSize.Rebin(10)
    mMFTTrackMeanClusterSize.GetXaxis().SetTitle("<Cl size>")
    mMFTTrackMeanClusterSize.GetXaxis().SetRangeUser(0, 8)
    SetHistStyle(mMFTTrackMeanClusterSize)

    mMFTTrackClusterSize = fInData.Get("mMFTTrackClusterSize")
    mMFTTrackClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrackClusterSize.SetMarkerColor(ROOT.kBlack)
    mMFTTrackClusterSize.SetMarkerStyle(20)
    mMFTTrackClusterSize.GetXaxis().SetTitle("Cl size")
    SetHistStyle(mMFTTrackClusterSize)

    fInMC = ROOT.TFile.Open("../output/MFTAssessmentMC.root", "READ")
    mMFTTrueTrackMeanClusterSize = fInMC.Get("mMFTTrueTrackMeanClusterSize")
    mMFTTrueTrackMeanClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackMeanClusterSize.Rebin(10)
    mMFTTrueTrackMeanClusterSize.GetXaxis().SetTitle("<Cl size>")
    mMFTTrueTrackMeanClusterSize.GetXaxis().SetRangeUser(0, 8)
    SetHistStyle(mMFTTrueTrackMeanClusterSize)

    mMFTTrueTrackClusterSize = fInMC.Get("mMFTTrueTrackClusterSize")
    mMFTTrueTrackClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackClusterSize.GetXaxis().SetTitle("Cl size")
    SetHistStyle(mMFTTrueTrackClusterSize)

    mMFTTrueTrackClusterSizeNuclei = fInMC.Get("mMFTTrueTrackClusterSizeNuclei")
    mMFTTrueTrackClusterSizeNuclei.SetLineColor(ROOT.kRed+1)
    mMFTTrueTrackClusterSizeNuclei.SetFillColorAlpha(ROOT.kRed+1, 0.5)
    mMFTTrueTrackClusterSizeNuclei.GetXaxis().SetTitle("Cl size")
    SetHistStyle(mMFTTrueTrackClusterSizeNuclei)

    mMFTTrueTrackClusterSizeNonNuclei = fInMC.Get("mMFTTrueTrackClusterSizeNonNuclei")
    mMFTTrueTrackClusterSizeNonNuclei.SetLineColor(ROOT.kAzure+2)
    mMFTTrueTrackClusterSizeNonNuclei.SetFillColorAlpha(ROOT.kAzure+2, 0.5)
    mMFTTrueTrackClusterSizeNonNuclei.GetXaxis().SetTitle("Cl size")
    SetHistStyle(mMFTTrueTrackClusterSizeNonNuclei)

    mMFTTrueTrackMeanClusterSizeNuclei = fInMC.Get("mMFTTrueTrackMeanClusterSizeNuclei")
    mMFTTrueTrackMeanClusterSizeNuclei.SetLineColor(ROOT.kRed+1)
    mMFTTrueTrackMeanClusterSizeNuclei.SetFillColorAlpha(ROOT.kRed+1, 0.5)
    mMFTTrueTrackMeanClusterSizeNuclei.GetXaxis().SetTitle("<Cl size>")
    mMFTTrueTrackMeanClusterSizeNuclei.GetXaxis().SetRangeUser(0, 8)
    mMFTTrueTrackMeanClusterSizeNuclei.Rebin(10)

    mMFTTrueTrackMeanClusterSizeNonNuclei = fInMC.Get("mMFTTrueTrackMeanClusterSizeNonNuclei")
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetLineColor(ROOT.kAzure+2)
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetFillColorAlpha(ROOT.kAzure+2, 0.5)
    mMFTTrueTrackMeanClusterSizeNonNuclei.GetXaxis().SetTitle("<Cl size>")
    mMFTTrueTrackMeanClusterSizeNonNuclei.GetXaxis().SetRangeUser(0, 8)
    mMFTTrueTrackMeanClusterSizeNonNuclei.Rebin(10)

    histRatioNuclei = mMFTTrueTrackMeanClusterSizeNuclei.Clone("histRatioNuclei")
    histRatioNuclei.Divide(mMFTTrueTrackMeanClusterSize)

    histRatioNonNuclei = mMFTTrueTrackMeanClusterSizeNonNuclei.Clone("histRatioNonNuclei")
    histRatioNonNuclei.Divide(mMFTTrueTrackMeanClusterSize)

    canvasRatioNuclei = ROOT.TCanvas("canvasRatioNuclei", "canvasRatioNuclei", 800, 600)
    histRatioNuclei.Draw()
    histRatioNonNuclei.Draw("same")
    canvasRatioNuclei.Update()

    mMFTTrackMeanClusterSizeNorm = mMFTTrackMeanClusterSize.Clone("mMFTTrackMeanClusterSizeNorm")
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm = mMFTTrueTrackMeanClusterSizeNonNuclei.Clone("mMFTTrueTrackMeanClusterSizeNonNucleiNorm")

    mMFTTrackMeanClusterSizeNorm.Scale(1. / mMFTTrackMeanClusterSizeNorm.GetMaximum())
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Scale(1. / mMFTTrueTrackMeanClusterSizeNonNucleiNorm.GetMaximum())

    mMFTTrackClusterSizeNorm = mMFTTrackClusterSize.Clone("mMFTTrackClusterSizeNorm")
    mMFTTrueTrackClusterSizeNorm = mMFTTrueTrackClusterSize.Clone("mMFTTrueTrackClusterSizeNorm")
    mMFTTrueTrackClusterSizeNucleiNorm = mMFTTrueTrackClusterSizeNuclei.Clone("mMFTTrueTrackClusterSizeNucleiNorm")
    mMFTTrueTrackClusterSizeNonNucleiNorm = mMFTTrueTrackClusterSizeNonNuclei.Clone("mMFTTrueTrackClusterSizeNonNucleiNorm")

    mMFTTrackClusterSizeNorm.Scale(1. / mMFTTrackClusterSize.Integral())
    mMFTTrueTrackClusterSizeNorm.Scale(1. / mMFTTrueTrackClusterSize.Integral())
    mMFTTrueTrackClusterSizeNucleiNorm.Scale(1. / mMFTTrueTrackClusterSizeNuclei.Integral())
    mMFTTrueTrackClusterSizeNonNucleiNorm.Scale(1. / mMFTTrueTrackClusterSizeNonNuclei.Integral())

    legendResults = ROOT.TLegend(0.65, 0.70, 0.85, 0.89, " ", "brNDC")
    SetLegend(legendResults)
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSize, "All", "F")
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSizeNuclei, "Nuclei (Z = 2)", "F")
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSizeNonNuclei, "Non-Nuclei (Z = 1)", "F")

    canvasMeanCLsizeMC = ROOT.TCanvas("canvasMeanCLsizeMC", "canvasMeanCLsizeMC", 800, 600)
    mMFTTrueTrackMeanClusterSize.Draw("H")
    mMFTTrueTrackMeanClusterSizeNuclei.Draw("Hsame")
    mMFTTrueTrackMeanClusterSizeNonNuclei.Draw("Hsame")
    legendResults.Draw("same")
    canvasMeanCLsizeMC.Update()
    canvasMeanCLsizeMC.SaveAs("mean_clsize_per_track_mc.pdf")

    legendComparison = ROOT.TLegend(0.65, 0.70, 0.85, 0.89, " ", "brNDC")
    SetLegend(legendComparison)
    legendComparison.AddEntry(mMFTTrueTrackMeanClusterSize, "Data", "PL")
    legendComparison.AddEntry(mMFTTrueTrackMeanClusterSizeNonNuclei, "MC (Z = 1)", "F")

    canvasMeanCLsizeComp = ROOT.TCanvas("canvasMeanCLsizeComp", "canvasMeanCLsizeComp", 800, 600)
    mMFTTrackMeanClusterSizeNorm.Draw("EP")
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Draw("Hsame")
    legendComparison.Draw("same")
    canvasMeanCLsizeComp.Update()
    canvasMeanCLsizeComp.SaveAs("mean_clsize_data_mc_comparison.pdf")

    legendComparisonCLsize = ROOT.TLegend(0.65, 0.70, 0.85, 0.89, " ", "brNDC")
    SetLegend(legendComparisonCLsize)
    legendComparisonCLsize.AddEntry(mMFTTrackClusterSize, "Data", "LP")
    legendComparisonCLsize.AddEntry(mMFTTrueTrackClusterSize, "MC", "L")


    canvasCLsizeMC = ROOT.TCanvas("canvasCLsizeMC", "canvasCLsizeMC", 800, 600)
    gPad.SetLogx(1)
    gPad.SetLogy(1)
    mMFTTrueTrackClusterSize.Draw("Hsame")
    mMFTTrueTrackClusterSizeNuclei.Draw("Hsame")
    mMFTTrueTrackClusterSizeNonNuclei.Draw("Hsame")
    legendResults.Draw("same")
    canvasCLsizeMC.Update()
    canvasCLsizeMC.SaveAs("clsize_per_track_mc.pdf")

    canvasCLsizeComp = ROOT.TCanvas("canvasCLsizeComp", "canvasCLsizeComp", 800, 600)
    gPad.SetLogx(1)
    gPad.SetLogy(1)
    mMFTTrackClusterSizeNorm.Draw("Hsame")
    mMFTTrueTrackClusterSizeNonNucleiNorm.Draw("Hsame")
    legendComparison.Draw("same")
    canvasCLsizeComp.Update()
    canvasCLsizeComp.SaveAs("clsize_data_mc_comparison.pdf")

    input()


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--analysis", help="plot results", action="store_true")
    args = parser.parse_args()

    if args.analysis:
        analysis()

main()