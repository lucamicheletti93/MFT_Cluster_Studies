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
    SetHistStyle(mMFTTrackMeanClusterSize)

    fInMC = ROOT.TFile.Open("../output/MFTAssessmentMC.root", "READ")
    mMFTTrueTrackMeanClusterSize = fInMC.Get("mMFTTrueTrackMeanClusterSize")
    mMFTTrueTrackMeanClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackMeanClusterSize.Rebin(10)
    SetHistStyle(mMFTTrueTrackMeanClusterSize)

    mMFTTrueTrackMeanClusterSizeNuclei = fInMC.Get("mMFTTrueTrackMeanClusterSizeNuclei")
    mMFTTrueTrackMeanClusterSizeNuclei.SetLineColor(ROOT.kRed+1)
    mMFTTrueTrackMeanClusterSizeNuclei.SetFillColorAlpha(ROOT.kRed+1, 0.5)
    mMFTTrueTrackMeanClusterSizeNuclei.Rebin(10)

    mMFTTrueTrackMeanClusterSizeNonNuclei = fInMC.Get("mMFTTrueTrackMeanClusterSizeNonNuclei")
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetLineColor(ROOT.kAzure+2)
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetFillColorAlpha(ROOT.kAzure+2, 0.5)
    mMFTTrueTrackMeanClusterSizeNonNuclei.Rebin(10)

    histRatioNuclei = mMFTTrueTrackMeanClusterSizeNuclei.Clone("histRatioNuclei")
    histRatioNuclei.Divide(mMFTTrueTrackMeanClusterSize)

    histRatioNonNuclei = mMFTTrueTrackMeanClusterSizeNonNuclei.Clone("histRatioNonNuclei")
    histRatioNonNuclei.Divide(mMFTTrueTrackMeanClusterSize)

    canvasRatioNuclei = ROOT.TCanvas("canvasRatioNuclei", "canvasRatioNuclei", 800, 600)
    histRatioNuclei.Draw()
    histRatioNonNuclei.Draw("same")
    canvasRatioNuclei.Update()

    legendResults = ROOT.TLegend(0.65, 0.70, 0.85, 0.89, " ", "brNDC")
    SetLegend(legendResults)
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSize, "All", "F")
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSizeNuclei, "Nuclei (Z = 2)", "F")
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSizeNonNuclei, "Non-Nuclei (Z = 1)", "F")

    canvasResultsMC = ROOT.TCanvas("canvasResultsMC", "canvasResultsMC", 800, 600)
    mMFTTrueTrackMeanClusterSize.Draw("H")
    mMFTTrueTrackMeanClusterSizeNuclei.Draw("Hsame")
    mMFTTrueTrackMeanClusterSizeNonNuclei.Draw("Hsame")
    legendResults.Draw("same")
    canvasResultsMC.Update()
    canvasResultsMC.SaveAs("mean_cluster_size_per_track_mc.pdf")

    mMFTTrackMeanClusterSizeNorm = mMFTTrackMeanClusterSize.Clone("mMFTTrackMeanClusterSizeNorm")
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm = mMFTTrueTrackMeanClusterSizeNonNuclei.Clone("mMFTTrueTrackMeanClusterSizeNonNucleiNorm")

    mMFTTrackMeanClusterSizeNorm.Scale(1. / mMFTTrackMeanClusterSizeNorm.GetMaximum())
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Scale(1. / mMFTTrueTrackMeanClusterSizeNonNucleiNorm.GetMaximum())

    legendComparison = ROOT.TLegend(0.65, 0.70, 0.85, 0.89, " ", "brNDC")
    SetLegend(legendComparison)
    legendComparison.AddEntry(mMFTTrueTrackMeanClusterSize, "Data", "F")
    legendComparison.AddEntry(mMFTTrueTrackMeanClusterSizeNonNuclei, "MC (Z = 1)", "F")

    canvasComparison = ROOT.TCanvas("canvasComparison", "canvasComparison", 800, 600)
    mMFTTrackMeanClusterSizeNorm.Draw("EP")
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Draw("Hsame")
    legendComparison.Draw("same")
    canvasComparison.Update()
    canvasComparison.SaveAs("data_mc_comparison.pdf")

    input()


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--analysis", help="plot results", action="store_true")
    args = parser.parse_args()

    if args.analysis:
        analysis()

main()