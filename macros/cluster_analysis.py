import os
import sys
import argparse
import yaml
import random
import ROOT
from ROOT import gROOT, gBenchmark, gPad, gStyle, kTRUE, kFALSE
sys.path.append('../utils')
from plot_library import LoadStyle, SetLegend, SetHistStyle, SetLatex


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
    mMFTTrackMeanClusterSize.GetXaxis().SetTitle("<cluster size>")
    mMFTTrackMeanClusterSize.GetXaxis().SetRangeUser(0, 8)
    SetHistStyle(mMFTTrackMeanClusterSize)

    mMFTTrackClusterSize = fInData.Get("mMFTTrackClusterSize")
    mMFTTrackClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrackClusterSize.SetMarkerColor(ROOT.kBlack)
    mMFTTrackClusterSize.SetMarkerStyle(20)
    mMFTTrackClusterSize.GetXaxis().SetTitle("cluster size")
    SetHistStyle(mMFTTrackClusterSize)

    fInMC = ROOT.TFile.Open("../output/MFTAssessmentMC.root", "READ")
    mMFTTrueTrackMeanClusterSize = fInMC.Get("mMFTTrueTrackMeanClusterSize")
    mMFTTrueTrackMeanClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackMeanClusterSize.SetLineWidth(2)
    mMFTTrueTrackMeanClusterSize.Rebin(10)
    mMFTTrueTrackMeanClusterSize.GetXaxis().SetTitle("<cluster size>")
    mMFTTrueTrackMeanClusterSize.GetXaxis().SetRangeUser(0, 8)
    SetHistStyle(mMFTTrueTrackMeanClusterSize)

    mMFTTrueTrackClusterSize = fInMC.Get("mMFTTrueTrackClusterSize")
    mMFTTrueTrackClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackClusterSize.SetLineWidth(2)
    mMFTTrueTrackClusterSize.GetXaxis().SetTitle("cluster size")
    SetHistStyle(mMFTTrueTrackClusterSize)

    mMFTTrueTrackClusterSizeNuclei = fInMC.Get("mMFTTrueTrackClusterSizeNuclei")
    mMFTTrueTrackClusterSizeNuclei.SetLineColor(ROOT.kRed+1)
    mMFTTrueTrackClusterSizeNuclei.SetLineWidth(2)
    mMFTTrueTrackClusterSizeNuclei.SetFillColorAlpha(ROOT.kRed+1, 0.6)
    mMFTTrueTrackClusterSizeNuclei.GetXaxis().SetTitle("cluster size")
    SetHistStyle(mMFTTrueTrackClusterSizeNuclei)

    mMFTTrueTrackClusterSizeNonNuclei = fInMC.Get("mMFTTrueTrackClusterSizeNonNuclei")
    mMFTTrueTrackClusterSizeNonNuclei.SetLineColor(ROOT.kAzure+2)
    mMFTTrueTrackClusterSizeNonNuclei.SetLineWidth(2)
    mMFTTrueTrackClusterSizeNonNuclei.SetFillColorAlpha(ROOT.kAzure+2, 0.6)
    mMFTTrueTrackClusterSizeNonNuclei.GetXaxis().SetTitle("cluster size")
    SetHistStyle(mMFTTrueTrackClusterSizeNonNuclei)

    mMFTTrueTrackMeanClusterSize3He = fInMC.Get("mMFTTrueTrackMeanClusterSizeNuclei")
    mMFTTrueTrackMeanClusterSizeAnti3He = fInMC.Get("mMFTTrueTrackMeanClusterSizeAntiNuclei")

    mMFTTrueTrackMeanClusterSizeNuclei = mMFTTrueTrackMeanClusterSize3He.Clone("mMFTTrueTrackMeanClusterSizeNuclei")
    mMFTTrueTrackMeanClusterSizeNuclei.Add(mMFTTrueTrackMeanClusterSizeAnti3He)

    #mMFTTrueTrackMeanClusterSizeNuclei = fInMC.Get("mMFTTrueTrackMeanClusterSizeNuclei")
    mMFTTrueTrackMeanClusterSizeNuclei.SetLineColor(ROOT.kRed+1)
    mMFTTrueTrackMeanClusterSizeNuclei.SetLineWidth(2)
    mMFTTrueTrackMeanClusterSizeNuclei.SetFillColorAlpha(ROOT.kRed+1, 0.6)
    mMFTTrueTrackMeanClusterSizeNuclei.GetXaxis().SetTitle("<cluster size>")
    mMFTTrueTrackMeanClusterSizeNuclei.GetXaxis().SetRangeUser(0, 8)
    mMFTTrueTrackMeanClusterSizeNuclei.Rebin(10)
    SetHistStyle(mMFTTrueTrackMeanClusterSizeNuclei)

    mMFTTrueTrackMeanClusterSizeNonNuclei = fInMC.Get("mMFTTrueTrackMeanClusterSizeNonNuclei")
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetLineColor(ROOT.kAzure+2)
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetLineWidth(2)
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetFillColorAlpha(ROOT.kAzure+2, 0.6)
    mMFTTrueTrackMeanClusterSizeNonNuclei.GetXaxis().SetTitle("<cluster size>")
    mMFTTrueTrackMeanClusterSizeNonNuclei.GetXaxis().SetRangeUser(0, 8)
    mMFTTrueTrackMeanClusterSizeNonNuclei.Rebin(10)
    SetHistStyle(mMFTTrueTrackMeanClusterSizeNonNuclei)

    mMFTTrueTrackMeanClusterSizeLF = fInMC.Get("mMFTTrueTrackMeanClusterSizeLF")
    mMFTTrueTrackMeanClusterSizeLF.SetLineColor(ROOT.kAzure+2)
    mMFTTrueTrackMeanClusterSizeLF.SetLineWidth(2)
    mMFTTrueTrackMeanClusterSizeLF.SetFillColorAlpha(ROOT.kAzure+2, 0.6)
    mMFTTrueTrackMeanClusterSizeLF.GetXaxis().SetTitle("<cluster size>")
    mMFTTrueTrackMeanClusterSizeLF.GetXaxis().SetRangeUser(0, 8)
    mMFTTrueTrackMeanClusterSizeLF.Rebin(10)
    SetHistStyle(mMFTTrueTrackMeanClusterSizeLF)

    histRatioNuclei = mMFTTrueTrackMeanClusterSizeNuclei.Clone("histRatioNuclei")
    histRatioNuclei.Divide(mMFTTrueTrackMeanClusterSize)

    histRatioNonNuclei = mMFTTrueTrackMeanClusterSizeNonNuclei.Clone("histRatioNonNuclei")
    histRatioNonNuclei.Divide(mMFTTrueTrackMeanClusterSize)

    canvasRatioNuclei = ROOT.TCanvas("canvasRatioNuclei", "canvasRatioNuclei", 800, 600)
    histRatioNuclei.Draw()
    histRatioNonNuclei.Draw("same")
    canvasRatioNuclei.Update()

    mMFTTrueTrackMeanClusterSizeSum = mMFTTrueTrackMeanClusterSizeLF.Clone("mMFTTrueTrackMeanClusterSizeSum")
    mMFTTrueTrackMeanClusterSizeSum.Add(mMFTTrueTrackMeanClusterSizeNuclei)
    mMFTTrueTrackMeanClusterSizeSum.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackMeanClusterSizeSum.SetFillColor(0)
    mMFTTrueTrackMeanClusterSizeSum.SetLineWidth(2)
    mMFTTrueTrackMeanClusterSizeSum.GetXaxis().SetTitle("<cluster size>")
    mMFTTrueTrackMeanClusterSizeSum.GetXaxis().SetRangeUser(0, 8)
    SetHistStyle(mMFTTrueTrackMeanClusterSizeSum)

    mMFTTrackMeanClusterSizeNorm = mMFTTrackMeanClusterSize.Clone("mMFTTrackMeanClusterSizeNorm")
    mMFTTrueTrackMeanClusterSizeNorm = mMFTTrueTrackMeanClusterSize.Clone("mMFTTrueTrackMeanClusterSizeNorm")
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm = mMFTTrueTrackMeanClusterSizeNonNuclei.Clone("mMFTTrueTrackMeanClusterSizeNonNucleiNorm")
    mMFTTrueTrackMeanClusterSizeNucleiNorm = mMFTTrueTrackMeanClusterSizeNuclei.Clone("mMFTTrueTrackMeanClusterSizeNucleiNorm")
    mMFTTrueTrackMeanClusterSizeLFNorm = mMFTTrueTrackMeanClusterSizeLF.Clone("mMFTTrueTrackMeanClusterSizeLF")

    mMFTTrackMeanClusterSizeNorm.Scale(1. / mMFTTrackMeanClusterSizeNorm.GetMaximum())
    mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Scale(1. / mMFTTrueTrackMeanClusterSizeNonNucleiNorm.GetMaximum())
    mMFTTrueTrackMeanClusterSizeNucleiNorm.Scale(1. / mMFTTrueTrackMeanClusterSizeNucleiNorm.GetMaximum())
    mMFTTrueTrackMeanClusterSizeLFNorm.Scale(1. / mMFTTrueTrackMeanClusterSizeLFNorm.GetMaximum())

    mMFTTrackClusterSizeNorm = mMFTTrackClusterSize.Clone("mMFTTrackClusterSizeNorm")
    mMFTTrueTrackClusterSizeNorm = mMFTTrueTrackClusterSize.Clone("mMFTTrueTrackClusterSizeNorm")
    mMFTTrueTrackClusterSizeNucleiNorm = mMFTTrueTrackClusterSizeNuclei.Clone("mMFTTrueTrackClusterSizeNucleiNorm")
    mMFTTrueTrackClusterSizeNonNucleiNorm = mMFTTrueTrackClusterSizeNonNuclei.Clone("mMFTTrueTrackClusterSizeNonNucleiNorm")

    mMFTTrackClusterSizeNorm.Scale(1. / mMFTTrackClusterSize.Integral())
    mMFTTrueTrackClusterSizeNorm.Scale(1. / mMFTTrueTrackClusterSize.Integral())
    mMFTTrueTrackClusterSizeNucleiNorm.Scale(1. / mMFTTrueTrackClusterSizeNuclei.Integral())
    mMFTTrueTrackClusterSizeNonNucleiNorm.Scale(1. / mMFTTrueTrackClusterSizeNonNuclei.Integral())

    latexTitleResults = ROOT.TLatex()
    latexTitleResults.SetTextSize(0.050)
    latexTitleResults.SetNDC()
    latexTitleResults.SetTextFont(42)

    legendResults = ROOT.TLegend(0.70, 0.50, 0.89, 0.73, " ", "brNDC")
    SetLegend(legendResults)
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSize, "All", "F")
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSizeNuclei, "^{3}He / ^{3}#bar{He}", "F")
    legendResults.AddEntry(mMFTTrueTrackMeanClusterSizeNonNuclei, "Non-Nuclei", "F")

    # Mean CLsize per track MC 
    canvasMeanCLsizeMC = ROOT.TCanvas("canvasMeanCLsizeMC", "canvasMeanCLsizeMC", 800, 600)
    mMFTTrueTrackMeanClusterSize.Draw("HE")
    mMFTTrueTrackMeanClusterSizeNonNuclei.Draw("HEsame")
    mMFTTrueTrackMeanClusterSizeNuclei.Draw("HEsame")
    legendResults.Draw("same")
    latexTitleResults.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitleResults.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    latexTitleResults.DrawLatex(0.6, 0.73, "Run3 simulation")
    canvasMeanCLsizeMC.Update()
    canvasMeanCLsizeMC.SaveAs("mean_clsize_per_track_mc.pdf")

    # Mean CLsize per track MC 
    legendResults_pi_k_p = ROOT.TLegend(0.70, 0.50, 0.89, 0.73, " ", "brNDC")
    SetLegend(legendResults_pi_k_p)
    legendResults_pi_k_p.AddEntry(mMFTTrueTrackMeanClusterSize, "All", "F")
    legendResults_pi_k_p.AddEntry(mMFTTrueTrackMeanClusterSizeNuclei, "^{3}He / ^{3}#bar{He}", "F")
    legendResults_pi_k_p.AddEntry(mMFTTrueTrackMeanClusterSizeLF, "#pi, K, p", "F")

    canvasMeanCLsizeMC_pi_k_p = ROOT.TCanvas("canvasMeanCLsizeMC_pi_k_p", "canvasMeanCLsizeMC_pi_k_p", 800, 600)
    #mMFTTrueTrackMeanClusterSize.Draw("HE")
    mMFTTrueTrackMeanClusterSizeSum.Draw("HE")
    mMFTTrueTrackMeanClusterSizeLF.Draw("HEsame")
    mMFTTrueTrackMeanClusterSizeNuclei.Draw("HEsame")
    legendResults_pi_k_p.Draw("same")
    latexTitleResults.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitleResults.DrawLatex(0.6, 0.79, "-3.6 < #eta < -2.4")
    latexTitleResults.DrawLatex(0.6, 0.73, "Run3 simulation")
    canvasMeanCLsizeMC_pi_k_p.Update()
    canvasMeanCLsizeMC_pi_k_p.SaveAs("mean_clsize_per_track_mc_pi_k_p.pdf")

    # Mean CLsize per track MC Normalized
    legendResultsMCNorm = ROOT.TLegend(0.70, 0.50, 0.89, 0.73, " ", "brNDC")
    SetLegend(legendResultsMCNorm)
    legendResultsMCNorm.AddEntry(mMFTTrueTrackMeanClusterSizeNucleiNorm, "^{3}He / ^{3}#bar{He}", "F")
    legendResultsMCNorm.AddEntry(mMFTTrueTrackMeanClusterSizeLFNorm, "#pi, K, p", "F")

    canvasMeanCLsizeMCNorm = ROOT.TCanvas("canvasMeanCLsizeMCNorm", "canvasMeanCLsizeMCNorm", 800, 600)
    mMFTTrueTrackMeanClusterSizeLFNorm.GetXaxis().SetRangeUser(0, 12)
    mMFTTrueTrackMeanClusterSizeLFNorm.GetYaxis().SetRangeUser(0, 1.1)
    #mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Draw("H")
    mMFTTrueTrackMeanClusterSizeLFNorm.Draw("Hsame")
    mMFTTrueTrackMeanClusterSizeNucleiNorm.Draw("Hsame")
    legendResultsMCNorm.Draw("same")
    latexTitleResults.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitleResults.DrawLatex(0.6, 0.80, "-3.6 < #eta < -2.4")
    latexTitleResults.DrawLatex(0.6, 0.75, "Run3 simulation")
    canvasMeanCLsizeMCNorm.Update()
    canvasMeanCLsizeMCNorm.SaveAs("mean_clsize_per_track_mc_norm.pdf")

    #mMFTTrueTrackMeanClusterSizeNucleiNorm.Fit("gaus")
    #mMFTTrueTrackMeanClusterSizeNonNucleiNorm.Fit("gaus")

    legendComparison = ROOT.TLegend(0.65, 0.70, 0.85, 0.89, " ", "brNDC")
    SetLegend(legendComparison)
    legendComparison.AddEntry(mMFTTrueTrackMeanClusterSize, "Data", "PL")
    legendComparison.AddEntry(mMFTTrueTrackMeanClusterSizeNonNuclei, "MC (Z = 1)", "F")

    # Data - MC Mean CLsize comparison
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

    # MC - MC Mean CLsize comparison
    canvasCLsizeMC = ROOT.TCanvas("canvasCLsizeMC", "canvasCLsizeMC", 800, 600)
    gPad.SetLogx(1)
    gPad.SetLogy(1)
    mMFTTrueTrackClusterSize.Draw("Hsame")
    mMFTTrueTrackClusterSizeNuclei.Draw("Hsame")
    mMFTTrueTrackClusterSizeNonNuclei.Draw("Hsame")
    legendResults.Draw("same")
    canvasCLsizeMC.Update()
    canvasCLsizeMC.SaveAs("clsize_per_track_mc.pdf")

    # Data - MC CLsize comparison
    canvasCLsizeComp = ROOT.TCanvas("canvasCLsizeComp", "canvasCLsizeComp", 800, 600)
    gPad.SetLogx(1)
    gPad.SetLogy(1)
    mMFTTrackClusterSizeNorm.Draw("Hsame")
    mMFTTrueTrackClusterSizeNonNucleiNorm.Draw("Hsame")
    legendComparison.Draw("same")
    canvasCLsizeComp.Update()
    canvasCLsizeComp.SaveAs("clsize_data_mc_comparison.pdf")
    

    # Charge efficiency selection
    treeReaderInput = ROOT.TTreeReader("treeTrackClusterSize", fInMC)
    mChargeRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeRec")
    mPtRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPtRec")
    mEtaRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mEtaRec")
    mPhiRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPhiRec")
    #mClsizeRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mClsizeRec")
    #mLayerRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mLayerRec")
    #mPattIdRec = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPattIdRec")
    mMeanClsizePerTrackRec = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mMeanClsizePerTrackRec")
    mChargeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mChargeGen")
    mPtGen = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPtGen")
    mEtaGen = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mEtaGen")
    mPhiGen = ROOT.TTreeReaderValue(ROOT.double)(treeReaderInput, "mPhiGen")
    mPdgCodeGen = ROOT.TTreeReaderValue(ROOT.int)(treeReaderInput, "mPdgCodeGen")

    effChargeMatch = ROOT.TEfficiency("effChargeMatch", "Charge Match;p_t [GeV];#epsilon", 20, 0, 10)
    effChargeMatchNuclei = ROOT.TEfficiency("effChargeMatchNuclei", "Charge Match;p_t [GeV];#epsilon", 20, 0, 10)

    print("Processing {}...".format(fInMC.GetName()))
    while treeReaderInput.Next():
        deltaCharge = mChargeRec.Get()[0] - mChargeGen.Get()[0]
        effChargeMatch.Fill(not deltaCharge, mPtGen.Get()[0])

    histTotalChargeMatch = effChargeMatch.GetCopyTotalHisto()
    histPassedChargeMatch = effChargeMatch.GetCopyPassedHisto()
    histEffChargeMatch = ROOT.TH1D("histEffChargeMatch", "", 20, 0, 10)
    histEffChargeMatch.Divide(histPassedChargeMatch, histTotalChargeMatch, 1, 1, "B")
    #histEffChargeMatch.SetMarkerStyle(20)
    #histEffChargeMatch.SetMarkerSize(0.8)
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

    latexTitle = ROOT.TLatex()
    latexTitle.SetTextSize(0.050)
    latexTitle.SetNDC()
    latexTitle.SetTextFont(42)

    canvasChargeMatch = ROOT.TCanvas("canvasChargeMatch", "", 800, 600)
    #gPad.SetLogy()
    histEffChargeMatch.Draw("E")
    lineUnity.Draw("same")
    latexTitle.DrawLatex(0.6, 0.85, "pp #sqrt{s} = 13 TeV")
    latexTitle.DrawLatex(0.6, 0.80, "Run3 simulation")
    canvasChargeMatch.Update()
    canvasChargeMatch.SaveAs("charge_match_eff.pdf")


    input()


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument("--analysis", help="plot results", action="store_true")
    args = parser.parse_args()

    if args.analysis:
        analysis()

main()