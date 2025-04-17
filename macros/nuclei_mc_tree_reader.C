#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TKey.h"
#include "THashList.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include <ROOT/RDataFrame.hxx>

#endif

void LoadStyle();

inline void SetHist(auto *hist, Color_t mkrCol, int mkrSty, double mkrSz, Color_t lnCol, int lnWidth, int fillSty, double alpha = 1) {
    hist -> SetMarkerColorAlpha(mkrCol, alpha);
    hist -> SetMarkerStyle(mkrSty);
    hist -> SetMarkerSize(mkrSz);
    hist -> SetLineColorAlpha(lnCol, alpha);
    hist -> SetLineWidth(lnWidth);
    hist -> SetFillStyle(fillSty);
}

void nuclei_mc_tree_reader() {
    LoadStyle();

    float fPosX, fPosY, fPosZ, fPt, fEta, fPhi, fFwdDcaX, fFwdDcaY, fChi2MatchMCHMID, fChi2MatchMCHMFT = -99999;
    uint64_t fMftClusterSizesAndTrackFlags;
    uint16_t fMcDecision;
    UChar_t fTrackType;
    int fSign;

    int mTrackType = -999;
    float mPt = -999;
    float mEta = -999;
    float mPhi = -999;
    int mSign = -999;
    float mMch2MCHMID = -999;
    float mMch2MFTMCH = -999;
    float mFwdDcaX = -999;
    float mFwdDcaY = -999;
    float mFwdDcaZ = -999;
    int mClsizeLayer0 = -999;
    int mClsizeLayer1 = -999;
    int mClsizeLayer2 = -999;
    int mClsizeLayer3 = -999;
    int mClsizeLayer4 = -999;
    int mClsizeLayer5 = -999;
    int mClsizeLayer6 = -999;
    int mClsizeLayer7 = -999;
    int mClsizeLayer8 = -999;
    int mClsizeLayer9 = -999;
    double mMeanClsizePerTrack = -999;

    TH1D *histMeanClSize = new TH1D("histMeanClSize", "; <Cluster Size>", 60, 0, 30); SetHist(histMeanClSize, kBlack, 20, 0.8, kBlack, 1, 0, 1);
    TH1D *histEta = new TH1D("histEta", "; #eta", 1000, -5, 5); SetHist(histEta, kBlack, 20, 0.8, kBlack, 1, 0, 1);
    TH1D *histFwdDcaX = new TH1D("histFwdDcaX", "Fwd. DCA_{x} ; DCA_{x} (cm)", 200, -10, 10); SetHist(histFwdDcaX, kBlack, 20, 0.8, kBlack, 1, 0, 1);
    TH1D *histFwdDcaY = new TH1D("histFwdDcaY", "Fwd. DCA_{y} ; DCA_{y} (cm)", 200, -10, 10); SetHist(histFwdDcaY, kBlack, 20, 0.8, kBlack, 1, 0, 1);
    TH1D *histFwdDcaZ = new TH1D("histFwdDcaZ", "Fwd. DCA_{z} ; DCA_{z} (cm)", 200, -1, 19); SetHist(histFwdDcaZ, kBlack, 20, 0.8, kBlack, 1, 0, 1);

    const int nMcSignals = 7;
    string mcSignals[] = {"He^{3}_{primary}", "He^{3}_{from transport}", "e", "#mu", "#pi", "K", "p"};
    Color_t mcSignalsColor[] = {kRed+1, kOrange+7, kBlue+1, kGreen+2, kAzure+4, kViolet+7, kMagenta};
    TH1D *histMcSignals = new TH1D("histMcSignals", "; MC signals", nMcSignals+1, 0, nMcSignals+1);

    TH1D *histMeanClSizeMcSignals[nMcSignals];
    TH1D *histEtaMcSignals[nMcSignals];
    TH1D *histFwdDcaXMcSignals[nMcSignals];
    TH1D *histFwdDcaYMcSignals[nMcSignals];
    TH1D *histFwdDcaZMcSignals[nMcSignals];
    TH1D *histNormFwdDcaXMcSignals[nMcSignals];
    TH1D *histNormFwdDcaYMcSignals[nMcSignals];
    TH1D *histNormFwdDcaZMcSignals[nMcSignals];

    histMcSignals -> GetXaxis() -> SetBinLabel(1, "All");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histMcSignals -> GetXaxis() -> SetBinLabel(2+iMcSignal, mcSignals[iMcSignal].c_str());
        histMeanClSizeMcSignals[iMcSignal] = new TH1D(Form("histMeanClSize_%s", mcSignals[iMcSignal].c_str()), "; <Cluster Size>", 60, 0, 30);
        histEtaMcSignals[iMcSignal] = new TH1D(Form("histEta_%s", mcSignals[iMcSignal].c_str()), "; #eta", 1000, -5, 5);
        histFwdDcaXMcSignals[iMcSignal] = new TH1D(Form("histFwdDcaXMcSignals_%s", mcSignals[iMcSignal].c_str()), "; DCA_{x}", 200, -10, 10);
        histFwdDcaYMcSignals[iMcSignal] = new TH1D(Form("histFwdDcaYMcSignals_%s", mcSignals[iMcSignal].c_str()), "; DCA_{y}", 200, -10, 10);
        histFwdDcaZMcSignals[iMcSignal] = new TH1D(Form("histFwdDcaZMcSignals_%s", mcSignals[iMcSignal].c_str()), "; DCA_{z}", 200, -1, 19);

        SetHist(histMeanClSizeMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histEtaMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histFwdDcaXMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histFwdDcaYMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histFwdDcaZMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histMeanClSizeMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histFwdDcaXMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histFwdDcaYMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);
        SetHist(histFwdDcaZMcSignals[iMcSignal], mcSignalsColor[iMcSignal], 20, 0, mcSignalsColor[iMcSignal], 1, 0, 0.7);

        if (mcSignals[iMcSignal] == "He^{3}_{primary}") {
            histMeanClSizeMcSignals[iMcSignal] -> SetFillStyle(3352); histMeanClSizeMcSignals[iMcSignal] -> SetFillColor(kRed+1);
            histEtaMcSignals[iMcSignal] -> SetFillStyle(3352); histEtaMcSignals[iMcSignal] -> SetFillColor(kRed+1);
            histFwdDcaXMcSignals[iMcSignal] -> SetFillStyle(3352); histFwdDcaXMcSignals[iMcSignal] -> SetFillColor(kRed+1);
            histFwdDcaYMcSignals[iMcSignal] -> SetFillStyle(3352); histFwdDcaYMcSignals[iMcSignal] -> SetFillColor(kRed+1);
            histFwdDcaZMcSignals[iMcSignal] -> SetFillStyle(3352); histFwdDcaZMcSignals[iMcSignal] -> SetFillColor(kRed+1);
        }
    }

    TFile *fIn = new TFile("/Users/lucamicheletti/alice/local_train_test_mc_nuclei/reducedAO2D.root", "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

        TTree *treeIn = (TTree*) fIn -> Get(Form("%s/O2rtfwdpidall", dirName.Data()));
        treeIn -> SetBranchAddress("fPosX", &fPosX);
        treeIn -> SetBranchAddress("fPosY", &fPosY);
        treeIn -> SetBranchAddress("fPosZ", &fPosZ);
        treeIn -> SetBranchAddress("fPt", &fPt);
        treeIn -> SetBranchAddress("fEta", &fEta);
        treeIn -> SetBranchAddress("fPhi", &fPhi);
        treeIn -> SetBranchAddress("fSign", &fSign);
        treeIn -> SetBranchAddress("fFwdDcaX", &fFwdDcaX);
        treeIn -> SetBranchAddress("fFwdDcaY", &fFwdDcaY);
        treeIn -> SetBranchAddress("fChi2MatchMCHMID", &fChi2MatchMCHMID);
        treeIn -> SetBranchAddress("fChi2MatchMCHMFT", &fChi2MatchMCHMFT);
        treeIn -> SetBranchAddress("fMftClusterSizesAndTrackFlags", &fMftClusterSizesAndTrackFlags);
        treeIn -> SetBranchAddress("fTrackType", &fTrackType);
        treeIn -> SetBranchAddress("fMcDecision", &fMcDecision);

        for (int iEntry = 0;iEntry < treeIn -> GetEntries();iEntry++) {
            treeIn -> GetEntry(iEntry);
            if (fEta > -2.5 || fEta < -3.6) continue;
            double meanClSize = 0;
            int nClusters = 0;
            for (int i = 0; i < 10; ++i) {
                double size = (fMftClusterSizesAndTrackFlags >> (i * 6)) & 0x3fULL;
                if (size > 62) {
                    continue;
                }
                if (size > 0) {
                    meanClSize += size;
                    nClusters += 1;
                } else {
                    size = -999;
                }

                switch(i) {
                    case 0:
                        mClsizeLayer0 = size;
                        break;
                    case 1:
                        mClsizeLayer1 = size;
                        break;
                    case 2:
                        mClsizeLayer2 = size;
                        break;
                    case 3:
                        mClsizeLayer3 = size;
                        break;
                    case 4:
                        mClsizeLayer4 = size;
                        break;
                    case 5:
                        mClsizeLayer5 = size;
                        break;
                    case 6:
                        mClsizeLayer6 = size;
                        break;
                    case 7:
                        mClsizeLayer7 = size;
                        break;
                    case 8:
                        mClsizeLayer8 = size;
                        break;
                    case 9:
                        mClsizeLayer9 = size;
                        break;
                }
            }
            meanClSize /= nClusters;

            mTrackType = fTrackType;
            mPt = fPt;
            mEta = fEta;
            mPhi = fPhi;
            mSign = fSign;
            mMch2MCHMID = fChi2MatchMCHMID;
            mMch2MFTMCH = fChi2MatchMCHMFT;
            mFwdDcaX = fFwdDcaX;
            mFwdDcaY = fFwdDcaY;
            mFwdDcaZ = - TMath::SinH(mEta) * TMath::Sqrt(mFwdDcaX*mFwdDcaX + mFwdDcaY*mFwdDcaY);
            histMeanClSize -> Fill(meanClSize);
            histEta -> Fill(mEta);
            histFwdDcaX -> Fill(mFwdDcaX);
            histFwdDcaY -> Fill(mFwdDcaY);
            histFwdDcaZ -> Fill(mFwdDcaZ);
            histMcSignals -> AddBinContent(1, 1);

            for (int i = 0; i < nMcSignals; ++i) {
                bool flag = (fMcDecision >> i) & 1;
                if (flag) {
                    histMeanClSizeMcSignals[i] -> Fill(meanClSize);
                    histEtaMcSignals[i] -> Fill(mEta);
                    histMcSignals -> AddBinContent(2+i, 1);
                    histFwdDcaXMcSignals[i] -> Fill(mFwdDcaX);
                    histFwdDcaYMcSignals[i] -> Fill(mFwdDcaY);
                    histFwdDcaZMcSignals[i] -> Fill(mFwdDcaZ);
                }
            }
        }
    }
    fIn -> Close();

    TCanvas *canvasMcSignals = new TCanvas("canvasMcSignals", "");
    gPad -> SetLogy(true);
    histMcSignals -> Draw("H");

    TCanvas *canvasMftCluseterSizePerTrack = new TCanvas("canvasMftCluseterSizePerTrack", "");
    gPad -> SetLogy(true);
    histMeanClSize -> Draw("EP");

    TLegend *legend = new TLegend(0.70, 0.20, 0.85, 0.90, " ", "brNDC");
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
    legend -> AddEntry(histMeanClSize, "All", "L");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histMeanClSizeMcSignals[iMcSignal] -> Draw("H SAME");
        legend -> AddEntry(histMeanClSizeMcSignals[iMcSignal], mcSignals[iMcSignal].c_str(), "L");
    }
    legend -> Draw();

    TCanvas *canvasEta = new TCanvas("canvasEta", "");
    gPad -> SetLogy(true);
    histEta -> GetYaxis() -> SetRangeUser(1, 1e7);
    histEta -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histEtaMcSignals[iMcSignal] -> Draw("H SAME");
    }

    TCanvas *canvasFwdDca = new TCanvas("canvasFwdDca", "", 1200, 800);
    canvasFwdDca -> Divide(3, 2);
    canvasFwdDca -> cd(1);
    gPad -> SetLogy(true);
    histFwdDcaX -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histFwdDcaXMcSignals[iMcSignal] -> Draw("H SAME");
    }

    canvasFwdDca -> cd(2);
    gPad -> SetLogy(true);
    histFwdDcaY -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histFwdDcaYMcSignals[iMcSignal] -> Draw("H SAME");
    }

    canvasFwdDca -> cd(3);
    gPad -> SetLogy(true);
    histFwdDcaZ -> GetYaxis() -> SetRangeUser(1, 1e7);
    histFwdDcaZ -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histFwdDcaZMcSignals[iMcSignal] -> Draw("H SAME");
    }

    canvasFwdDca -> cd(4);
    gPad -> SetLogy(true);
    TH1D *histNormFwdDcaX = (TH1D*) histFwdDcaX -> Clone("histNormFwdDcaX");
    histNormFwdDcaX -> SetTitle("Normalized Fwd. DCA_{x}");
    histNormFwdDcaX -> Scale(1. / histNormFwdDcaX -> Integral());
    histNormFwdDcaX -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histNormFwdDcaXMcSignals[iMcSignal] = (TH1D*) histFwdDcaXMcSignals[iMcSignal] -> Clone(Form("histNormFwdDcaX_%s", mcSignals[iMcSignal].c_str()));
        histNormFwdDcaXMcSignals[iMcSignal] -> Scale(1. / histNormFwdDcaXMcSignals[iMcSignal] -> Integral());
        histNormFwdDcaXMcSignals[iMcSignal] -> Draw("H SAME");
    }

    canvasFwdDca -> cd(5);
    gPad -> SetLogy(true);
    TH1D *histNormFwdDcaY = (TH1D*) histFwdDcaY -> Clone("histNormFwdDcaY");
    histNormFwdDcaY -> SetTitle("Normalized Fwd. DCA_{y}");
    histNormFwdDcaY -> Scale(1. / histNormFwdDcaY -> Integral());
    histNormFwdDcaY -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histNormFwdDcaYMcSignals[iMcSignal] = (TH1D*) histFwdDcaYMcSignals[iMcSignal] -> Clone(Form("histNormFwdDcaY_%s", mcSignals[iMcSignal].c_str()));
        histNormFwdDcaYMcSignals[iMcSignal] -> Scale(1. / histNormFwdDcaYMcSignals[iMcSignal] -> Integral());
        histNormFwdDcaYMcSignals[iMcSignal] -> Draw("H SAME");
    }

    canvasFwdDca -> cd(6);
    gPad -> SetLogy(true);
    TH1D *histNormFwdDcaZ = (TH1D*) histFwdDcaZ -> Clone("histNormFwdDcaZ");
    histNormFwdDcaZ -> SetTitle("Normalized Fwd. DCA_{z}");
    histNormFwdDcaZ -> Scale(1. / histNormFwdDcaZ -> Integral());
    histNormFwdDcaZ -> GetYaxis() -> SetRangeUser(1e-4, 1e-1);
    histNormFwdDcaZ -> Draw("EP");
    for (int iMcSignal = 0;iMcSignal < nMcSignals;iMcSignal++) {
        histNormFwdDcaZMcSignals[iMcSignal] = (TH1D*) histFwdDcaZMcSignals[iMcSignal] -> Clone(Form("histNormFwdDcaZ_%s", mcSignals[iMcSignal].c_str()));
        histNormFwdDcaZMcSignals[iMcSignal] -> Scale(1. / histNormFwdDcaZMcSignals[iMcSignal] -> Integral());
        histNormFwdDcaZMcSignals[iMcSignal] -> Draw("H SAME");
    }

    canvasMcSignals -> SaveAs("qa_figures/McSignalsSummay.pdf");
    canvasMftCluseterSizePerTrack -> SaveAs("qa_figures/MftCluseterSizePerTrackSummay.pdf");
    canvasEta -> SaveAs("qa_figures/EtaSummay.pdf");
    canvasFwdDca -> SaveAs("qa_figures/FwdDcaSummay.pdf");

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.05,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.05,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetPadBottomMargin(0.15);
    gStyle -> SetPadTopMargin(0.05);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetEndErrorSize(0.0);
    gStyle -> SetTitleSize(0.05,"X");
    gStyle -> SetTitleSize(0.045,"Y");
    gStyle -> SetLabelSize(0.045,"X");
    gStyle -> SetLabelSize(0.045,"Y");
    gStyle -> SetTitleOffset(1.2,"X");
    gStyle -> SetTitleOffset(1.35,"Y");
}