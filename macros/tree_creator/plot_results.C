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
void SetLegend(TLegend *);

void plot_results() {
    LoadStyle();

    const int nDatasets = 3;
    string datasets[] = {"LHC23k6c_pass1", "LHC23zs_pass2", "LHC23zzh_pass1_small"};
    int colors[] = {633, 864, 417};

    TH1F *histMeanClsizePerTrack[nDatasets];
    TH1F *histClsizePerLayer[nDatasets][10];
    
    for (int iDataset = 0;iDataset < nDatasets;iDataset++) {
        histMeanClsizePerTrack[iDataset] = new TH1F(Form("histMeanClsizePerTrack_%i", iDataset), "", 48, 0, 15);
        histMeanClsizePerTrack[iDataset] -> SetLineColor(colors[iDataset]);
        histMeanClsizePerTrack[iDataset] -> SetLineWidth(2);
        histMeanClsizePerTrack[iDataset] -> GetXaxis() -> SetRangeUser(0, 12);
        histMeanClsizePerTrack[iDataset] -> GetXaxis() -> SetTitle("<Cluser size>");
        histMeanClsizePerTrack[iDataset] -> GetYaxis() -> SetRangeUser(1e-5, 1);
        histMeanClsizePerTrack[iDataset] -> GetYaxis() -> SetTitle("Counts");

        for (int i = 0;i < 10;i++) {
            histClsizePerLayer[iDataset][i] = new TH1F(Form("histClsizePerLayer_%i_%i", iDataset, i), "", 100, -0.5, 99.5);
            histClsizePerLayer[iDataset][i] -> SetLineColor(colors[iDataset]);
        }
    }

    int mTrackType = -999;
    float mPt = -999;
    float mEta = -999;
    float mPhi = -999;
    int mSign = -999;
    float mMch2MCHMID = -999;
    float mMch2MFTMCH = -999;
    float mFwdDcaX = -999;
    float mFwdDcaY = -999;
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

    for (int iDataset = 0;iDataset < nDatasets;iDataset++) {
        Printf("[info] Processing /Users/lucamicheletti/cernbox/MTF_Cluster_Data/%s/treeMftPidTraining.root", datasets[iDataset].c_str());
        string pathToFile = Form("/Users/lucamicheletti/cernbox/MTF_Cluster_Data/%s/treeMftPidTraining.root", datasets[iDataset].c_str());
        TFile *fIn = new TFile(pathToFile.c_str(), "READ");

        TTree *treeIn = (TTree*) fIn -> Get("treeTrackClusterSizeData");
        treeIn -> SetBranchAddress("mTrackType", &mTrackType);
        treeIn -> SetBranchAddress("mPt", &mPt);
        treeIn -> SetBranchAddress("mEta", &mEta);
        treeIn -> SetBranchAddress("mPhi", &mPhi);
        treeIn -> SetBranchAddress("mSign", &mSign);
        treeIn -> SetBranchAddress("mMch2MCHMID", &mMch2MCHMID);
        treeIn -> SetBranchAddress("mMch2MFTMCH", &mMch2MFTMCH);
        treeIn -> SetBranchAddress("mFwdDcaX", &mFwdDcaX);
        treeIn -> SetBranchAddress("mFwdDcaY", &mFwdDcaY);
        treeIn -> SetBranchAddress("mClsizeLayer0", &mClsizeLayer0);
        treeIn -> SetBranchAddress("mClsizeLayer1", &mClsizeLayer1);
        treeIn -> SetBranchAddress("mClsizeLayer2", &mClsizeLayer2);
        treeIn -> SetBranchAddress("mClsizeLayer3", &mClsizeLayer3);
        treeIn -> SetBranchAddress("mClsizeLayer4", &mClsizeLayer4);
        treeIn -> SetBranchAddress("mClsizeLayer5", &mClsizeLayer5);
        treeIn -> SetBranchAddress("mClsizeLayer6", &mClsizeLayer6);
        treeIn -> SetBranchAddress("mClsizeLayer7", &mClsizeLayer7);
        treeIn -> SetBranchAddress("mClsizeLayer8", &mClsizeLayer8);
        treeIn -> SetBranchAddress("mClsizeLayer9", &mClsizeLayer9);
        treeIn -> SetBranchAddress("mMeanClsizePerTrack", &mMeanClsizePerTrack);

        for (int iEntry = 0;iEntry < treeIn -> GetEntries();iEntry++) {
            treeIn -> GetEntry(iEntry);
            histMeanClsizePerTrack[iDataset] -> Fill(mMeanClsizePerTrack);

            histClsizePerLayer[iDataset][0] -> Fill(mClsizeLayer0);
            histClsizePerLayer[iDataset][1] -> Fill(mClsizeLayer1);
            histClsizePerLayer[iDataset][2] -> Fill(mClsizeLayer2);
            histClsizePerLayer[iDataset][3] -> Fill(mClsizeLayer3);
            histClsizePerLayer[iDataset][4] -> Fill(mClsizeLayer4);
            histClsizePerLayer[iDataset][5] -> Fill(mClsizeLayer5);
            histClsizePerLayer[iDataset][6] -> Fill(mClsizeLayer6);
            histClsizePerLayer[iDataset][7] -> Fill(mClsizeLayer7);
            histClsizePerLayer[iDataset][8] -> Fill(mClsizeLayer8);
            histClsizePerLayer[iDataset][9] -> Fill(mClsizeLayer9);
        }
    }

    TCanvas *canvasMeanClsizePerTrack = new TCanvas("canvasMeanClsizePerTrack", "", 800, 600);
    canvasMeanClsizePerTrack -> SetFillColor(0);
    canvasMeanClsizePerTrack -> SetBorderMode(0);
    canvasMeanClsizePerTrack -> SetBorderSize(0);
    canvasMeanClsizePerTrack -> SetTickx(1);
    canvasMeanClsizePerTrack -> SetTicky(1);
    canvasMeanClsizePerTrack -> SetLeftMargin(0.15);
    canvasMeanClsizePerTrack -> SetBottomMargin(0.1518219);
    canvasMeanClsizePerTrack -> SetFrameBorderMode(0);
    canvasMeanClsizePerTrack -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);

    for (int iDataset = 0;iDataset < nDatasets;iDataset++) {
        histMeanClsizePerTrack[iDataset] -> Scale(1. / histMeanClsizePerTrack[iDataset] -> Integral());
        histMeanClsizePerTrack[iDataset] -> Draw("EH SAME");
    }

    TLatex *letexTitle = new TLatex();
    letexTitle -> SetTextSize(0.045);
    letexTitle -> SetNDC();
    letexTitle -> SetTextFont(42);
    //letexTitle -> DrawLatex(0.38, 0.82, "LHC23zs pass2, pp, #sqrt{#it{s}} = 13.6 TeV");
    letexTitle -> DrawLatex(0.45, 0.80, "MFT-MCH-MID tracks, 2.5 < #it{y} < 4");
    //letexTitle -> DrawLatex(0.58, 0.65, Form("#mu = %3.2f", histMeanClsizePerTrack -> GetMean()));
    //letexTitle -> DrawLatex(0.58, 0.60, Form("#sigma = %3.2f", histMeanClsizePerTrack -> GetStdDev()));

    TLegend *legendMeanClsizePerTrack = new TLegend(0.58,0.5,0.75,0.72," ","brNDC");
    SetLegend(legendMeanClsizePerTrack);
    for (int iDataset = 0;iDataset < nDatasets;iDataset++) {
        legendMeanClsizePerTrack -> AddEntry(histMeanClsizePerTrack[iDataset], Form("%s", datasets[iDataset].c_str()), "L");
    }
    legendMeanClsizePerTrack -> Draw("SAME");

    TCanvas *canvasClsizePerLayer = new TCanvas("canvasClsizePerLayer", "", 6000, 1200);
    canvasClsizePerLayer -> SetFillColor(0);
    canvasClsizePerLayer -> SetBorderMode(0);
    canvasClsizePerLayer -> SetBorderSize(0);
    canvasClsizePerLayer -> SetTickx(1);
    canvasClsizePerLayer -> SetTicky(1);
    canvasClsizePerLayer -> SetLeftMargin(0.15);
    canvasClsizePerLayer -> SetBottomMargin(0.1518219);
    canvasClsizePerLayer -> SetFrameBorderMode(0);
    canvasClsizePerLayer -> SetFrameBorderMode(0);
    canvasClsizePerLayer -> Divide(5, 2);

    for (int i = 0;i < 10;i++) {
        canvasClsizePerLayer -> cd(i+1);
        gPad -> SetLogy(1);
        for (int iDataset = 0;iDataset < nDatasets;iDataset++) {
            histClsizePerLayer[iDataset][i] -> Scale(1. / histClsizePerLayer[iDataset][i] -> Integral());
            histClsizePerLayer[iDataset][i] -> Draw("HE SAME");
        }
    }

    canvasClsizePerLayer -> SaveAs("clsize_per_layer.pdf");
    canvasMeanClsizePerTrack -> SaveAs("mean_clsize_per_track.pdf");



}
////////////////////////////////////////////////////////////////////////////////
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
}
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}