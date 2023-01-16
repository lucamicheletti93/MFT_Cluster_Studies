#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>

#include <gsl/gsl>
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
#include "TStopwatch.h"

#endif

void process_output_fast(){

    double mMeanClsizePerTrackRec = 0;
    int mClsizeRecLayer0 = 0;
    int mClsizeRecLayer1 = 0;
    int mClsizeRecLayer2 = 0;
    int mClsizeRecLayer3 = 0;
    int mClsizeRecLayer4 = 0;
    int mClsizeRecLayer5 = 0;
    int mClsizeRecLayer6 = 0;
    int mClsizeRecLayer7 = 0;
    int mClsizeRecLayer8 = 0;
    int mClsizeRecLayer9 = 0;
    int mPattIdRecLayer0 = 0;
    int mPattIdRecLayer1 = 0;
    int mPattIdRecLayer2 = 0;
    int mPattIdRecLayer3 = 0;
    int mPattIdRecLayer4 = 0;
    int mPattIdRecLayer5 = 0;
    int mPattIdRecLayer6 = 0;
    int mPattIdRecLayer7 = 0;
    int mPattIdRecLayer8 = 0;
    int mPattIdRecLayer9 = 0;
    double mEtaRec = 0;
    double mPhiRec = 0;

    int mPdgCodeGen = 0;
    bool mIsMIP = kFALSE;

    string fInNames[] = {"../output/MFTAssessmentData_LHC22r_529341_1000.root"};
    string fOutNames[] = {"../output/MFTAssessmentData_LHC22r_529341_1000_skimmed.root"};
    string treeNames[] = {"treeTrackClusterSizeData"};
    bool isMC[] = {kFALSE};
    string prodLabel[] = {"Data"};

    const int nProds = sizeof(treeNames) / sizeof(treeNames[0]);
    const int nLayers = 10;
    int iProd = 0;

    TH1F *histMeanClsizePerTrackRec[nProds];
    TH1F *histEtaRec[nProds];
    TH1F *histPhiRec[nProds];
    TH1F *histClsizePerLayer[nProds][nLayers];
    TH1F *histPattIdPerLayer[nProds][nLayers];
    TH2F *histEtaVsMeanClsizePerTrackRec[nProds];

    TStopwatch timer;
    for (auto fInName : fInNames) {
        TFile *fIn = new TFile(fInName.c_str(), "READ");
        iProd = find(fInNames, fInNames + nProds, fInName) - fInNames;
        std::cout << fInName << " " << iProd << std::endl;

        histMeanClsizePerTrackRec[iProd] = new TH1F("histMeanClsizePerTrackRec", "", 800, 0, 8);
        histMeanClsizePerTrackRec[iProd] -> SetDirectory(0);

        histEtaRec[iProd] = new TH1F("histEtaRec", "", 200, -10, 10);
        histEtaRec[iProd] -> SetDirectory(0);

        histPhiRec[iProd] = new TH1F("histPhiRec", "", 100, -2*TMath::Pi(), 2*TMath::Pi());
        histPhiRec[iProd] -> SetDirectory(0);

        histEtaVsMeanClsizePerTrackRec[iProd] = new TH2F("histEtaVsMeanClsizePerTrackRec", "", 800, 0, 8, 200, -10, 10);
        histEtaVsMeanClsizePerTrackRec[iProd] -> SetDirectory(0);

        for (int iLayer = 0;iLayer < 10;iLayer++) {
            histClsizePerLayer[iProd][iLayer] = new TH1F(Form("histClsizeRecLayer%i", iLayer), "", 1000, -0.5, 999.5);
            histClsizePerLayer[iProd][iLayer] -> SetDirectory(0);

            histPattIdPerLayer[iProd][iLayer] = new TH1F(Form("histPattIdRecLayer%i", iLayer), "", 10000, -0.5, 9999.5);
            histPattIdPerLayer[iProd][iLayer] -> SetDirectory(0);
        }

        TTree *tree = (TTree*) fIn -> Get(treeNames[iProd].c_str());
        tree -> SetBranchAddress("mMeanClsizePerTrackRec", &mMeanClsizePerTrackRec);
        tree -> SetBranchAddress("mClsizeRecLayer0", &mClsizeRecLayer0);
        tree -> SetBranchAddress("mClsizeRecLayer1", &mClsizeRecLayer1);
        tree -> SetBranchAddress("mClsizeRecLayer2", &mClsizeRecLayer2);
        tree -> SetBranchAddress("mClsizeRecLayer3", &mClsizeRecLayer3);
        tree -> SetBranchAddress("mClsizeRecLayer4", &mClsizeRecLayer4);
        tree -> SetBranchAddress("mClsizeRecLayer5", &mClsizeRecLayer5);
        tree -> SetBranchAddress("mClsizeRecLayer6", &mClsizeRecLayer6);
        tree -> SetBranchAddress("mClsizeRecLayer7", &mClsizeRecLayer7);
        tree -> SetBranchAddress("mClsizeRecLayer8", &mClsizeRecLayer8);
        tree -> SetBranchAddress("mClsizeRecLayer9", &mClsizeRecLayer9);
        tree -> SetBranchAddress("mPattIdRecLayer0", &mPattIdRecLayer0);
        tree -> SetBranchAddress("mPattIdRecLayer1", &mPattIdRecLayer1);
        tree -> SetBranchAddress("mPattIdRecLayer2", &mPattIdRecLayer2);
        tree -> SetBranchAddress("mPattIdRecLayer3", &mPattIdRecLayer3);
        tree -> SetBranchAddress("mPattIdRecLayer4", &mPattIdRecLayer4);
        tree -> SetBranchAddress("mPattIdRecLayer5", &mPattIdRecLayer5);
        tree -> SetBranchAddress("mPattIdRecLayer6", &mPattIdRecLayer6);
        tree -> SetBranchAddress("mPattIdRecLayer7", &mPattIdRecLayer7);
        tree -> SetBranchAddress("mPattIdRecLayer8", &mPattIdRecLayer8);
        tree -> SetBranchAddress("mPattIdRecLayer9", &mPattIdRecLayer9);
        tree -> SetBranchAddress("mEtaRec", &mEtaRec);
        tree -> SetBranchAddress("mPhiRec", &mPhiRec);

        if (isMC[iProd]) {
            tree -> SetBranchAddress("mPdgCodeGen", &mPdgCodeGen);
            tree -> SetBranchAddress("mIsMIP", &mIsMIP);
        }

        TStopwatch timer;
        timer.Start();

        int nEv = tree -> GetEntries();
        for (int i = 0;i < nEv;i++) {
            tree -> GetEntry(i);
            printf("Percentage %f \r", (double) i / (double) nEv);

            histEtaVsMeanClsizePerTrackRec[iProd] -> Fill(mMeanClsizePerTrackRec, mEtaRec);
            if (mEtaRec < -4. || mEtaRec > -2.5) {
                continue;
            }
            if (isMC[iProd]) {
                if (mIsMIP) {
                    histMeanClsizePerTrackRec[iProd] -> Fill(mMeanClsizePerTrackRec);
                }
            } else {
                histMeanClsizePerTrackRec[iProd] -> Fill(mMeanClsizePerTrackRec);
                histEtaRec[iProd] -> Fill(mEtaRec);
                histPhiRec[iProd] -> Fill(mPhiRec);
            }

            histClsizePerLayer[iProd][0] -> Fill(mClsizeRecLayer0);
            histClsizePerLayer[iProd][1] -> Fill(mClsizeRecLayer1);
            histClsizePerLayer[iProd][2] -> Fill(mClsizeRecLayer2);
            histClsizePerLayer[iProd][3] -> Fill(mClsizeRecLayer3);
            histClsizePerLayer[iProd][4] -> Fill(mClsizeRecLayer4);
            histClsizePerLayer[iProd][5] -> Fill(mClsizeRecLayer5);
            histClsizePerLayer[iProd][6] -> Fill(mClsizeRecLayer6);
            histClsizePerLayer[iProd][7] -> Fill(mClsizeRecLayer7);
            histClsizePerLayer[iProd][8] -> Fill(mClsizeRecLayer8);
            histClsizePerLayer[iProd][9] -> Fill(mClsizeRecLayer9);
            histPattIdPerLayer[iProd][0] -> Fill(mPattIdRecLayer0);
            histPattIdPerLayer[iProd][1] -> Fill(mPattIdRecLayer1);
            histPattIdPerLayer[iProd][2] -> Fill(mPattIdRecLayer2);
            histPattIdPerLayer[iProd][3] -> Fill(mPattIdRecLayer3);
            histPattIdPerLayer[iProd][4] -> Fill(mPattIdRecLayer4);
            histPattIdPerLayer[iProd][5] -> Fill(mPattIdRecLayer5);
            histPattIdPerLayer[iProd][6] -> Fill(mPattIdRecLayer6);
            histPattIdPerLayer[iProd][7] -> Fill(mPattIdRecLayer7);
            histPattIdPerLayer[iProd][8] -> Fill(mPattIdRecLayer8);
            histPattIdPerLayer[iProd][9] -> Fill(mPattIdRecLayer9);
        }
        timer.Stop();
        printf("RT=%7.3f s, Cpu=%7.3f s\n", timer.RealTime(), timer.CpuTime());
        fIn -> Close();
    }

    for (auto fOutName : fOutNames) {
        TFile *fOut = new TFile(fOutName.c_str(), "RECREATE");
        for (auto fOutName : fOutNames) {
            iProd = find(fOutNames, fOutNames + nProds, fOutName) - fOutNames;
            histMeanClsizePerTrackRec[iProd] -> Write();
            histEtaRec[iProd] -> Write();
            histPhiRec[iProd] -> Write();
            histEtaVsMeanClsizePerTrackRec[iProd] -> Write();

            for (int iLayer = 0;iLayer < 10;iLayer++) {
                histClsizePerLayer[iProd][iLayer] -> Write();
                histPattIdPerLayer[iProd][iLayer] -> Write();
            }
        }
        fOut -> Close();
    }
}

/*
void plot_results(){
    gStyle -> SetOptStat(0);
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

    TFile *fIn1 = new TFile("../output/MFTAssessmentMC.root", "READ");
    TH1F *mMFTTrackClusterSize1 = (TH1F*) fIn1 -> Get("mMFTTrackClusterSize");
    mMFTTrackClusterSize1 -> Scale(1. / mMFTTrackClusterSize1 -> Integral());
    mMFTTrackClusterSize1 -> SetLineColor(kRed);

    TFile *fIn2 = new TFile("../output/MFTAssessmentData_LHC22r_529341_1000.root", "READ");
    TH1F *mMFTTrackClusterSize2 = (TH1F*) fIn2 -> Get("mMFTTrackClusterSize");
    mMFTTrackClusterSize2 -> Scale(1. / mMFTTrackClusterSize2 -> Integral());
    mMFTTrackClusterSize2 -> SetLineColor(kBlack);

    TCanvas *canvasTest = new TCanvas("canvasTest", "", 600, 600);
    gPad -> SetLogx();
    gPad -> SetLogy();
    mMFTTrackClusterSize1 -> Draw("H");
    mMFTTrackClusterSize2 -> Draw("EPsame");


    TFile *fIn = new TFile("output.root", "READ");

    TH1F *histMeanClsizePerTrackRec_Data = (TH1F*) fIn -> Get("MeanClsizePerTrackRec_Data");
    histMeanClsizePerTrackRec_Data -> Rebin(4);
    histMeanClsizePerTrackRec_Data -> SetMarkerStyle(24);
    histMeanClsizePerTrackRec_Data -> SetMarkerColor(kBlack);
    histMeanClsizePerTrackRec_Data -> SetLineColor(kBlack);
    histMeanClsizePerTrackRec_Data -> Scale(1. / histMeanClsizePerTrackRec_Data -> Integral());

    TH1F *histClsizeLayer_Data[10];
    for (int iLayer = 0;iLayer < 10;iLayer++) {
        histClsizeLayer_Data[iLayer] = (TH1F*) fIn -> Get(Form("ClsizeLayer_%i_Data", iLayer));
        histClsizeLayer_Data[iLayer] -> Rebin(2);
        histClsizeLayer_Data[iLayer] -> SetMarkerStyle(24);
        histClsizeLayer_Data[iLayer] -> SetMarkerSize(0.8);
        histClsizeLayer_Data[iLayer] -> SetMarkerColor(kBlack);
        histClsizeLayer_Data[iLayer] -> SetLineColor(kBlack);
        histClsizeLayer_Data[iLayer] -> Scale(1. / histClsizeLayer_Data[iLayer] -> Integral());
    }

    std::cout << histMeanClsizePerTrackRec_Data -> GetMean() << std::endl;

    TH1F *histMeanClsizePerTrackRec_MC = (TH1F*) fIn -> Get("MeanClsizePerTrackRec_MC");
    histMeanClsizePerTrackRec_MC -> Rebin(4);
    histMeanClsizePerTrackRec_MC -> SetLineColor(kRed+1);
    histMeanClsizePerTrackRec_MC -> SetLineWidth(2);
    histMeanClsizePerTrackRec_MC -> Scale(1. / histMeanClsizePerTrackRec_MC -> Integral());

    TH1F *histClsizeLayer_MC[10];
    for (int iLayer = 0;iLayer < 10;iLayer++) {
        histClsizeLayer_MC[iLayer] = (TH1F*) fIn -> Get(Form("ClsizeLayer_%i_MC", iLayer));
        histClsizeLayer_MC[iLayer] -> Rebin(2);
        histClsizeLayer_MC[iLayer] -> SetLineColor(kRed+1);
        histClsizeLayer_MC[iLayer] -> SetLineWidth(2);
        histClsizeLayer_MC[iLayer] -> Scale(1. / histClsizeLayer_MC[iLayer] -> Integral());
    }

    std::cout << histMeanClsizePerTrackRec_MC -> GetMean() << std::endl;

    TCanvas *canvasMeanClsizePerTrackRec = new TCanvas("canvasMeanClsizePerTrackRec", "", 800, 600);
    gPad -> SetLogy();
    histMeanClsizePerTrackRec_Data -> Draw("EP");
    histMeanClsizePerTrackRec_MC -> Draw("Hsame");

    TCanvas *canvasClsizeLayer = new TCanvas("canvasClsizeLayer", "", 3000, 1200);
    canvasClsizeLayer -> Divide(5, 2);

    for (int iLayer = 0;iLayer < 10;iLayer++) {
        canvasClsizeLayer -> cd(iLayer+1);
        gPad -> SetLogx();
        gPad -> SetLogy();
        histClsizeLayer_Data[iLayer] -> Draw("EP");
        histClsizeLayer_MC[iLayer] -> Draw("Hsame");
        std::cout << histClsizeLayer_Data[iLayer] -> GetMean() << std::endl;
        std::cout << histClsizeLayer_MC[iLayer] -> GetMean() << std::endl;
    }

}
*/