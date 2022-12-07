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

void read_tree(){

    double mMeanClsizePerTrackRec = 0;
    int mPdgCodeGen = 0;

    string fileNames[] = {"../output/MFTAssessmentData_LHC22n_526195_0430.root", "../output/MFTAssessmentMC.root"};
    string treeNames[] = {"treeTrackClusterSizeData", "treeTrackClusterSize"};
    bool isMC[] = {kFALSE, kTRUE};

    int n = sizeof(treeNames) / sizeof(treeNames[0]);
    int index = 0;

    TH1F *histList[n];

    TStopwatch timer;
    for (auto fileName : fileNames) {
        TFile *fIn = new TFile(fileName.c_str(), "READ");
        index = find(fileNames, fileNames + n, fileName) - fileNames;
        std::cout << fileName << " " << index << std::endl;

        histList[index] = new TH1F(treeNames[index].c_str(), "", 80, 0, 8);
        histList[index] -> SetDirectory(0);

        TTree *tree = (TTree*) fIn -> Get(treeNames[index].c_str());
        tree -> SetBranchAddress("mMeanClsizePerTrackRec", &mMeanClsizePerTrackRec);
        if (isMC[index]) {
            tree -> SetBranchAddress("mPdgCodeGen", &mPdgCodeGen);
        }

        TStopwatch timer;
        timer.Start();

        int nEv = tree -> GetEntries();
        for (int i = 0;i < nEv;i++) {
            tree -> GetEntry(i);
            //histMeanClsizePerTrackRec -> Fill(mMeanClsizePerTrackRec);
            printf("Percentage %f \r", (double) i / (double) nEv);
            if (isMC[index]) {
                if (TMath::Abs(mPdgCodeGen) == 1000020030) {
                    histList[index] -> Fill(mMeanClsizePerTrackRec);
                }
            } else {
                histList[index] -> Fill(mMeanClsizePerTrackRec);
            }
        }
        timer.Stop();
        printf("RT=%7.3f s, Cpu=%7.3f s\n", timer.RealTime(), timer.CpuTime());
        fIn -> Close();
    }

    TFile *fOut = new TFile("output.root", "RECREATE");
    for (auto fileName : fileNames) {
        index = find(fileNames, fileNames + n, fileName) - fileNames;
        histList[index] -> Write();
    }
    fOut -> Close();
}

void plot_results(){
    TFile *fIn = new TFile("output.root", "READ");

    TH1F *histData = (TH1F*) fIn -> Get("treeTrackClusterSizeData");
    histData -> SetMarkerStyle(20);
    histData -> SetMarkerColor(kBlack);
    histData -> SetLineColor(kBlack);

    TH1F *histMC = (TH1F*) fIn -> Get("treeTrackClusterSize");
    histMC -> SetLineColor(kRed+1);

    TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
    histData -> Draw("EP");
    histData -> Draw("Hsame");

}