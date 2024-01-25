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

void read_tree_fast() {
    //string dataset = "LHC23k6c_pass1";
    //string dataset = "LHC23zs_pass2";
    string dataset = "LHC23zzh_pass1_small";
    
    float fPosX, fPosY, fPosZ, fPt, fEta, fPhi, fFwdDcaX, fFwdDcaY, fChi2MatchMCHMID, fChi2MatchMCHMFT = -99999;
    uint64_t fMftClusterSizesAndTrackFlags;
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

    TFile *fOut = new TFile(Form("/Users/lucamicheletti/cernbox/MTF_Cluster_Data/%s/treeMftPidTraining.root", dataset.c_str()), "RECREATE");
    
    TTree *treeOut = new TTree("treeTrackClusterSizeData", "Tree with ClusterSize Data Info");
    treeOut -> Branch("mTrackType", &mTrackType, "mTrackType/I");
    treeOut -> Branch("mPt", &mPt, "mPt/F");
    treeOut -> Branch("mEta", &mEta, "mEta/F");
    treeOut -> Branch("mPhi", &mPhi, "mPhi/F");
    treeOut -> Branch("mSign", &mSign, "mSign/I");
    treeOut -> Branch("mMch2MCHMID", &mMch2MCHMID, "mMch2MCHMID/F");
    treeOut -> Branch("mMch2MFTMCH", &mMch2MFTMCH, "mMch2MFTMCH/F");
    treeOut -> Branch("mFwdDcaX", &mFwdDcaX, "mFwdDcaX/F");
    treeOut -> Branch("mFwdDcaY", &mFwdDcaY, "mFwdDcaY/F");
    treeOut -> Branch("mClsizeLayer0", &mClsizeLayer0);
    treeOut -> Branch("mClsizeLayer1", &mClsizeLayer1);
    treeOut -> Branch("mClsizeLayer2", &mClsizeLayer2);
    treeOut -> Branch("mClsizeLayer3", &mClsizeLayer3);
    treeOut -> Branch("mClsizeLayer4", &mClsizeLayer4);
    treeOut -> Branch("mClsizeLayer5", &mClsizeLayer5);
    treeOut -> Branch("mClsizeLayer6", &mClsizeLayer6);
    treeOut -> Branch("mClsizeLayer7", &mClsizeLayer7);
    treeOut -> Branch("mClsizeLayer8", &mClsizeLayer8);
    treeOut -> Branch("mClsizeLayer9", &mClsizeLayer9);
    treeOut -> Branch("mMeanClsizePerTrack", &mMeanClsizePerTrack);

    Printf("[info] Processing /Users/lucamicheletti/cernbox/MTF_Cluster_Data/%s/AO2D.root ...", dataset.c_str());
    string pathToFile = Form("/Users/lucamicheletti/cernbox/MTF_Cluster_Data/%s/AO2D.root", dataset.c_str());
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
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

        for (int iEntry = 0;iEntry < treeIn -> GetEntries();iEntry++) {
            treeIn -> GetEntry(iEntry);


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
            mMeanClsizePerTrack = meanClSize;

            treeOut -> Fill();
        }
    }
    fIn -> Close();

    Printf("[info] Saving output in /Users/lucamicheletti/cernbox/MTF_Cluster_Data/%s/treeMftPidTraining.root", dataset.c_str());
    fOut -> cd();
    treeOut -> Write();
    fOut -> Close();
}