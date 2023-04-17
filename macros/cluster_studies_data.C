#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>

#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/GlobalFwdTrack.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DataFormatsMCH/TrackMCH.h"

#include "ITStracking/IOUtils.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

#include "Steer/MCKinematicsReader.h"

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
#include "CommonDataFormat/RangeReference.h"
//#include "DetectorsVertexing/DCAFitterN.h"

//#include "CCDB/BasicCCDBManager.h"
//#include "CCDB/CCDBTimeStampUtils.h"
//#include "Framework/CCDBParamSpec.h"

#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using o2::itsmft::CompClusterExt;
using CompClusterExt = o2::itsmft::CompClusterExt;
using ITSCluster = o2::BaseCluster<float>;
using Vec3 = ROOT::Math::SVector<double, 3>;
using MCLabCont = o2::dataformats::MCTruthContainer<o2::MCCompLabel>;
o2::itsmft::ChipMappingMFT mMFTChipMapper;


void cluster_studies_data(){

    //string pathIn = "cross_check/data";
    //string pathIn = "cross_check/data/reco_test/reco_files";
    string pathIn = "data/523308_apass3_relval_cpu2_debug_0050";

    TH1F *mTrackClusterSize = new TH1F("mMFTTrackClusterSize", "Cluster Size Per Track; # clusters; # entries", 100, 0.5, 100.5);
    TH1F *mTrackMeanClusterSize = new TH1F("mMFTTrackMeanClusterSize", "Mean Cluster Size Per Track; # clusters; # entries", 50, 0.5, 10.5);
    TH1F *mTrackEta = new TH1F("mTrackEta", "Eta Per Track; # Eta; # entries", 50, -4, -2);

    TH1F *mGlobalTrackMeanClusterSize = new TH1F("mGlobalTrackMeanClusterSize", "Mean Cluster Size Per Global Track; # clusters; # entries", 50, 0.5, 10.5);
    TH1F *mGlobalTrackPt = new TH1F("mGlobalTrackPt", "Pt Per Global Track; # Pt; # entries", 100, 0, 20);
    TH1F *mGlobalTrackEta = new TH1F("mGlobalTrackEta", "Eta Per Global Track; # Eta; # entries", 50, -4, -2);
    TH2F *mGlobalTrackEtaMFTMCHMatchingChi2 = new TH2F("mGlobalTrackEtaMFTMCHMatchingChi2", "Eta vs MFT-MCH #chi^{2}; # Eta; # MFT-MCH #chi^{2}", 50, -4, -2, 100, 0., 100e4);
    TH1F *mGlobalTrackEtaMFT = new TH1F("mGlobalTrackEtaMFT", "Eta MFT Per Global Track; # Eta; # entries", 50, -4, -2);
    TH1F *mGlobalTrackEtaMCH = new TH1F("mGlobalTrackEtaMCH", "Eta MCH Per Global Track; # Eta; # entries", 50, -4, -2);

    TH1F *mGlobalTrackMeanClusterSizeMuon = new TH1F("mGlobalTrackMeanClusterSizeMuon", "Mean Cluster Size Per Global Track [muon]; # clusters; # entries", 50, 0.5, 10.5);
    TH1F *mGlobalTrackPtMuon = new TH1F("mGlobalTrackPtMuon", "Pt Per Global Track [muon]; # Pt; # entries", 100, 0, 20);
    TH1F *mGlobalTrackEtaMuon = new TH1F("mGlobalTrackEtaMuon", "Eta Per Global Track [muon]; # Eta; # entries", 50, -4, -2);
    TH2F *mGlobalTrackEtaMFTMCHMatchingChi2Muon = new TH2F("mGlobalTrackEtaMFTMCHMatchingChi2Muon", "Eta vs MFT-MCH #chi^{2} [muon]; # Eta; # MFT-MCH #chi^{2}", 50, -4, -2, 100, 0., 100e4);
    TH1F *mGlobalTrackEtaMFTMuon = new TH1F("mGlobalTrackEtaMFTMuon", "Eta MFT Per Global Track [muon]; # Eta; # entries", 50, -4, -2);
    TH1F *mGlobalTrackEtaMCHMuon = new TH1F("mGlobalTrackEtaMCHMuon", "Eta MCH Per Global Track [muon]; # Eta; # entries", 50, -4, -2);

    TH1F *mMIDMatchingChi2 = new TH1F("mMIDMatchingChi2", "MID matching Chi2; # MID matching Chi2; # entries", 100, -5, 20.);
    TH1F *mMFTMCHMatchingChi2 = new TH1F("mMFTMCHMatchingChi2", "MFT-MCH matching Chi2; # MFT-MCH matching Chi2; # entries", 5000, 0., 100e4);

    //auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_1320237.root"); // March 2022
    auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_0747987.root"); // December 2022
    o2::itsmft::TopologyDictionary mDictionary = *(reinterpret_cast<o2::itsmft::TopologyDictionary *>(fIn.Get("ccdb_object")));

    printf("Reading the input files\n");
    TFile *fMFTClusters = new TFile(Form("%s/mftclusters.root", pathIn.c_str()), "READ");
    TFile *fMFTTracks = new TFile(Form("%s/mfttracks.root", pathIn.c_str()), "READ");
    TFile *fMCHTracks = new TFile(Form("%s/mchtracks.root", pathIn.c_str()), "READ");
    TFile *fGlobalTracks = new TFile(Form("%s/globalfwdtracks.root", pathIn.c_str()), "READ");
    
    TTree *treeMFTClusters = (TTree*) fMFTClusters -> Get("o2sim");
    TTree *treeMFTTracks = (TTree*) fMFTTracks -> Get("o2sim");
    TTree *treeMCHTracks = (TTree*) fMCHTracks -> Get("o2sim");
    TTree *treeGlobalTracks = (TTree*) fGlobalTracks -> Get("GlobalFwdTracks");

    TFile *fOut = new TFile(Form("%s/cross_check_data.root", pathIn.c_str()), "RECREATE");
    // Global tracks
    int mChargeGlobalRec = -999;
    float mPtGlobalRec = -999;
    float mEtaGlobalRec = -999;
    float mPhiGlobalRec = -999;
    // MCH tracks
    int mChargeMCHRec = -999;
    float mPtMCHRec = -999;
    float mEtaMCHRec = -999;
    float mPhiMCHRec = -999;
    // MFT tracks
    int mChargeMFTRec = -999;
    float mPtMFTRec = -999;
    float mEtaMFTRec = -999;
    float mPhiMFTRec = -999;
    // Matching quality
    float mMch2MCHMIDRec = -999;
    float mMch2MFTMCHRec = -999;

    int mClsizeRecLayer0 = -999;
    int mClsizeRecLayer1 = -999;
    int mClsizeRecLayer2 = -999;
    int mClsizeRecLayer3 = -999;
    int mClsizeRecLayer4 = -999;
    int mClsizeRecLayer5 = -999;
    int mClsizeRecLayer6 = -999;
    int mClsizeRecLayer7 = -999;
    int mClsizeRecLayer8 = -999;
    int mClsizeRecLayer9 = -999;
    int mPattIdRecLayer0 = -999;
    int mPattIdRecLayer1 = -999;
    int mPattIdRecLayer2 = -999;
    int mPattIdRecLayer3 = -999;
    int mPattIdRecLayer4 = -999;
    int mPattIdRecLayer5 = -999;
    int mPattIdRecLayer6 = -999;
    int mPattIdRecLayer7 = -999;
    int mPattIdRecLayer8 = -999;
    int mPattIdRecLayer9 = -999;
    float mMeanClsizePerTrackRec = -999;

    TTree *mTreeTrackClusterSizeData = new TTree("treeTrackClusterSizeData", "Tree with ClusterSize Data Info");
    mTreeTrackClusterSizeData -> Branch("mChargeGlobalRec", &mChargeGlobalRec, "mChargeGlobalRec/I");
    mTreeTrackClusterSizeData -> Branch("mPtGlobalRec", &mPtGlobalRec, "mPtGlobalRec/F");
    mTreeTrackClusterSizeData -> Branch("mEtaGlobalRec", &mEtaGlobalRec, "mEtaGlobalRec/F");
    mTreeTrackClusterSizeData -> Branch("mPhiGlobalRec", &mPhiGlobalRec, "mPhiGlobalRec/F");

    mTreeTrackClusterSizeData -> Branch("mChargeMCHRec", &mChargeMCHRec, "mChargeMCHRec/I");
    mTreeTrackClusterSizeData -> Branch("mPtMCHRec", &mPtMCHRec, "mPtMCHRec/F");
    mTreeTrackClusterSizeData -> Branch("mEtaMCHRec", &mEtaMCHRec, "mEtaMCHRec/F");
    mTreeTrackClusterSizeData -> Branch("mPhiMCHRec", &mPhiMCHRec, "mPhiMCHRec/F");

    mTreeTrackClusterSizeData -> Branch("mChargeMFTRec", &mChargeMFTRec, "mChargeMFTRec/I");
    mTreeTrackClusterSizeData -> Branch("mPtMFTRec", &mPtMFTRec, "mPtMFTRec/F");
    mTreeTrackClusterSizeData -> Branch("mEtaMFTRec", &mEtaMFTRec, "mEtaMFTRec/F");
    mTreeTrackClusterSizeData -> Branch("mPhiMFTRec", &mPhiMFTRec, "mPhiMFTRec/F");

    mTreeTrackClusterSizeData -> Branch("mMch2MCHMIDRec", &mMch2MCHMIDRec, "mMch2MCHMIDRec/F");
    mTreeTrackClusterSizeData -> Branch("mMch2MFTMCHRec", &mMch2MFTMCHRec, "mMch2MFTMCHRec/F");

    // Clsize per layer
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer0", &mClsizeRecLayer0);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer1", &mClsizeRecLayer1);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer2", &mClsizeRecLayer2);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer3", &mClsizeRecLayer3);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer4", &mClsizeRecLayer4);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer5", &mClsizeRecLayer5);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer6", &mClsizeRecLayer6);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer7", &mClsizeRecLayer7);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer8", &mClsizeRecLayer8);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer9", &mClsizeRecLayer9);
    // Pattern ID per layer
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer0", &mPattIdRecLayer0);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer1", &mPattIdRecLayer1);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer2", &mPattIdRecLayer2);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer3", &mPattIdRecLayer3);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer4", &mPattIdRecLayer4);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer5", &mPattIdRecLayer5);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer6", &mPattIdRecLayer6);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer7", &mPattIdRecLayer7);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer8", &mPattIdRecLayer8);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer9", &mPattIdRecLayer9);
    // Mean cluster size
    mTreeTrackClusterSizeData -> Branch("mMeanClsizePerTrackRec", &mMeanClsizePerTrackRec);

    // get MFT tracks
    std::vector<o2::mft::TrackMFT> *mMFTTracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTTracksROF = nullptr;
    std::vector<int> *mMFTTrackClusIdx = nullptr;

    // get MFT clusters
    std::vector<o2::itsmft::CompClusterExt> *mMFTClusters = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTClustersROF = nullptr;
    std::vector<unsigned char> *mMFTClusterPatterns = nullptr;

    std::vector<int> ClTrackID;
    std::vector<int> ClEvID;
    std::vector<int> ClSiezes;

    // get MCH tracks
    std::vector<o2::mch::TrackMCH> *mMCHTracks = nullptr;

    // get GlobalFwdTracks
    std::vector<o2::dataformats::GlobalFwdTrack> *mGlobalFwdTracks = nullptr;

    treeMFTTracks -> SetBranchAddress("MFTTrack", &mMFTTracks);
    treeMFTTracks -> SetBranchAddress("MFTTracksROF", &mMFTTracksROF);
    treeMFTTracks -> SetBranchAddress("MFTTrackClusIdx", &mMFTTrackClusIdx);
    
    treeMFTClusters -> SetBranchAddress("MFTClusterComp", &mMFTClusters);
    treeMFTClusters -> SetBranchAddress("MFTClustersROF", &mMFTClustersROF);
    treeMFTClusters -> SetBranchAddress("MFTClusterPatt", &mMFTClusterPatterns);

    treeMCHTracks -> SetBranchAddress("tracks", &mMCHTracks);

    treeGlobalTracks -> SetBranchAddress("fwdtracks", &mGlobalFwdTracks);
    
    std::cout << treeMFTClusters -> GetEntriesFast() << std::endl;
    for (int frame = 0; frame < treeMFTClusters -> GetEntriesFast(); frame++) {
        if (!treeMFTClusters -> GetEvent(frame) || !treeMFTTracks -> GetEvent(frame) || !treeMCHTracks -> GetEvent(frame) || !treeGlobalTracks -> GetEvent(frame)) {
        ///if (!treeMFTClusters -> GetEvent(frame) || !treeMFTTracks -> GetEvent(frame)) {
            LOG(info) << "Skipping frame: " << frame;
            continue;
        }
        //cout << "N clusters = " << mMFTClusters -> size() << std::endl;

        std::vector<o2::itsmft::ClusterPattern> pattVec;
        pattVec.reserve(mMFTClusters -> size());
        auto pattItMFT = mMFTClusterPatterns -> begin();

        for (unsigned int iClus{0}; iClus < mMFTClusters -> size(); ++iClus) {
            auto &clus = mMFTClusters -> at(iClus);
            auto pattID = clus.getPatternID();
            int npix;
            o2::itsmft::ClusterPattern patt;
            if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mDictionary.isGroup(pattID)) {
                patt.acquirePattern(pattItMFT);
                npix = patt.getNPixels();
            } else {
                npix = mDictionary.getNpixels(pattID);
                patt = mDictionary.getPattern(pattID);
            }
            mTrackClusterSize -> Fill(npix);
            pattVec.push_back(patt);
        }

        for (unsigned int iTrack{0}; iTrack < mMFTTracks -> size(); ++iTrack) {
            auto &oneTrack = mMFTTracks -> at(iTrack);
            auto ncls = oneTrack.getNumberOfPoints();
            auto offset = oneTrack.getExternalClusterIndexOffset();

            float mean = 0, norm = 0;

            for (int icls = 0; icls < ncls; ++icls) {
                auto clsEntry = mMFTTrackClusIdx -> at(offset + icls);
                auto &oneCluster = mMFTClusters -> at(clsEntry);
                auto &patt = pattVec.at(clsEntry);
                int npix = patt.getNPixels();
                //std::cout << npix << " + ";

                mean += npix;
                norm += 1;
            }
            mean /= norm;
            //std::cout << " -> MEAN CLs = " << mMFTTracks -> size() << std::endl;
            mTrackMeanClusterSize -> Fill(mean);
            mTrackEta -> Fill(oneTrack.getEta());
        }

        std::cout << " -> Global tracks size = " << mGlobalFwdTracks -> size() << std::endl;


        for (unsigned int iGlobalFwdTrack{0}; iGlobalFwdTrack < mGlobalFwdTracks -> size(); ++iGlobalFwdTrack) {
            auto &oneGlobalFwdTrack = mGlobalFwdTracks -> at(iGlobalFwdTrack);
            int oneMftTrackId = oneGlobalFwdTrack.getMFTTrackID();
            int oneMchTrackId = oneGlobalFwdTrack.getMCHTrackID();

            auto &oneMftTrack = mMFTTracks -> at(oneMftTrackId);
            auto &oneMchTrack = mMCHTracks -> at(oneMchTrackId);

            auto ncls = oneMftTrack.getNumberOfPoints();
            auto offset = oneMftTrack.getExternalClusterIndexOffset();
            //cout << "TRACK ID = " << oneMftTrackId << "; N points = " << ncls << endl;

            float mean = 0, norm = 0;
            int chargeGlobal = oneGlobalFwdTrack.getCharge();
            float ptGlobal = oneGlobalFwdTrack.getPt();
            float etaGlobal = oneGlobalFwdTrack.getEta();
            float phiGlobal = oneGlobalFwdTrack.getPhi();

            float xMCH = oneMchTrack.getX();
            float yMCH = oneMchTrack.getY();
            float pMCH = oneMchTrack.getP();
            float pxMCH = oneMchTrack.getPz();
            float pyMCH = oneMchTrack.getPz();
            float pzMCH = oneMchTrack.getPz();

            int chargeMCH = oneMchTrack.getSign();
            float ptMCH = TMath::Sqrt(pxMCH * pxMCH + pyMCH * pyMCH);
            float etaMCH = 0.5 * TMath::Log((pMCH + pzMCH) / (pMCH - pzMCH));
            float phiMCH = TMath::ATan2(yMCH, xMCH);

            int chargeMFT = oneMftTrack.getCharge();
            float ptMFT = oneMftTrack.getPt();
            float etaMFT = oneMftTrack.getEta();
            float phiMFT = oneMftTrack.getPhi();

            float ch2MCHMID = oneGlobalFwdTrack.getMIDMatchingChi2();
            float ch2MFTMCH = oneGlobalFwdTrack.getMFTMCHMatchingChi2();

            mChargeGlobalRec = chargeGlobal;
            mPtGlobalRec = ptGlobal;
            mEtaGlobalRec = etaGlobal;
            mPhiGlobalRec = phiGlobal;

            mChargeMCHRec = chargeMCH;
            mPtMCHRec = ptMCH;
            mEtaMCHRec = etaMCH;
            mPhiMCHRec = phiMCH;

            mChargeMFTRec = chargeMFT;
            mPtMFTRec = ptMFT;
            mEtaMFTRec = etaMFT;
            mPhiMFTRec = phiMFT;

            mMch2MCHMIDRec = ch2MCHMID;
            mMch2MFTMCHRec = ch2MFTMCH;

            //std::cout << mChargeGlobalRec << " " << mPtGlobalRec << " " << mEtaGlobalRec << " " << mPhiGlobalRec << " " << mMch2MCHMIDRec << " " << mMch2MFTMCHRec << std::endl;

            for (int icls = 0; icls < ncls; ++icls) {
                auto clsEntry = mMFTTrackClusIdx -> at(offset + icls);
                auto &oneCluster = mMFTClusters -> at(clsEntry);
                int clsLayer = mMFTChipMapper.chip2Layer(oneCluster.getChipID());
                auto &patt = pattVec.at(clsEntry);
                auto pattID = oneCluster.getPatternID();
                int npix = patt.getNPixels();
                //std::cout << npix << " (" << clsLayer << ") + ";

                switch(clsLayer) {
                    case 0:
                        mClsizeRecLayer0 = npix;
                        mPattIdRecLayer0 = pattID;
                        break;
                    case 1:
                        mClsizeRecLayer1 = npix;
                        mPattIdRecLayer1 = pattID;
                        break;
                    case 2:
                        mClsizeRecLayer2 = npix;
                        mPattIdRecLayer2 = pattID;
                        break;
                    case 3:
                        mClsizeRecLayer3 = npix;
                        mPattIdRecLayer3 = pattID;
                        break;
                    case 4:
                        mClsizeRecLayer4 = npix;
                        mPattIdRecLayer4 = pattID;
                        break;
                    case 5:
                        mClsizeRecLayer5 = npix;
                        mPattIdRecLayer5 = pattID;
                        break;
                    case 6:
                        mClsizeRecLayer6 = npix;
                        mPattIdRecLayer6 = pattID;
                        break;
                    case 7:
                        mClsizeRecLayer7 = npix;
                        mPattIdRecLayer7 = pattID;
                        break;
                    case 8:
                        mClsizeRecLayer8 = npix;
                        mPattIdRecLayer8 = pattID;
                        break;
                    case 9:
                        mClsizeRecLayer9 = npix;
                        mPattIdRecLayer9 = pattID;
                        break;
                }

                mean += npix;
                norm += 1;
            }
            mean /= norm;
            mMeanClsizePerTrackRec = mean;

            mMIDMatchingChi2 -> Fill(ch2MCHMID);
            mMFTMCHMatchingChi2 -> Fill(ch2MFTMCH);
            
            mGlobalTrackMeanClusterSize -> Fill(mean);
            mGlobalTrackPt -> Fill(ptGlobal);
            mGlobalTrackEta -> Fill(etaGlobal);
            mGlobalTrackEtaMFT -> Fill(etaMFT);
            mGlobalTrackEtaMCH -> Fill(etaMCH);
            mGlobalTrackEtaMFTMCHMatchingChi2 -> Fill(etaGlobal, ch2MFTMCH);

            if (ch2MCHMID > 0 && TMath::Abs(etaMCH) > 2.5 && TMath::Abs(etaMCH) < 4) {
                mGlobalTrackMeanClusterSizeMuon -> Fill(mean);
                mGlobalTrackPtMuon -> Fill(ptGlobal);
                mGlobalTrackEtaMuon -> Fill(etaGlobal);
                mGlobalTrackEtaMFTMuon -> Fill(etaMFT);
                mGlobalTrackEtaMCHMuon -> Fill(etaMCH);
                mGlobalTrackEtaMFTMCHMatchingChi2Muon -> Fill(etaGlobal, ch2MFTMCH);
            }

            // Fill the tree
            mTreeTrackClusterSizeData -> Fill();

            mChargeGlobalRec = -999;
            mPtGlobalRec = -999;
            mEtaGlobalRec = -999;
            mPhiGlobalRec = -999;
            
            mChargeMCHRec = -999;
            mPtMCHRec = -999;
            mEtaMCHRec = -999;
            mPhiMCHRec = -999;

            mChargeMFTRec = -999;
            mPtMFTRec = -999;
            mEtaMFTRec = -999;
            mPhiMFTRec = -999;
            // Clsize per layer
            mClsizeRecLayer0 = -999;
            mClsizeRecLayer1 = -999;
            mClsizeRecLayer2 = -999;
            mClsizeRecLayer3 = -999;
            mClsizeRecLayer4 = -999;
            mClsizeRecLayer5 = -999;
            mClsizeRecLayer6 = -999;
            mClsizeRecLayer7 = -999;
            mClsizeRecLayer8 = -999;
            mClsizeRecLayer9 = -999;
            // Pattern ID per layer
            mPattIdRecLayer0 = -999;
            mPattIdRecLayer1 = -999;
            mPattIdRecLayer2 = -999;
            mPattIdRecLayer3 = -999;
            mPattIdRecLayer4 = -999;
            mPattIdRecLayer5 = -999;
            mPattIdRecLayer6 = -999;
            mPattIdRecLayer7 = -999;
            mPattIdRecLayer8 = -999;
            mPattIdRecLayer9 = -999;
        }
    }

    fOut -> cd();
    mTrackClusterSize -> Write();
    mTrackMeanClusterSize -> Write();
    mTrackEta -> Write();

    mGlobalTrackMeanClusterSize -> Write();
    mGlobalTrackPt -> Write();
    mGlobalTrackEta -> Write();
    mGlobalTrackEtaMFT -> Write();
    mGlobalTrackEtaMCH -> Write();
    mGlobalTrackEtaMFTMCHMatchingChi2 -> Write();

    mGlobalTrackMeanClusterSizeMuon -> Write();
    mGlobalTrackPtMuon -> Write();
    mGlobalTrackEtaMuon -> Write();
    mGlobalTrackEtaMFTMuon -> Write();
    mGlobalTrackEtaMCHMuon -> Write();
    mGlobalTrackEtaMFTMCHMatchingChi2Muon -> Write();

    mMIDMatchingChi2 -> Write();
    mMFTMCHMatchingChi2 -> Write();

    mTreeTrackClusterSizeData -> Write();
    fOut -> Close();
}









void mft_cluster_studies_data(){

    //string pathIn = "cross_check/data";
    //string pathIn = "cross_check/data/reco_test/reco_files";
    string pathIn = "data/523308_apass3_relval_cpu2_debug_0040";

    TH1F *mTrackClusterSize = new TH1F("mMFTTrackClusterSize", "Cluster Size Per Track; # clusters; # entries", 100, 0.5, 100.5);
    TH1F *mTrackMeanClusterSize = new TH1F("mMFTTrackMeanClusterSize", "Mean Cluster Size Per Track; # clusters; # entries", 50, 0.5, 10.5);
    TH1F *mTrackEta = new TH1F("mTrackEta", "Eta Per Track; # Eta; # entries", 50, -4, -2);

    TH1F *mGlobalTrackMeanClusterSize = new TH1F("mGlobalTrackMeanClusterSize", "Mean Cluster Size Per Global Track; # clusters; # entries", 50, 0.5, 10.5);
    TH1F *mGlobalTrackPt = new TH1F("mGlobalTrackPt", "Pt Per Global Track; # Pt; # entries", 100, 0, 20);
    TH1F *mGlobalTrackEta = new TH1F("mGlobalTrackEta", "Eta Per Global Track; # Eta; # entries", 50, -4, -2);
    TH2F *mGlobalTrackEtaMFTMCHMatchingChi2 = new TH2F("mGlobalTrackEtaMFTMCHMatchingChi2", "Eta vs MFT-MCH #chi^{2}; # Eta; # MFT-MCH #chi^{2}", 50, -4, -2, 100, 0., 100e4);
    TH1F *mGlobalTrackEtaMFT = new TH1F("mGlobalTrackEtaMFT", "Eta MFT Per Global Track; # Eta; # entries", 50, -4, -2);
    TH1F *mGlobalTrackEtaMCH = new TH1F("mGlobalTrackEtaMCH", "Eta MCH Per Global Track; # Eta; # entries", 50, -4, -2);

    TH1F *mGlobalTrackMeanClusterSizeMuon = new TH1F("mGlobalTrackMeanClusterSizeMuon", "Mean Cluster Size Per Global Track [muon]; # clusters; # entries", 50, 0.5, 10.5);
    TH1F *mGlobalTrackPtMuon = new TH1F("mGlobalTrackPtMuon", "Pt Per Global Track [muon]; # Pt; # entries", 100, 0, 20);
    TH1F *mGlobalTrackEtaMuon = new TH1F("mGlobalTrackEtaMuon", "Eta Per Global Track [muon]; # Eta; # entries", 50, -4, -2);
    TH2F *mGlobalTrackEtaMFTMCHMatchingChi2Muon = new TH2F("mGlobalTrackEtaMFTMCHMatchingChi2Muon", "Eta vs MFT-MCH #chi^{2} [muon]; # Eta; # MFT-MCH #chi^{2}", 50, -4, -2, 100, 0., 100e4);
    TH1F *mGlobalTrackEtaMFTMuon = new TH1F("mGlobalTrackEtaMFTMuon", "Eta MFT Per Global Track [muon]; # Eta; # entries", 50, -4, -2);
    TH1F *mGlobalTrackEtaMCHMuon = new TH1F("mGlobalTrackEtaMCHMuon", "Eta MCH Per Global Track [muon]; # Eta; # entries", 50, -4, -2);

    TH1F *mMIDMatchingChi2 = new TH1F("mMIDMatchingChi2", "MID matching Chi2; # MID matching Chi2; # entries", 100, -5, 20.);
    TH1F *mMFTMCHMatchingChi2 = new TH1F("mMFTMCHMatchingChi2", "MFT-MCH matching Chi2; # MFT-MCH matching Chi2; # entries", 5000, 0., 100e4);

    //auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_1320237.root"); // March 2022
    auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_0747987.root"); // December 2022
    o2::itsmft::TopologyDictionary mDictionary = *(reinterpret_cast<o2::itsmft::TopologyDictionary *>(fIn.Get("ccdb_object")));

    printf("Reading the input files\n");
    TFile *fMFTClusters = new TFile(Form("%s/mftclusters.root", pathIn.c_str()), "READ");
    TFile *fMFTTracks = new TFile(Form("%s/mfttracks.root", pathIn.c_str()), "READ");
    
    TTree *treeMFTClusters = (TTree*) fMFTClusters -> Get("o2sim");
    TTree *treeMFTTracks = (TTree*) fMFTTracks -> Get("o2sim");

    TFile *fOut = new TFile(Form("%s/cross_check_data_mft_only.root", pathIn.c_str()), "RECREATE");
    // Global tracks
    int mChargeGlobalRec = -999;
    float mPtGlobalRec = -999;
    float mEtaGlobalRec = -999;
    float mPhiGlobalRec = -999;
    // MCH tracks
    int mChargeMCHRec = -999;
    float mPtMCHRec = -999;
    float mEtaMCHRec = -999;
    float mPhiMCHRec = -999;
    // MFT tracks
    int mChargeMFTRec = -999;
    float mPtMFTRec = -999;
    float mEtaMFTRec = -999;
    float mPhiMFTRec = -999;
    // Matching quality
    float mMch2MCHMIDRec = -999;
    float mMch2MFTMCHRec = -999;

    int mClsizeRecLayer0 = -999;
    int mClsizeRecLayer1 = -999;
    int mClsizeRecLayer2 = -999;
    int mClsizeRecLayer3 = -999;
    int mClsizeRecLayer4 = -999;
    int mClsizeRecLayer5 = -999;
    int mClsizeRecLayer6 = -999;
    int mClsizeRecLayer7 = -999;
    int mClsizeRecLayer8 = -999;
    int mClsizeRecLayer9 = -999;
    int mPattIdRecLayer0 = -999;
    int mPattIdRecLayer1 = -999;
    int mPattIdRecLayer2 = -999;
    int mPattIdRecLayer3 = -999;
    int mPattIdRecLayer4 = -999;
    int mPattIdRecLayer5 = -999;
    int mPattIdRecLayer6 = -999;
    int mPattIdRecLayer7 = -999;
    int mPattIdRecLayer8 = -999;
    int mPattIdRecLayer9 = -999;
    float mMeanClsizePerTrackRec = -999;

    TTree *mTreeTrackClusterSizeData = new TTree("treeTrackClusterSizeData", "Tree with ClusterSize Data Info");
    mTreeTrackClusterSizeData -> Branch("mChargeGlobalRec", &mChargeGlobalRec, "mChargeGlobalRec/I");
    mTreeTrackClusterSizeData -> Branch("mPtGlobalRec", &mPtGlobalRec, "mPtGlobalRec/F");
    mTreeTrackClusterSizeData -> Branch("mEtaGlobalRec", &mEtaGlobalRec, "mEtaGlobalRec/F");
    mTreeTrackClusterSizeData -> Branch("mPhiGlobalRec", &mPhiGlobalRec, "mPhiGlobalRec/F");

    mTreeTrackClusterSizeData -> Branch("mChargeMCHRec", &mChargeMCHRec, "mChargeMCHRec/I");
    mTreeTrackClusterSizeData -> Branch("mPtMCHRec", &mPtMCHRec, "mPtMCHRec/F");
    mTreeTrackClusterSizeData -> Branch("mEtaMCHRec", &mEtaMCHRec, "mEtaMCHRec/F");
    mTreeTrackClusterSizeData -> Branch("mPhiMCHRec", &mPhiMCHRec, "mPhiMCHRec/F");

    mTreeTrackClusterSizeData -> Branch("mChargeMFTRec", &mChargeMFTRec, "mChargeMFTRec/I");
    mTreeTrackClusterSizeData -> Branch("mPtMFTRec", &mPtMFTRec, "mPtMFTRec/F");
    mTreeTrackClusterSizeData -> Branch("mEtaMFTRec", &mEtaMFTRec, "mEtaMFTRec/F");
    mTreeTrackClusterSizeData -> Branch("mPhiMFTRec", &mPhiMFTRec, "mPhiMFTRec/F");

    mTreeTrackClusterSizeData -> Branch("mMch2MCHMIDRec", &mMch2MCHMIDRec, "mMch2MCHMIDRec/F");
    mTreeTrackClusterSizeData -> Branch("mMch2MFTMCHRec", &mMch2MFTMCHRec, "mMch2MFTMCHRec/F");

    // Clsize per layer
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer0", &mClsizeRecLayer0);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer1", &mClsizeRecLayer1);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer2", &mClsizeRecLayer2);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer3", &mClsizeRecLayer3);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer4", &mClsizeRecLayer4);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer5", &mClsizeRecLayer5);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer6", &mClsizeRecLayer6);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer7", &mClsizeRecLayer7);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer8", &mClsizeRecLayer8);
    mTreeTrackClusterSizeData -> Branch("mClsizeRecLayer9", &mClsizeRecLayer9);
    // Pattern ID per layer
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer0", &mPattIdRecLayer0);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer1", &mPattIdRecLayer1);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer2", &mPattIdRecLayer2);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer3", &mPattIdRecLayer3);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer4", &mPattIdRecLayer4);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer5", &mPattIdRecLayer5);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer6", &mPattIdRecLayer6);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer7", &mPattIdRecLayer7);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer8", &mPattIdRecLayer8);
    mTreeTrackClusterSizeData -> Branch("mPattIdRecLayer9", &mPattIdRecLayer9);
    // Mean cluster size
    mTreeTrackClusterSizeData -> Branch("mMeanClsizePerTrackRec", &mMeanClsizePerTrackRec);

    // get MFT tracks
    std::vector<o2::mft::TrackMFT> *mMFTTracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTTracksROF = nullptr;
    std::vector<int> *mMFTTrackClusIdx = nullptr;

    // get MFT clusters
    std::vector<o2::itsmft::CompClusterExt> *mMFTClusters = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTClustersROF = nullptr;
    std::vector<unsigned char> *mMFTClusterPatterns = nullptr;

    std::vector<int> ClTrackID;
    std::vector<int> ClEvID;
    std::vector<int> ClSiezes;

    treeMFTTracks -> SetBranchAddress("MFTTrack", &mMFTTracks);
    treeMFTTracks -> SetBranchAddress("MFTTracksROF", &mMFTTracksROF);
    treeMFTTracks -> SetBranchAddress("MFTTrackClusIdx", &mMFTTrackClusIdx);
    
    treeMFTClusters -> SetBranchAddress("MFTClusterComp", &mMFTClusters);
    treeMFTClusters -> SetBranchAddress("MFTClustersROF", &mMFTClustersROF);
    treeMFTClusters -> SetBranchAddress("MFTClusterPatt", &mMFTClusterPatterns);
    
    std::cout << treeMFTClusters -> GetEntriesFast() << std::endl;
    for (int frame = 0; frame < treeMFTClusters -> GetEntriesFast(); frame++) {
        if (!treeMFTClusters -> GetEvent(frame) || !treeMFTTracks -> GetEvent(frame)) {
            LOG(info) << "Skipping frame: " << frame;
            continue;
        }

        std::vector<o2::itsmft::ClusterPattern> pattVec;
        pattVec.reserve(mMFTClusters -> size());
        auto pattItMFT = mMFTClusterPatterns -> begin();

        for (unsigned int iClus{0}; iClus < mMFTClusters -> size(); ++iClus) {
            auto &clus = mMFTClusters -> at(iClus);
            auto pattID = clus.getPatternID();
            int npix;
            o2::itsmft::ClusterPattern patt;
            if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mDictionary.isGroup(pattID)) {
                patt.acquirePattern(pattItMFT);
                npix = patt.getNPixels();
            } else {
                npix = mDictionary.getNpixels(pattID);
                patt = mDictionary.getPattern(pattID);
            }
            mTrackClusterSize -> Fill(npix);
            pattVec.push_back(patt);
        }

        for (unsigned int iTrack{0}; iTrack < mMFTTracks -> size(); ++iTrack) {
            auto &oneTrack = mMFTTracks -> at(iTrack);
            auto ncls = oneTrack.getNumberOfPoints();
            auto offset = oneTrack.getExternalClusterIndexOffset();

            float mean = 0, norm = 0;

            mChargeMFTRec = oneTrack.getCharge();
            mPtMFTRec = oneTrack.getPt();
            mEtaMFTRec = oneTrack.getEta();
            mPhiMFTRec = oneTrack.getPhi();

            for (int icls = 0; icls < ncls; ++icls) {
                auto clsEntry = mMFTTrackClusIdx -> at(offset + icls);
                auto &oneCluster = mMFTClusters -> at(clsEntry);
                int clsLayer = mMFTChipMapper.chip2Layer(oneCluster.getChipID());
                auto &patt = pattVec.at(clsEntry);
                auto pattID = oneCluster.getPatternID();
                int npix = patt.getNPixels();
                //std::cout << npix << " + ";

                switch(clsLayer) {
                    case 0:
                        mClsizeRecLayer0 = npix;
                        mPattIdRecLayer0 = pattID;
                        break;
                    case 1:
                        mClsizeRecLayer1 = npix;
                        mPattIdRecLayer1 = pattID;
                        break;
                    case 2:
                        mClsizeRecLayer2 = npix;
                        mPattIdRecLayer2 = pattID;
                        break;
                    case 3:
                        mClsizeRecLayer3 = npix;
                        mPattIdRecLayer3 = pattID;
                        break;
                    case 4:
                        mClsizeRecLayer4 = npix;
                        mPattIdRecLayer4 = pattID;
                        break;
                    case 5:
                        mClsizeRecLayer5 = npix;
                        mPattIdRecLayer5 = pattID;
                        break;
                    case 6:
                        mClsizeRecLayer6 = npix;
                        mPattIdRecLayer6 = pattID;
                        break;
                    case 7:
                        mClsizeRecLayer7 = npix;
                        mPattIdRecLayer7 = pattID;
                        break;
                    case 8:
                        mClsizeRecLayer8 = npix;
                        mPattIdRecLayer8 = pattID;
                        break;
                    case 9:
                        mClsizeRecLayer9 = npix;
                        mPattIdRecLayer9 = pattID;
                        break;
                }

                mean += npix;
                norm += 1;
            }
            mean /= norm;
            mMeanClsizePerTrackRec = mean;
            
            //std::cout << " -> MEAN CLs = " << mMFTTracks -> size() << std::endl;
            mTrackMeanClusterSize -> Fill(mean);
            mTrackEta -> Fill(oneTrack.getEta());

            // Fill the tree
            //std::cout << chargeMFT << " " << ptMFT << " " << etaMFT << " " << phiMFT << std::endl;
            mTreeTrackClusterSizeData -> Fill();

            mChargeGlobalRec = -999;
            mPtGlobalRec = -999;
            mEtaGlobalRec = -999;
            mPhiGlobalRec = -999;
            
            mChargeMCHRec = -999;
            mPtMCHRec = -999;
            mEtaMCHRec = -999;
            mPhiMCHRec = -999;

            mChargeMFTRec = -999;
            mPtMFTRec = -999;
            mEtaMFTRec = -999;
            mPhiMFTRec = -999;
            // Clsize per layer
            mClsizeRecLayer0 = -999;
            mClsizeRecLayer1 = -999;
            mClsizeRecLayer2 = -999;
            mClsizeRecLayer3 = -999;
            mClsizeRecLayer4 = -999;
            mClsizeRecLayer5 = -999;
            mClsizeRecLayer6 = -999;
            mClsizeRecLayer7 = -999;
            mClsizeRecLayer8 = -999;
            mClsizeRecLayer9 = -999;
            // Pattern ID per layer
            mPattIdRecLayer0 = -999;
            mPattIdRecLayer1 = -999;
            mPattIdRecLayer2 = -999;
            mPattIdRecLayer3 = -999;
            mPattIdRecLayer4 = -999;
            mPattIdRecLayer5 = -999;
            mPattIdRecLayer6 = -999;
            mPattIdRecLayer7 = -999;
            mPattIdRecLayer8 = -999;
            mPattIdRecLayer9 = -999;
        }
    }

    fOut -> cd();
    mTrackClusterSize -> Write();
    mTrackMeanClusterSize -> Write();
    mTrackEta -> Write();

    mGlobalTrackMeanClusterSize -> Write();
    mGlobalTrackPt -> Write();
    mGlobalTrackEta -> Write();
    mGlobalTrackEtaMFT -> Write();
    mGlobalTrackEtaMCH -> Write();
    mGlobalTrackEtaMFTMCHMatchingChi2 -> Write();

    mGlobalTrackMeanClusterSizeMuon -> Write();
    mGlobalTrackPtMuon -> Write();
    mGlobalTrackEtaMuon -> Write();
    mGlobalTrackEtaMFTMuon -> Write();
    mGlobalTrackEtaMCHMuon -> Write();
    mGlobalTrackEtaMFTMCHMatchingChi2Muon -> Write();

    mMIDMatchingChi2 -> Write();
    mMFTMCHMatchingChi2 -> Write();

    mTreeTrackClusterSizeData -> Write();
    fOut -> Close();
}