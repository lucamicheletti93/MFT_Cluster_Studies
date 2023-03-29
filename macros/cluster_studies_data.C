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
#include "DetectorsVertexing/DCAFitterN.h"

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
            double ch2MCHMID = oneGlobalFwdTrack.getMIDMatchingChi2();
            double ch2MFTMCH = oneGlobalFwdTrack.getMFTMCHMatchingChi2();
            double pt = oneGlobalFwdTrack.getPt();
            double eta = oneGlobalFwdTrack.getEta();
            double etaMFT = oneMftTrack.getEta();
            
            double pMCH = oneMchTrack.getP();
            double pzMCH = oneMchTrack.getPz();
            double etaMCH = 0.5 * TMath::Log((pMCH + pzMCH) / (pMCH - pzMCH));

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

            mMIDMatchingChi2 -> Fill(ch2MCHMID);
            mMFTMCHMatchingChi2 -> Fill(ch2MFTMCH);
            
            mGlobalTrackMeanClusterSize -> Fill(mean);
            mGlobalTrackPt -> Fill(pt);
            mGlobalTrackEta -> Fill(eta);
            mGlobalTrackEtaMFT -> Fill(etaMFT);
            mGlobalTrackEtaMCH -> Fill(etaMCH);
            mGlobalTrackEtaMFTMCHMatchingChi2 -> Fill(eta, ch2MFTMCH);

            if (ch2MCHMID > 0 && TMath::Abs(etaMCH) > 2.5 && TMath::Abs(etaMCH) < 4) {
                mGlobalTrackMeanClusterSizeMuon -> Fill(mean);
                mGlobalTrackPtMuon -> Fill(pt);
                mGlobalTrackEtaMuon -> Fill(eta);
                mGlobalTrackEtaMFTMuon -> Fill(etaMFT);
                mGlobalTrackEtaMCHMuon -> Fill(etaMCH);
                mGlobalTrackEtaMFTMCHMatchingChi2Muon -> Fill(eta, ch2MFTMCH);
            }
        }

        TFile *fOut = new TFile(Form("%s/cross_check_data.root", pathIn.c_str()), "RECREATE");
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
        fOut -> Close();

    }
}