#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>

#include "CommonDataFormat/RangeReference.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TrkClusRef.h"
#include "DataFormatsMFT/TrackMFT.h"

#include "ITStracking/IOUtils.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

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

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CCDBTimeStampUtils.h"

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


void cluster_studies(){
    TH1F *histClsize = new TH1F("histClsize", "", 100, -0.5, 99.5);

    //LOG(info) << "Loading LATEST dictionary: if you are analysing data older than JUNE check out the dictionary";
    auto fIn = TFile("o2_itsmft_TopologyDictionary_1653153873993.root");
    o2::itsmft::TopologyDictionary mDictionary = *(reinterpret_cast<o2::itsmft::TopologyDictionary *>(fIn.Get("ccdb_object")));

    printf("Running a first test\n");
    TFile *fMFTClusters = new TFile("mfttracks.root", "READ");
    fMFTClusters -> ls();
    TFile *fMFTTracks = new TFile("mftclusters.root", "READ");
    fMFTTracks -> ls();
    //auto fMCTracks = TFile::Open("sgn_1_Kine.root");

    TTree *treeMFTClusters = (TTree*) fMFTClusters -> Get("o2sim");
    TTree *treeMFTTracks = (TTree*) fMFTTracks -> Get("o2sim");
    //TTree *treeMCTracks = (TTree*) fMCTracks -> Get("o2sim");

    // get tracks
    std::vector<o2::mft::TrackMFT> *mMFTTracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTTracksROF = nullptr;
    std::vector<int> *mMFTTrackClusIdx = nullptr;

    // get clusters
    std::vector<o2::itsmft::CompClusterExt> *mMFTClusters = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTClustersROF = nullptr;
    std::vector<unsigned char> *mMFTClusterPatterns = nullptr;

    //std::vector<const unsigned char> *pattIt = nullptr;
    //std::vector<o2::BaseCluster<float>> *mMFTClustersGlobal;

    treeMFTClusters -> SetBranchAddress("MFTTrack", &mMFTTracks);
    treeMFTClusters -> SetBranchAddress("MFTTracksROF", &mMFTTracksROF);
    treeMFTClusters -> SetBranchAddress("MFTTrackClusIdx", &mMFTTrackClusIdx);

    treeMFTTracks -> SetBranchAddress("MFTClusterComp", &mMFTClusters);
    treeMFTTracks -> SetBranchAddress("MFTClustersROF", &mMFTClustersROF);
    treeMFTTracks -> SetBranchAddress("MFTClusterPatt", &mMFTClusterPatterns);

    std::cout << treeMFTClusters -> GetEntriesFast() << std::endl;
    for (int frame = 0; frame < treeMFTClusters -> GetEntriesFast(); frame++) {
        // LOOP OVER FRAMES  
        //LOG(info) << "Skipping frame: " << frame;
        //o2::mft::ioutils::convertCompactClusters(mMFTClusters, pattIt, mMFTClustersGlobal, mDictionary);
        // LOOP OVER FRAMES  
        if (!treeMFTClusters -> GetEvent(frame) || !treeMFTTracks -> GetEvent(frame)) {
            LOG(info) << "Skipping frame: " << frame;
            continue;
        }
        cout << mMFTClusters -> size() << std::endl;
        //pattIt = mMFTClusterPatterns -> begin();
        //mMFTClustersGlobal -> clear();
        //mMFTClustersGlobal -> reserve(mMFTClusters -> size());
        //o2::mft::ioutils::convertCompactClusters(mMFTClusters, pattIt, mMFTClustersGlobal, mDictionary);

        std::vector<o2::itsmft::ClusterPattern> pattVec;
        pattVec.reserve(mMFTClusters -> size());
        auto pattItMFT = mMFTClusterPatterns -> begin();

        for (unsigned int iClus{0}; iClus < mMFTClusters -> size(); ++iClus) {
            //auto &clus = mMFTClusters[iClus];
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
            histClsize -> Fill(npix);
            pattVec.push_back(patt);
        }

        std::cout << mMFTTracks -> size() << std::endl;
        for (unsigned int iTrack{0}; iTrack < mMFTTracks -> size(); ++iTrack) {
            auto &oneTrack = mMFTTracks -> at(iTrack);
            auto ncls = oneTrack.getNumberOfPoints();
            auto offset = oneTrack.getExternalClusterIndexOffset();

            for (int icls = 0; icls < ncls; ++icls) {
                auto clsEntry = mMFTTrackClusIdx -> at(offset + icls);
                auto &oneCluster = mMFTClusters -> at(clsEntry);
                auto &patt = pattVec.at(clsEntry);
                //auto globalCluster = mMFTClustersGlobal[clsEntry];
                int npix = patt.getNPixels();
                std::cout << npix << " , ";
            }
            std::cout << std::endl;
        }

    }

    TFile *fOut = new TFile("test_cluster_studies.root", "RECREATE");
    histClsize -> Write();
    fOut -> Close();

    /*
    mMFTTracks = ctx.inputs().get<gsl::span<o2::mft::TrackMFT>>("tracks");
    mMFTTracksROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("tracksrofs");
    mMFTTrackClusIdx = ctx.inputs().get<gsl::span<int>>("trackClIdx");

    
    mMFTClusters = ctx.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compClusters");
    mMFTClustersROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("clustersrofs");
    mMFTClusterPatterns = ctx.inputs().get<gsl::span<unsigned char>>("patterns");
    */

    /*
    std::vector<o2::itsmft::ClusterPattern> pattVec;
    pattVec.reserve(mMFTClusters -> size());
    auto pattItMFT = mMFTClusterPatterns -> begin();
    for (unsigned int iClus{0}; iClus < mMFTClusters -> size(); ++iClus)
    {
        std::cout << "ciao" << std::endl;
    }
    */
}