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

    std::vector<o2::BaseCluster<float>> mMFTClustersGlobal;

    treeMFTClusters -> SetBranchAddress("MFTTrack", &mMFTTracks);
    treeMFTClusters -> SetBranchAddress("MFTTracksROF", &mMFTTracksROF);
    treeMFTClusters -> SetBranchAddress("MFTTrackClusIdx", &mMFTTrackClusIdx);

    treeMFTTracks -> SetBranchAddress("MFTClusterComp", &mMFTClusters);
    treeMFTTracks -> SetBranchAddress("MFTClustersROF", &mMFTClustersROF);
    treeMFTTracks -> SetBranchAddress("MFTClusterPatt", &mMFTClusterPatterns);

    for (int frame = 0; frame < treeMFTClusters -> GetEntriesFast(); frame++) {
        // LOOP OVER FRAMES  
        //LOG(info) << "Skipping frame: " << frame;
        //o2::mft::ioutils::convertCompactClusters(mMFTClusters, pattIt, mMFTClustersGlobal, mDictionary);
        // LOOP OVER FRAMES  
        if (!treeMFTClusters -> GetEvent(frame) || !treeMFTTracks -> GetEvent(frame)) {
            LOG(info) << "Skipping frame: " << frame;
            continue;
        }
        for (unsigned int iClus{0}; iClus < ITSclus->size(); iClus++) {

        }
    }

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