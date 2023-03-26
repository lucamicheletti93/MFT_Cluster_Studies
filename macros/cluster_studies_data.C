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

    string pathIn = "cross_check/data";

    TH1F *mTrackClusterSize = new TH1F("mMFTTrackClusterSize", "Cluster Size Per Track; # clusters; # entries", 100, 0.5, 100.5);
    TH1F *mTrackMeanClusterSize = new TH1F("mMFTTrackMeanClusterSize", "Mean Cluster Size Per Track; # clusters; # entries", 200, 0.5, 10.5);

    auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_1320237.root"); // March 2022
    //auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_0747987.root"); // December 2022
    o2::itsmft::TopologyDictionary mDictionary = *(reinterpret_cast<o2::itsmft::TopologyDictionary *>(fIn.Get("ccdb_object")));

    printf("Running a first test\n");
    TFile *fMFTClusters = new TFile(Form("%s/mftclusters.root", pathIn.c_str()), "READ");
    fMFTClusters -> ls();
    TFile *fMFTTracks = new TFile(Form("%s/mfttracks.root", pathIn.c_str()), "READ");
    fMFTTracks -> ls();

    
    TTree *treeMFTClusters = (TTree*) fMFTClusters -> Get("o2sim");
    TTree *treeMFTTracks = (TTree*) fMFTTracks -> Get("o2sim");

    // get tracks
    std::vector<o2::mft::TrackMFT> *mMFTTracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTTracksROF = nullptr;
    std::vector<int> *mMFTTrackClusIdx = nullptr;

    // get clusters
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
        cout << "N clusters = " << mMFTClusters -> size() << std::endl;

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
                std::cout << npix << " + ";

                mean += npix;
                norm += 1;
            }
            mean /= norm;
            std::cout << " -> MEAN CLs = " << mean << std::endl;
            mTrackMeanClusterSize -> Fill(mean);
        }

        TFile *fOut = new TFile(Form("%s/cross_check_data.root", pathIn.c_str()), "RECREATE");
        mTrackClusterSize -> Write();
        mTrackMeanClusterSize -> Write();
        fOut -> Close();

    }
}