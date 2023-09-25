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


void cluster_studies_mc(){

    //string pathIn = "cross_check/mc";
    string pathIn = "/Users/lucamicheletti/cernbox/summer_student_mft_pid/mc/24_07_23_prod";

    TH1F *mTrackClusterSize = new TH1F("mMFTTrackClusterSize", "Cluster Size Per Track; # clusters; # entries", 100, 0.5, 100.5);
    TH1F *mTrackClusterSizeMC = new TH1F("mMFTTrackClusterSizeMC", "Cluster Size Per Track MC; # clusters; # entries", 100, 0.5, 100.5);
    TH1F *mTrackMeanClusterSize = new TH1F("mMFTTrackMeanClusterSize", "Mean Cluster Size Per Track; # clusters; # entries", 100, 0.5, 20.5);
    TH1F *mTrackMeanClusterSize3He = new TH1F("mTrackMeanClusterSize3He", "Mean Cluster Size Per Track 3He; # clusters; # entries", 100, 0.5, 20.5);
    TH1F *mTrackMeanClusterSizeMu = new TH1F("mTrackMeanClusterSizeMu", "Mean Cluster Size Per Track Mu; # clusters; # entries", 100, 0.5, 20.5);
    TH1F *mTrackMeanClusterSizePi = new TH1F("mTrackMeanClusterSizePi", "Mean Cluster Size Per Track Pi; # clusters; # entries", 100, 0.5, 20.5);
    TH1F *mTrackMeanClusterSizeKaon = new TH1F("mTrackMeanClusterSizeKaon", "Mean Cluster Size Per Track Kaon; # clusters; # entries", 100, 0.5, 20.5);
    TH1F *mTrackMeanClusterSizeProton = new TH1F("mTrackMeanClusterSizeProton", "Mean Cluster Size Per Track Proton; # clusters; # entries", 100, 0.5, 20.5);

    auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_1320237.root"); // March 2022
    //auto fIn = TFile("ccdb/o2_itsmft_TopologyDictionary_0747987.root"); // December 2022
    o2::itsmft::TopologyDictionary mDictionary = *(reinterpret_cast<o2::itsmft::TopologyDictionary *>(fIn.Get("ccdb_object")));

    printf("Running a first test\n");
    TFile *fMFTClusters = new TFile(Form("%s/mftclusters.root", pathIn.c_str()), "READ");
    fMFTClusters -> ls();
    TFile *fMFTTracks = new TFile(Form("%s/mfttracks.root", pathIn.c_str()), "READ");
    fMFTTracks -> ls();
    TFile *fMCTracks = new TFile(Form("%s/sgn_1_Kine.root", pathIn.c_str()), "READ");
    fMCTracks -> ls();

    
    TTree *treeMFTClusters = (TTree*) fMFTClusters -> Get("o2sim");
    TTree *treeMFTTracks = (TTree*) fMFTTracks -> Get("o2sim");
    TTree *treeMCTracks = treeMCTracks = (TTree*) fMCTracks -> Get("o2sim");
    //TTree *treeMCTracks = nullptr;

    // get tracks
    std::vector<o2::mft::TrackMFT> *mMFTTracks = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTTracksROF = nullptr;
    std::vector<int> *mMFTTrackClusIdx = nullptr;

    // get clusters
    std::vector<o2::itsmft::CompClusterExt> *mMFTClusters = nullptr;
    std::vector<o2::itsmft::ROFRecord> *mMFTClustersROF = nullptr;
    std::vector<unsigned char> *mMFTClusterPatterns = nullptr;

    // get MC tracks
    std::vector<o2::MCTrack> *mMCtracks = nullptr;
    std::vector<std::vector<std::vector<int>>> mTrueTracksMap;
    o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mMFTClusterLabels = nullptr;
    vector<o2::MCCompLabel> *mMFTTrackLabels = nullptr;

    //std::vector<const unsigned char> *pattIt = nullptr;
    //std::vector<o2::BaseCluster<float>> *mMFTClustersGlobal;

    std::vector<int> ClTrackID;
    std::vector<int> ClEvID;
    std::vector<int> ClSiezes;

    treeMFTTracks -> SetBranchAddress("MFTTrack", &mMFTTracks);
    treeMFTTracks -> SetBranchAddress("MFTTracksROF", &mMFTTracksROF);
    treeMFTTracks -> SetBranchAddress("MFTTrackClusIdx", &mMFTTrackClusIdx);
    treeMFTTracks -> SetBranchAddress("MFTTrackMCTruth", &mMFTTrackLabels);

    treeMFTClusters -> SetBranchAddress("MFTClusterComp", &mMFTClusters);
    treeMFTClusters -> SetBranchAddress("MFTClustersROF", &mMFTClustersROF);
    treeMFTClusters -> SetBranchAddress("MFTClusterPatt", &mMFTClusterPatterns);
    treeMFTClusters -> SetBranchAddress("MFTClusterMCTruth", &mMFTClusterLabels);

    treeMCTracks->SetBranchAddress("MCTrack", &mMCtracks);
    
    // Alternative way to read the kinematics file -> not finalized
    //o2::steer::MCKinematicsReader mMCReader(Form("%s/sgn_2", pathIn.c_str()), o2::steer::MCKinematicsReader::Mode::kMCKine);
    
    int myCounter = 0;
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
        cout << "N clusters = " << mMFTClusters -> size() << std::endl;
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
            mTrackClusterSize -> Fill(npix);
            pattVec.push_back(patt);

            ///
            auto &labCls = (mMFTClusterLabels -> getLabels(iClus))[0];
            int  trackID, evID, srcID;
            bool fake;
            labCls.get(trackID, evID, srcID, fake);
            //if (!labCls.isNoise() && labCls.isValid() && labCls.isCorrect() && !labCls.isFake()) {
                ClTrackID.push_back(trackID);
                ClEvID.push_back(evID);
                ClSiezes.push_back(npix);
            //}
            ///
        }

        //std::cout << "ENTERING THE CODE FOR MC ASSOCIATION" << std::endl;
        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
        auto nev = treeMCTracks -> GetEntriesFast();
        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++) {
            treeMCTracks -> GetEvent(n);
            mcTracksMatrix[n].resize(mMCtracks -> size());
            for (unsigned int mcI{0}; mcI < mMCtracks -> size(); ++mcI) {
                mcTracksMatrix[n][mcI] = mMCtracks -> at(mcI);
            }
        }

        //std::cout << mMFTTracks -> size() << std::endl;
        for (unsigned int iTrack{0}; iTrack < mMFTTracks -> size(); ++iTrack) {
            if (myCounter>100) {
				break;
			}
            auto &oneTrack = mMFTTracks -> at(iTrack);
            auto ncls = oneTrack.getNumberOfPoints();
            auto offset = oneTrack.getExternalClusterIndexOffset();

            float mean = 0, norm = 0;
            int trackPDG = -999999999;

            //auto &trackLabel = (mMFTTrackLabels -> getLabels(iTrack))[0];
            auto &trackLabel = mMFTTrackLabels -> at(iTrack);
            if (trackLabel.isCorrect()) {
                int  trackID, evID, srcID;
                bool fake;
                trackLabel.get(trackID, evID, srcID, fake);
                //LOG(info) << "(Labels info: trackID="<<trackID<<", eventID="<<evID<<", srcID="<<srcID;
                trackPDG = mcTracksMatrix[evID][trackID].GetPdgCode();
            } else {
                continue;
            }

            for (int icls = 0; icls < ncls; ++icls) {
                auto clsEntry = mMFTTrackClusIdx -> at(offset + icls);
                auto &oneCluster = mMFTClusters -> at(clsEntry);
                auto &patt = pattVec.at(clsEntry);
                //auto globalCluster = mMFTClustersGlobal[clsEntry];
                int npix = patt.getNPixels();
                std::cout << npix << " (" << clsEntry << ") + ";

                mean += npix;
                norm += 1;
            }
            mean /= norm;
            std::cout << " -> " << mean << " - " << trackPDG << std::endl;
            mTrackMeanClusterSize -> Fill(mean);
            if (TMath::Abs(trackPDG) == 1000020030) {
                mTrackMeanClusterSize3He -> Fill(mean);
            }

            if (TMath::Abs(trackPDG) == 13) {
                mTrackMeanClusterSizeMu -> Fill(mean);
            }

            if (TMath::Abs(trackPDG) == 211) {
                mTrackMeanClusterSizePi -> Fill(mean);
            }

            if (TMath::Abs(trackPDG) == 321) {
                mTrackMeanClusterSizeKaon -> Fill(mean);
            }

            if (TMath::Abs(trackPDG) == 2212) {
                mTrackMeanClusterSizeProton -> Fill(mean);
            }
            //std::cout << std::endl;
            myCounter++;
        }

        TFile *fOut = new TFile(Form("%s/cross_check_mc_test.root", pathIn.c_str()), "RECREATE");
        mTrackClusterSize -> Write();
        mTrackMeanClusterSize -> Write();
        mTrackClusterSizeMC -> Write();
        mTrackMeanClusterSize3He -> Write();
        mTrackMeanClusterSizeMu -> Write();
        mTrackMeanClusterSizePi -> Write();
        mTrackMeanClusterSizeKaon -> Write();
        mTrackMeanClusterSizeProton -> Write();
        fOut -> Close();

        return;
        //if (isMC) {
            /*int myCounter = 0;
            for (auto src = 0; src < mMCReader.getNSources(); src++) {
                for (Int_t event = 0; event < mMCReader.getNEvents(src); event++) {
                    auto evH = mMCReader.getMCEventHeader(src, event);

                    for (const auto& trueMFTTrackID : mTrueTracksMap[src][event]) {
                        auto mftTrack = mMFTTracks -> at(trueMFTTrackID);
                    }
                }
            }*/

            /*std::cout << "ENTERING THE CODE FOR MC ASSOCIATION" << std::endl;
            std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
            auto nev = treeMCTracks -> GetEntriesFast();
            mcTracksMatrix.resize(nev);
            for (int n = 0; n < nev; n++) {
                treeMCTracks -> GetEvent(n);
                mcTracksMatrix[n].resize(mMCtracks -> size());
                for (unsigned int mcI{0}; mcI < mMCtracks -> size(); ++mcI) {
                    mcTracksMatrix[n][mcI] = mMCtracks -> at(mcI);
                }
            }

            LOG(info) << "---- GETTING MC tracks info ----";
            for (int i = 0; i < ClEvID.size(); i++) {
                auto evID = ClEvID.at(i);
                auto trID = ClTrackID.at(i);
                auto trackPDG = mcTracksMatrix[evID][trID].GetPdgCode();
                mTrackClusterSizeMC -> Fill(ClSiezes.at(i));
            }
        }*/

    }
}