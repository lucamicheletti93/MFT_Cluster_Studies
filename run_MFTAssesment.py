import os
import sys
import argparse
import yaml
import random
import ROOT

def run_mft_assesment(inputCfg, mode):
    '''
    function to run the MFTAssesment workflow on time frames
    '''
    paths = []
    for file in os.listdir(inputCfg["input"]["input_path"]):
        d = os.path.join(inputCfg["input"]["input_path"], file)
        if os.path.isdir(d):
            paths.append(d)
            #print(d)

    counter = 0
    for path in paths:
        if mode == "test" and counter > 0:
            exit()
        print("----------------------   %s   ----------------------" % (path))
        os.chdir("%s" % (path))
        os.system("o2-mft-reco-workflow | o2-mft-assessment-workflow")
        counter = counter + 1
        print("----------------------> %s processed" % (path))

def merge_files(inputCfg):
    '''
    function for merging files
    '''
    paths = []
    counter = 0
    for file in os.listdir(inputCfg["input"]["input_path"]):
        d = os.path.join(inputCfg["input"]["input_path"], file)
        if os.path.isdir(d) and os.path.isfile("%s/MFTAssessment.root" % (d)):
            if counter == 0:
                mergeCommand = "hadd -f MFTAssessment.root "
            mergeCommand += "%s/MFTAssessment.root " % (d)
            counter = counter + 1
    
    print(mergeCommand)
    os.system(mergeCommand)

def plots():
    '''
    function for producing QC plots
    '''
    fIn = ROOT.TFile.Open("MFTAssessment.root", "READ")
    mMFTTrueTrackMeanClusterSize = fIn.Get("mMFTTrueTrackMeanClusterSize")
    mMFTTrueTrackMeanClusterSize.SetLineColor(ROOT.kBlack)
    mMFTTrueTrackMeanClusterSize.Rebin(10)
    mMFTTrueTrackMeanClusterSizeNuclei = fIn.Get("mMFTTrueTrackMeanClusterSizeNuclei")
    mMFTTrueTrackMeanClusterSizeNuclei.SetLineColor(ROOT.kRed)
    mMFTTrueTrackMeanClusterSizeNuclei.Rebin(10)
    mMFTTrueTrackMeanClusterSizeNonNuclei = fIn.Get("mMFTTrueTrackMeanClusterSizeNonNuclei")
    mMFTTrueTrackMeanClusterSizeNonNuclei.SetLineColor(ROOT.kBlue)
    mMFTTrueTrackMeanClusterSizeNonNuclei.Rebin(10)

    canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
    mMFTTrueTrackMeanClusterSize.Draw("H")
    mMFTTrueTrackMeanClusterSizeNuclei.Draw("Hsame")
    mMFTTrueTrackMeanClusterSizeNonNuclei.Draw("Hsame")
    canvas.Update()

    fOut = ROOT.TFile.Open("MyAnalysisResults.root", "RECREATE")
    canvas.Write()
    fOut.Close()

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--run_full", help="run over all the available files", action="store_true")
    parser.add_argument("--merge", help="merge all the files produced", action="store_true")
    parser.add_argument("--plot", help="plot results", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.run_full:
        run_mft_assesment(inputCfg, "full")
    if args.merge:
        merge_files(inputCfg)
    if args.plot:
        plots()

main()
