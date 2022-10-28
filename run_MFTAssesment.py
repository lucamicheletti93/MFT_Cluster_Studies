import os
import sys
import argparse
import yaml
import random
import ROOT
from ROOT import gROOT, gBenchmark, gPad, gStyle, kTRUE, kFALSE
from utils.plot_library import LoadStyle, SetLegend

def copy_from_grid(inputCfg):
    '''
    function for downloading files from alien grid (Run mft reco. on ctf)
    '''
    with open(inputCfg["input"]["run_list_file"]) as f:
        lines = f.readlines()
    
    for i in range(0, len(lines)):
        if not os.path.isdir("{}{}".format(inputCfg["input"]["ctf_input_path"], i)):
            os.system("mkdir -p {}{}".format(inputCfg["input"]["ctf_input_path"], i))
        os.system("alien_cp alien://%s/%s file:%s%i" % (inputCfg["input"]["alien_input_path"], lines[i].replace("\n", ""), inputCfg["input"]["ctf_input_path"], i))

    for i in range(0, len(lines)):
        os.chdir("%s%i" % (inputCfg["input"]["ctf_input_path"], i))
        if not os.path.isfile("mftclusters.root") or not os.path.isfile("mfttracks.root"):
            os.system("o2-ctf-reader-workflow --onlyDet MFT  --ctf-input %s | o2-mft-reco-workflow --clusters-from-upstream --disable-mc  -b" % (lines[i].replace("\n", "")))
        if not os.path.isfile("o2simdigitizerworkflow_configuration.ini"):
            os.system("cp /home/lmichele/alice/mft_data/CTF/tf0/o2simdigitizerworkflow_configuration.ini .")
        if not os.path.isfile("mftdigits.roots"):
            os.system("cp /home/lmichele/alice/mft_data/CTF/tf0/mftdigits.root .")


def run_mft_assesment(inputCfg, mode):
    '''
    function to run the MFTAssesment workflow on time frames
    '''
    paths = []
    for file in os.listdir(inputCfg["input"]["local_input_path"]):
        d = os.path.join(inputCfg["input"]["local_input_path"], file)
        if os.path.isdir(d):
            paths.append(d)
            print(d)

    counter = 0
    for path in paths:
        if mode == "test" and counter > 0:
            exit()
        print("----------------------   %s   ----------------------" % (path))
        os.chdir("%s" % (path))
        if inputCfg["input"]["prod_type"] == "data":
            os.system("o2-mft-reco-workflow | o2-mft-assessment-workflow --disable-mc")
        if inputCfg["input"]["prod_type"] == "mc":
            os.system("o2-mft-reco-workflow | o2-mft-assessment-workflow")
        counter = counter + 1
        print("----------------------> %s processed" % (path))

def merge_files(inputCfg):
    '''
    function for merging files
    '''
    paths = []
    counter = 0
    for file in os.listdir(inputCfg["input"]["local_input_path"]):
        d = os.path.join(inputCfg["input"]["local_input_path"], file)
        if os.path.isdir(d) and os.path.isfile("%s/MFTAssessment.root" % (d)):
            if counter == 0:
                mergeCommand = "hadd -f {} ".format(inputCfg["output"]["local_output_file"])
            mergeCommand += "%s/MFTAssessment.root " % (d)
            counter = counter + 1
    
    print(mergeCommand)
    os.system(mergeCommand)

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--download", help="download files from grid", action="store_true")
    parser.add_argument("--run_full", help="run over all the available files", action="store_true")
    parser.add_argument("--run_test", help="run over a part of the sample", action="store_true")
    parser.add_argument("--merge", help="merge all the files produced", action="store_true")
    args = parser.parse_args()


    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.download:
        copy_from_grid(inputCfg)
    if args.run_full:
        run_mft_assesment(inputCfg, "full")
    if args.run_test:
        run_mft_assesment(inputCfg, "test")
    if args.merge:
        merge_files(inputCfg)

main()
