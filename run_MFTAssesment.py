import os
import sys
import argparse
import yaml
import random

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

    #if not os.path.isdir(inputCfg["output"]["output_dir"]) :
        #print("the directory does not exist, creating %s" % (inputCfg["output"]["output_dir"]))
        #os.system("mkdir -p %s" % (inputCfg["output"]["output_dir"]))

    #if mode == "test" :
        #fIn  = open(inputCfg["input"]["run_list_file"], "r")

    #if mode == "full" :
        #fIn  = open(inputCfg["input"]["run_list_file"], "r")
        #fOut = open("%s/%s" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"]), "w")

    #if mode == "terminate" :
        #if not os.path.isfile("%s/%s" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"])) :
            #print('Submitted jobs file does not exist! Do --full before')
            #return
        #fIn  = open("%s/%s" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"]), "r")
        #fOut = open("%s/%s" % (inputCfg["output"]["output_dir"], inputCfg["output"]["terminated_output_file"]), "w")

    #for run in fIn:
        #print(fr"aliroot -b -q %s\(\"%s\",%i\)" % (inputCfg["input"]["macro_to_run"], mode, int(run)))
        #if mode == "test" :
            #break
        #fOut.write(run)

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--run", help="test your task", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.run:
        run_mft_assesment(inputCfg, "full")

main()
