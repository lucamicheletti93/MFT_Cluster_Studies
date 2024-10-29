import argparse
import numpy as np
import pandas as pd
import uproot
import awkward as ak
from alive_progress import alive_bar


def read_tree_fast(input_path, output_path, suffix=""):
    # Define output ROOT file path and .parquet path
    root_output_file = f"{output_path}/treeMftPidTraining{suffix}.root"
    parquet_output_file = f"{output_path}/treeMftPidTraining{suffix}.gzip.parquet"

    # Data structure for output
    output_data = {
        "mTrackType": [],
        "mPosX": [],
        "mPosY": [],
        "mPosZ": [],
        "mPt": [],
        "mEta": [],
        "mPhi": [],
        "mSign": [],
        "mMch2MCHMID": [],
        "mMch2MFTMCH": [],
        "mFwdDcaX": [],
        "mFwdDcaY": [],
        "mClsizeLayer": [],
        "mMeanClsizePerTrack": [],
    }

    # Open input file and iterate over trees
    print(f"[info] loading file form {input_path}")
    with uproot.open(input_path) as f_in:

        for key in f_in.keys():
            if "O2rtfwdpidall" in key:
                print(f"[warning] Skipping {key}")
                continue
            tree_in = f_in[f"{key}/O2rtfwdpidall/"]

            # Read branches as arrays
            print("[info] reading branches")
            fPt = tree_in["fPt"].array(library="np")
            fEta = tree_in["fEta"].array(library="np")
            fPhi = tree_in["fPhi"].array(library="np")
            fSign = tree_in["fSign"].array(library="np")
            fFwdDcaX = tree_in["fFwdDcaX"].array(library="np")
            fFwdDcaY = tree_in["fFwdDcaY"].array(library="np")
            fChi2MatchMCHMID = tree_in["fChi2MatchMCHMID"].array(library="np")
            fChi2MatchMCHMFT = tree_in["fChi2MatchMCHMFT"].array(library="np")
            fMftClusterSizesAndTrackFlags = tree_in[
                "fMftClusterSizesAndTrackFlags"
            ].array(library="np")
            fPosX = tree_in["fPosX"].array(library="np")
            fPosY = tree_in["fPosY"].array(library="np")
            fPosZ = tree_in["fPosZ"].array(library="np")
            fTrackType = tree_in["fTrackType"].array(library="np")

            # Process each entry with progress tracking using alive_bar
            with alive_bar(len(fPt), title=f"[info] processing {key}") as bar:
                for i in range(len(fPt)):
                    mean_clsize = 0
                    n_clusters = 0
                    clsize_layers = np.full(10, -999, dtype=np.int32)

                    for j in range(10):
                        size = (int(fMftClusterSizesAndTrackFlags[i]) >> (j * 6)) & 0x3F
                        clsize_layers[j] = size if size <= 62 else -999
                        if 0 < size <= 62:
                            mean_clsize += size
                            n_clusters += 1

                    # Calculate mean cluster size
                    mean_clsize_per_track = (
                        mean_clsize / n_clusters if n_clusters > 0 else -999
                    )

                    # Append data
                    output_data["mTrackType"].append(fTrackType[i])
                    output_data["mPosX"].append(fPosX[i])
                    output_data["mPosY"].append(fPosY[i])
                    output_data["mPosZ"].append(fPosZ[i])
                    output_data["mPt"].append(fPt[i])
                    output_data["mEta"].append(fEta[i])
                    output_data["mPhi"].append(fPhi[i])
                    output_data["mSign"].append(fSign[i])
                    output_data["mMch2MCHMID"].append(fChi2MatchMCHMID[i])
                    output_data["mMch2MFTMCH"].append(fChi2MatchMCHMFT[i])
                    output_data["mFwdDcaX"].append(fFwdDcaX[i])
                    output_data["mFwdDcaY"].append(fFwdDcaY[i])
                    output_data["mClsizeLayer"].append(clsize_layers)
                    output_data["mMeanClsizePerTrack"].append(mean_clsize_per_track)
                    bar()

    # Convert to DataFrame for parquet writing
    df_output = pd.DataFrame(
        {
            "mTrackType": output_data["mTrackType"],
            "mPosX": output_data["mPosX"],
            "mPosY": output_data["mPosY"],
            "mPosZ": output_data["mPosZ"],
            "mPt": output_data["mPt"],
            "mEta": output_data["mEta"],
            "mPhi": output_data["mPhi"],
            "mSign": output_data["mSign"],
            "mMch2MCHMID": output_data["mMch2MCHMID"],
            "mMch2MFTMCH": output_data["mMch2MFTMCH"],
            "mFwdDcaX": output_data["mFwdDcaX"],
            "mFwdDcaY": output_data["mFwdDcaY"],
            "mClsizeLayer": output_data["mClsizeLayer"],
            "mMeanClsizePerTrack": output_data["mMeanClsizePerTrack"],
        }
    )

    # Save as .gzip.parquet
    df_output.to_parquet(parquet_output_file, compression="gzip")
    print(f"[info] output saved to {parquet_output_file}")

    # Save as .root file
    with uproot.recreate(root_output_file) as f_out:
        f_out["treeTrackClusterSizeData"] = df_output
    print(f"[info] output saved to {root_output_file}")
    print()
    print("[info] finished")


if __name__ == "__main__":
    # Argument parser for command-line options
    parser = argparse.ArgumentParser(
        description="Process ROOT files to extract track cluster size data."
    )
    parser.add_argument("input_path", type=str, help="Path to the input ROOT file")
    parser.add_argument("output_path", type=str, help="Path to the output directory")
    parser.add_argument(
        "--suffix", type=str, default="", help="Optional suffix for the output file"
    )

    args = parser.parse_args()
    read_tree_fast(args.input_path, args.output_path, args.suffix)
