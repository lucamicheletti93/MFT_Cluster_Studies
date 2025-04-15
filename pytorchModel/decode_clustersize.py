import uproot

# Apri il file ROOT
file = uproot.open("/data/shared/mft_pid/data/AO2D_matched_tracks.root")

# Mostra tutti gli oggetti presenti nel file (alberi, istogrammi, ecc.)
file.keys()

# Apri il tree corretto (ignora il ;1 finale)
tree = file["DF_2388938487428736/O2rtfwdpidall"]

# Mostra i nomi delle colonne (branches)
tree.keys()


import numpy as np
import pandas as pd
import uproot
import os
from alive_progress import alive_bar

# Definizione dei percorsi
input_path = "/data/shared/mft_pid/data/AO2D_matched_tracks.root"
output_path = "/home/rmbelcas/converted"
suffix = "DF_2388938487428736"
tree_name = "DF_2388938487428736/O2rtfwdpidall" 

# Crea la directory di output se non esiste
os.makedirs(output_path, exist_ok=True)

def read_tree_fast(input_path, output_path, tree_name, suffix=""):
    root_output_file = f"{output_path}/treeMftPidTraining{suffix}.root"
    parquet_output_file = f"{output_path}/treeMftPidTraining{suffix}.gzip.parquet"

    print(f"[INFO] Loading file from {input_path}")
    with uproot.open(input_path) as f_in:
        # Leggi direttamente il tree specificato
        tree_in = f_in[tree_name]

        # Leggi i branch come array
        print("[INFO] Reading branches...")
        fPt = tree_in["fPt"].array(library="np")
        fEta = tree_in["fEta"].array(library="np")
        fPhi = tree_in["fPhi"].array(library="np")
        fSign = tree_in["fSign"].array(library="np")
        fFwdDcaX = tree_in["fFwdDcaX"].array(library="np")
        fFwdDcaY = tree_in["fFwdDcaY"].array(library="np")
        fChi2MatchMCHMID = tree_in["fChi2MatchMCHMID"].array(library="np")
        fChi2MatchMCHMFT = tree_in["fChi2MatchMCHMFT"].array(library="np")
        fMftClusterSizesAndTrackFlags = tree_in["fMftClusterSizesAndTrackFlags"].array(library="np")
        fPosX = tree_in["fPosX"].array(library="np")
        fPosY = tree_in["fPosY"].array(library="np")
        fPosZ = tree_in["fPosZ"].array(library="np")
        fTrackType = tree_in["fTrackType"].array(library="np")

        # La colonna supervisionata (fMcDecision) potrebbe non esistere in tutti i dataset
        try:
            fMcDecision = tree_in["fMcDecision"].array(library="np")
            has_fMcDecision = True
        except KeyError:
            print("[WARNING] Colonna 'fMcDecision' non trovata.")
            has_fMcDecision = False

        # Dizionario per raccogliere i dati di output
        output_data = {
            "mTrackType": fTrackType,
            "mPosX": fPosX,
            "mPosY": fPosY,
            "mPosZ": fPosZ,
            "mPt": fPt,
            "mEta": fEta,
            "mPhi": fPhi,
            "mSign": fSign,
            "mMch2MCHMID": fChi2MatchMCHMID,
            "mMch2MFTMCH": fChi2MatchMCHMFT,
            "mFwdDcaX": fFwdDcaX,
            "mFwdDcaY": fFwdDcaY,
        }

        if has_fMcDecision:
            output_data["mMcDecision"] = fMcDecision

        # Aggiunta colonne per i cluster
        for i in range(10):
            output_data[f"mClsizeLayer{i}"] = []

        output_data["mMeanClsizePerTrack"] = []

        print("[INFO] Processing data...")

        with alive_bar(len(fPt), title="[INFO] Processing events") as bar:
            for i in range(len(fPt)):
                mean_clsize = 0
                n_clusters = 0
                clsize_layers = np.full(10, -999, dtype=np.int32)

                # Estrazione delle dimensioni dei cluster
                for j in range(10):
                    size = (int(fMftClusterSizesAndTrackFlags[i]) >> (j * 6)) & 0x3F
                    clsize_layers[j] = size if size <= 62 else -999
                    if 0 < size <= 62:
                        mean_clsize += size
                        n_clusters += 1

                # Calcolo della media della dimensione dei cluster
                mean_clsize_per_track = mean_clsize / n_clusters if n_clusters > 0 else -999

                # Scrittura nei dizionari
                for j in range(10):
                    output_data[f"mClsizeLayer{j}"].append(clsize_layers[j])

                output_data["mMeanClsizePerTrack"].append(mean_clsize_per_track)

                bar()

        # Conversione in DataFrame
        df_output = pd.DataFrame(output_data)

        # Salvataggio su parquet
        df_output.to_parquet(parquet_output_file, compression="gzip")
        print(f"[INFO] Output saved to {parquet_output_file}")

        # Salvataggio su ROOT (opzionale)
        with uproot.recreate(root_output_file) as f_out:
            f_out["treeTrackClusterSizeData"] = df_output

        print(f"[INFO] Output saved to {root_output_file}")
        print("[INFO] Finished.")

# ✅ Lancia la funzione
read_tree_fast(input_path, output_path, tree_name, suffix)
