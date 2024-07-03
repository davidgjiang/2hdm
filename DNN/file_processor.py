# Import necessary libraries
import pandas as pd
import numpy as np
import uproot as up
import os
from tqdm import tqdm

# Functions
def get_root_files(directory):
    """Get all .root file paths from a directory."""
    return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".root")]

def load_signal_data(paths, description,load_branches):
    """Load signal data from a list of paths into a DataFrame. Adds a column named "Class" and sets all values to 1."""
    signal_data = []
    for path in tqdm(paths, desc=f"Loading {description}"):
        reco_tree = up.open(path)["reco"]
        reco_df = reco_tree.arrays(expressions=load_branches,library="pd")
        temp_df = reco_df[(reco_df["nBjets_NOSYS"] > 1) & (reco_df["nBjets_NOSYS"] + reco_df["nLjets_NOSYS"] > 3)]  # preselection cuts
        signal_data.append(temp_df)
    df = pd.concat(signal_data).reset_index(drop=True)
    df["Class"] = 1 
    return df
    
def load_ttbar_data(path,load_branches):
    """Load ttbar data from a list of paths into a DataFrame. Adds a column named "Class" and sets all values to 0."""
    reco_tree = up.open(path)["reco"]
    total_entries = reco_tree.num_entries
    chunk_size = 10000  
    ttbar_data = []

    with tqdm(total=total_entries, desc=f"Loading ttbar file", unit="entries") as pbar:
        for start in range(0, total_entries, chunk_size):
            end = min(start + chunk_size, total_entries)
            reco_df = reco_tree.arrays(entry_start=start, entry_stop=end, expressions=load_branches,library="pd")
            temp_df = reco_df[(reco_df["nBjets_NOSYS"] > 1) & (reco_df["nBjets_NOSYS"] + reco_df["nLjets_NOSYS"] > 3)]  # preselection cuts
            ttbar_data.append(temp_df)
            pbar.update(chunk_size)

    # Combine all chunks into a single DataFrame
    df = pd.concat(ttbar_data).reset_index(drop=True)
    df["Class"] = 0
    return df
    
def deltaR(eta1, phi1, eta2, phi2):
    """Returns the deltaR value."""
    d_eta = eta1 - eta2
    d_phi = (phi1 - phi2 + np.pi) % (2 * np.pi) - np.pi
    return np.sqrt(d_eta**2 + d_phi**2)

def calculate_deltaR(df, object1_prefix, object2_prefix):
    """Calculates deltaR values between object1 and object2 and adds it as a column to the DataFrame."""
    col_name = f"deltaR_{object1_prefix}_{object2_prefix}"
    eta1 = df.get(f"{object1_prefix}_eta_NOSYS", df.get(f"{object1_prefix}_eta_fitted_NOSYS", np.nan))
    phi1 = df.get(f"{object1_prefix}_phi_NOSYS", df.get(f"{object1_prefix}_phi_fitted_NOSYS", np.nan))
    eta2 = df.get(f"{object2_prefix}_eta_NOSYS", df.get(f"{object2_prefix}_eta_fitted_NOSYS", np.nan))
    phi2 = df.get(f"{object2_prefix}_phi_NOSYS", df.get(f"{object2_prefix}_phi_fitted_NOSYS", np.nan))
    mask = (eta1 == -999) | (eta2 == -999) | (phi1 == -999) | (phi2 == -999)
    delta_r = deltaR(eta1, phi1, eta2, phi2)
    delta_r[mask] = np.nan
    df[col_name] = delta_r

# File Locations
high_mass = "/data/dajiang/2HDM/ML/samples/signals/high_mass/"
medium_mass = "/data/dajiang/2HDM/ML/samples/signals/medium_mass/"
low_mass = "/data/dajiang/2HDM/ML/samples/signals/low_mass/"

# Save only the useful branches (4-vectors). Also include weights for later use in TRExFitter. Need to save jet number branches for selection cuts too.
branches = [
    # preselection cut branches
    "nBjets_NOSYS", "nLjets_NOSYS",
    # final state particle 4-vectors
    "bjet1_pt_NOSYS","bjet1_phi_NOSYS","bjet1_eta_NOSYS","bjet1_mass_NOSYS",
    "bjet2_pt_NOSYS","bjet2_phi_NOSYS","bjet2_eta_NOSYS","bjet2_mass_NOSYS",
    "ljet1_pt_NOSYS","ljet1_phi_NOSYS","ljet1_eta_NOSYS","ljet1_mass_NOSYS",
    "ljet2_pt_NOSYS","ljet2_phi_NOSYS","ljet2_eta_NOSYS","ljet2_mass_NOSYS",
    "ljet3_pt_NOSYS","ljet3_phi_NOSYS","ljet3_eta_NOSYS","ljet3_mass_NOSYS",
    "ljet4_pt_NOSYS","ljet4_phi_NOSYS","ljet4_eta_NOSYS","ljet4_mass_NOSYS",
    "lepton_pt_NOSYS","lepton_eta_NOSYS","lepton_phi_NOSYS","lepton_mass_NOSYS","lepton_charge_NOSYS","lepton_pdgid_NOSYS",
    "MET_NOSYS","MET_phi_NOSYS",
    # intermediate particle 4-vectors
    "A_pt_fitted_NOSYS","A_phi_fitted_NOSYS","A_eta_fitted_NOSYS","A_mass_fitted_NOSYS",
    "Hp_pt_fitted_NOSYS","Hp_phi_fitted_NOSYS","Hp_eta_fitted_NOSYS","Hp_mass_fitted_NOSYS",
    "top_pt_fitted_NOSYS","top_phi_fitted_NOSYS","top_eta_fitted_NOSYS","top_mass_fitted_NOSYS",
    "WfromH_pt_fitted_NOSYS","WfromH_phi_fitted_NOSYS","WfromH_eta_fitted_NOSYS","WfromH_mass_fitted_NOSYS",
    "WfromTop_pt_fitted_NOSYS","WfromTop_phi_fitted_NOSYS","WfromTop_eta_fitted_NOSYS","WfromTop_mass_fitted_NOSYS",
    "Wb_nonTop_pt_fitted_NOSYS","Wb_nonTop_phi_fitted_NOSYS","Wb_nonTop_eta_fitted_NOSYS","Wb_nonTop_mass_fitted_NOSYS",
    # weights
    "weight_btagSF_DL1dv01_FixedCutBEff_85_NOSYS","weight_beamspot","weight_mc_NOSYS", "weight_btagSF_DL1dv01_Continuous_NOSYS",
    "weight_pileup_NOSYS","weight_jvt_effSF_NOSYS"
    ]

# Compute all deltaR combinations: 78 total
combinations = [
    ("A", "Hp"), ("A", "top"), ("A", "WfromTop"), ("A", "WfromH"), ("A", "Wb_nonTop"),
    ("Hp", "top"), ("Hp", "WfromTop"), ("Hp", "WfromH"), ("Hp", "Wb_nonTop"),
    ("top", "WfromTop"), ("top", "WfromH"), ("top", "Wb_nonTop"),
    ("WfromTop", "WfromH"), ("WfromTop", "Wb_nonTop"),
    ("WfromH", "Wb_nonTop"),
    ("bjet1", "bjet2"), ("ljet1", "ljet2"), ("bjet1", "ljet1"), ("bjet1", "ljet2"),
    ("bjet2", "ljet1"), ("bjet2", "ljet2"), ("bjet1", "ljet3"), ("bjet1", "ljet4"),
    ("bjet2", "ljet3"), ("bjet2", "ljet4"), ("ljet1", "ljet3"), ("ljet1", "ljet4"),
    ("ljet2", "ljet3"), ("ljet2", "ljet4"), ("ljet3", "ljet4"),
    ("bjet1", "lepton"), ("bjet2", "lepton"), ("ljet1", "lepton"),
    ("ljet2", "lepton"), ("ljet3", "lepton"), ("ljet4", "lepton"),
    ("A", "bjet1"), ("A", "bjet2"), ("A", "ljet1"), ("A", "ljet2"),
    ("A", "ljet3"), ("A", "ljet4"), ("A", "lepton"),
    ("top", "bjet1"), ("top", "bjet2"), ("top", "ljet1"), ("top", "ljet2"),
    ("top", "ljet3"), ("top", "ljet4"), ("top", "lepton"),
    ("WfromTop", "bjet1"), ("WfromTop", "bjet2"), ("WfromTop", "ljet1"), ("WfromTop", "ljet2"),
    ("WfromTop", "ljet3"), ("WfromTop", "ljet4"),("WfromTop", "lepton"), 
    ("Wb_nonTop", "bjet1"), ("Wb_nonTop", "bjet2"), ("Wb_nonTop", "ljet1"), ("Wb_nonTop", "ljet2"),
    ("Wb_nonTop", "ljet3"), ("Wb_nonTop", "ljet4"), ("Wb_nonTop", "lepton"),
    ("Hp", "bjet1"), ("Hp", "bjet2"), ("Hp", "ljet1"), ("Hp", "ljet2"),
    ("Hp", "ljet3"), ("Hp", "ljet4"), ("Hp", "lepton"),
    ("WfromH", "bjet1"), ("WfromH", "bjet2"), ("WfromH", "ljet1"), ("WfromH", "ljet2"),
    ("WfromH", "ljet3"), ("WfromH", "ljet4"), ("WfromH", "lepton"),
]

if __name__ == "__main__":
    print("*** File Processor ***")
    print("Converts signal and ttbar NTuples to parquet files and saves relevant branches to use in Machine Learning algorithms.")

    # Get the file paths for the samples
    high_mass_paths = get_root_files(high_mass)
    medium_mass_paths = get_root_files(medium_mass)
    low_mass_paths = get_root_files(low_mass)
    ttbar_path = "/data/dajiang/2HDM/ML/samples/ttbar_nom/mc20a/user.rashbypi.410470.PhPy8EG.DAOD_PHYSLITE.e6337_s3681_r13144_p5855.2024-06-03-v0.2_output/user.rashbypi.39661650._000001.output.root"


    # Process data and store into DataFrames
    high_mass_df = load_signal_data(high_mass_paths,"high mass files",branches)
    medium_mass_df = load_signal_data(medium_mass_paths,"medium mass files",branches)
    low_mass_df = load_signal_data(low_mass_paths,"low mass files",branches)
    ttbar_df = load_ttbar_data(ttbar_path,branches)

    # Calculate the deltaR pairs and store them as extra columns in the DataFrame
    for comb in tqdm(combinations,desc="Calculating deltaR values"):
        calculate_deltaR(high_mass_df, comb[0], comb[1])
        calculate_deltaR(medium_mass_df, comb[0], comb[1])
        calculate_deltaR(low_mass_df, comb[0], comb[1])
        calculate_deltaR(ttbar_df, comb[0], comb[1])

    # Save the DataFrames to parquet files
    high_mass_df.to_parquet("/home/dajiang/2hdm/DNN/processed_files/high_mass/high_mass_sig_df.parquet", engine="fastparquet")
    medium_mass_df.to_parquet("/home/dajiang/2hdm/DNN/processed_files/medium_mass/medium_mass_sig_df.parquet", engine="fastparquet")
    low_mass_df.to_parquet("/home/dajiang/2hdm/DNN/processed_files/low_mass/low_mass_sig_df.parquet", engine="fastparquet")
    ttbar_df.to_parquet("/home/dajiang/2hdm/DNN/processed_files/ttbar/ttbar_df.parquet", engine="fastparquet")

    print("Processing done.")
