# import functions
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import joblib


threshR = 250 # minimum radius in meters
Xmaxfiltering = True # apply Xmax underground removal
NmuNorm = False # apply Nmu normalization


def load_data_as_df (p, e, sin2theta, threshR, Xmaxfiltering):

    # file path
    fMuon = f'../GroundTotalParticles/{p}_{e}_{sin2theta}.npz'
    fXmax = f'../Xmax/{p}_{e}_sin2_{sin2theta}.dat'
    fRadE = f'../radEnergy/norm_sintheta2/{p}_{e}_{sin2theta}.npz'
    fEdep = f'../TotalEdepScint/{threshR}m/{p}_{e}_{sin2theta}.npz'
    
    # load data
    fileRadE = np.load(fRadE, allow_pickle=True)
    fileMuon = np.load(fMuon, allow_pickle= True)
    fileEdep = np.load(fEdep, allow_pickle=True)
    fileXmax = np.loadtxt(fXmax)

    # determine parameters
    energy = fileEdep["primaryE"] #  primary particle energy in GeV
    Xmax  = fileXmax[:,6] # Xmax 
    Nep_Xmax  = fileXmax[:,5] # number of total +-e at Xmax
    Ne_ground = fileMuon['nEP'] # number of total +-mu at ground
    Nmu_ground = fileMuon['nMu'] # number of total +-mu at ground
    zenith = fileMuon['zenith'] # zenith angle in rad
    Edep = fileEdep["Edep_tot"] # Nscint number of total scintillation photons
    RadE = fileRadE['radE_filtered(eV)']/1e9 #Erad in GeV
    alpha = fileRadE['alpha'] # angle between shower and magntic field

    sinalpha = np.sin(alpha)
    costheta = np.cos(zenith)

    X0 = 697.6
    Xv = X0/ costheta

    mask = Xmax >= 0 # keep only positive Xmax
    if Xmaxfiltering == True:
        mask = (Xmax >=0) & (Xmax <= Xv)

    # store all parameters in data frame
    df = pd.DataFrame({
        'particle': p,                          # primary particle
        'energy': energy[mask],                 # primary energy
        'costheta': costheta[mask],             # cos(zenith)
        'Xmax': Xmax[mask],                     # Xmax
        'sinalpha': sinalpha[mask],             # sin(alpha), alpha = angle between shower and magntic field
        'Edep': Edep[mask],                     # number of total scintillation photons
        'Erad': RadE[mask],                     # radiation energy in eV
        'Nmu': Nmu_ground[mask],                # number of muon at ground
        'Ne' : Nep_Xmax[mask],                  # number of electron and positron at Xmax
        'Ne_ground': Ne_ground[mask]            # number of electron and positron at ground    
    })
    
    return df

particles = ['proton', 'iron', 'helium', 'oxygen']
energies =  [f"lgE_{x/10:.1f}" for x in range(160,181)]
sin2thetas = [x/10 for x in range(0, 10)]

all_df = []
print('starting to download files')
for p in particles:
    for e in energies:
        for z in sin2thetas:
            try:
                df_bin = load_data_as_df (p, e, z, threshR, Xmaxfiltering)
                all_df.append(df_bin)
            except FileNotFoundError:
                continue 

# gather all bins into one big data frame
final_df = pd.concat(all_df, ignore_index=True)
final_df.replace([np.inf, -np.inf], np.nan, inplace=True) #replace inf values with nan
final_df.dropna(inplace=True) # remove all nans
# final_df.to_parquet(f'data_for_regression_{threshR}_Xmax_{Xmaxfiltering}_NmuNorm_{NmuNorm}.parquet', index=False)
final_df.to_parquet(f'DataForRegression/data_for_regression_{threshR}_Xmax_{Xmaxfiltering}.parquet', index=False)
print(final_df.head())
print(f"Total rows (events): {len(final_df)}/168,000")
print('file is saved')
