#!/usr/bin/env python3
# Usage: python3 submit_ScintillatorResponse.py -p proton -e lgE_16.2 -z 0.0 
import argparse
import numpy as np
import matplotlib.pyplot as plt
from corsikaio import CorsikaParticleFile
import pandas as pd
import random
import os

############################# FUNCTIONS AND CONSTANTS #############################

# kinetic energy 
def Ekin(px,py,pz,m):
    return np.sqrt(px**2 + py**2 + pz**2 + m**2) -m 

def R(x,y):
    return np.sqrt(x**2 + y**2)
# constants
mumass = 0.105658 #muon mass in GeV/c^2
emass = 5.11e-4 #electron mass in GeV/c^2

############################# DATA PREPARATION #############################

maskR  = False
minR   = 0
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--particle", type=str, help="Primary particle type")
parser.add_argument("-e", "--energy", type=str, help="Energy of primary particle ex. lgE_16.0")
parser.add_argument("-z", "--zenith", type=str, help="sin2_zenith ex. 0.0")
args = parser.parse_args() 

primary = args.particle
energy = args.energy
sin2theta = args.zenith

print(f"processing deposited energy {primary} {energy} {sin2theta}")

# dictionary for storing data
all_data = {
            "primaryE": [],
            "Edep_e": [],
            "Edep_mu": [],
            "Edep_epm": [],
            "Edep_tot": []

        }

# download digitized data from Agnieszka's thesis
dfe = pd.read_csv('ScintillatorResponse/electron_0.0.csv')
dmu = pd.read_csv('ScintillatorResponse/muon_0.0.csv')
digitized_files = [dfe, dmu]
plt_title = [rf'$e^{{-}}$', rf'$\mu^{{\pm}}$']
savefile = ['electron', 'muon']


runmin, runmax = 0, 199

for run in range(runmin, runmax+1):

    file_path = (f'/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energy}/sin2_{sin2theta}/{run:06d}/DAT{run:06d}')
    file_input = (f"/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energy}/sin2_{sin2theta}/{run:06d}/SIM{run:06d}.inp")

    # print(run)

    try:
        with open(file_input) as f:
            for line in f:
                parts = line.split()
                if parts[0] == "THETAP":
                    thetap = float(parts[1]) # zenith angle (deg)
                if parts[0] == "ERANGE": 
                    primE = float(parts[1]) # energy of primary particle (GeV)
    except FileNotFoundError:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue

    try:
        with CorsikaParticleFile(file_path, thinning= True) as file:
            # we only have one event per file, we can grab it like this
            event = next(file)

    except (OSError, IOError, StopIteration, IndexError, ValueError) as err:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue

    
    all_data["primaryE"].append(primE)

    # get the particle info
    #print(event.particles.dtype.names)
    particle_id = event.particles['particle_description'] // 1000 # corsika particle ID
    x = event.particles['x'] # x coordinate
    y = event.particles['y'] # y coordinate
    px = event.particles['px'] # momentum component in x direction in GeV
    py = event.particles['py'] # momentum component in y direction in GeV
    pz = event.particles['pz'] # momentum component in z direction in GeV
    weight = event.particles['thinning_weight'] # particle weight 


    # get indices of muons
    idx_mu = np.where((particle_id == 5) | (particle_id == 6))
    idx_epm = np.where((particle_id == 2) | (particle_id == 3))
    idx_electron = np.where((particle_id == 3))

    # particle position
    x_mu, y_mu = x[idx_mu], y[idx_mu]
    x_epm, y_epm = x[idx_epm], y[idx_epm]
    x_e, y_e = x[idx_electron], y[idx_electron]

    R_mu = R(x_mu, y_mu) 
    R_epm = R(x_epm, y_epm)
    R_e = R(x_e, y_e)

    # Kinetic energy in GeV
    Ek_mu = Ekin(px[idx_mu], py[idx_mu], pz[idx_mu], mumass)
    Ek_epm = Ekin(px[idx_epm], py[idx_epm], pz[idx_epm], emass) 
    Ek_e = Ekin(px[idx_electron], py[idx_electron], pz[idx_electron], emass) 

    weight_mu = weight[idx_mu]
    weight_epm = weight[idx_epm]
    weight_e = weight[idx_electron]

    # masking Ek and weight R > 200 m 
    if maskR == True:
        mask_mu = (R_mu > minR*10**3) # unit in cm
        mask_epm = (R_epm > minR*10**3) # unit in cm 
        mask_e = (R_e > minR*10**3) # unit in cm 

        Ek_mu = Ek_mu[mask_mu]
        Ek_epm = Ek_epm[mask_epm]
        Ek_e = Ek_e[mask_e]

        weight_mu = weight_mu[mask_mu]
        weight_epm = weight_epm[mask_epm]
        weight_e = weight_e[mask_e]

    # zenith angle for normalization
    theta = np.deg2rad(thetap)
    
    ############################# SCINTILLATOR RESPONSE #############################   

    total_Edep_all = 0 

    # Define array of particle types
    Ek_array = [Ek_e, Ek_mu]
    weight_array = [weight_e, weight_mu]

    for i in range(len(digitized_files)):

        df = digitized_files[i]
        df_x = df['x']
        df_y = df[' y'] / np.cos(theta) # normalized by cos(zenith)

        # interpolate digitized data
        interpolate_x = np.linspace(min(df_x), max(df_x), num = 400)
        interpolate_y = np.interp(interpolate_x, df_x, df_y)

        # Deposited energy as a function of Ek from CORSIKA file
        logEk = np.log10(Ek_array[i])
        corsika_x = logEk
        corsika_y = np.interp(corsika_x, df_x, df_y)

        # Weight deposited energy and calculate the total amount
        weightedEdep = corsika_y * weight_array[i] # apply weight factor 
        total_Edep = sum(weightedEdep) # total deposited energy of every particles shared the same type
        total_Edep_all += total_Edep # total deposited energy of every particles types

        if i == 0: all_data["Edep_e"].append(total_Edep) 
        elif i == 1: all_data["Edep_mu"].append(total_Edep)
    
    all_data["Edep_tot"].append(total_Edep_all)
    
output_dir = f"/data/user/wkammeem/CORSIKA/TotalEdepScint/{minR}m"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_path = f"{output_dir}/{primary}_{energy}_{sin2theta}.npz"
np.savez_compressed(output_path, **all_data)
print(f"Saved files to {output_path}")

        
