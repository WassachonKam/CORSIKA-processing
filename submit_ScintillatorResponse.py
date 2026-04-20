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

# constants
mumass = 0.105658 #muon mass in GeV/c^2
emass = 5.11e-4 #electron mass in GeV/c^2

############################# DATA PREPARATION #############################


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

    except (OSError, IOError) as err:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue
    except StopIteration:
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

    # Kinetic energy in GeV
    Ek_mu = Ekin(px[idx_mu], py[idx_mu], pz[idx_mu], mumass)
    Ek_epm = Ekin(px[idx_epm], py[idx_epm], pz[idx_epm], emass) 
    Ek_electron = Ekin(px[idx_electron], py[idx_electron], pz[idx_electron], emass) 

    # zenith angle for normalization
    theta = np.deg2rad(thetap)
    
    ############################# SCINTILLATOR RESPONSE #############################   

    total_Edep_all = 0 

    # Define array of particle types
    Ek_array = [Ek_electron, Ek_mu]
    weight_array = [weight[idx_electron], weight[idx_mu]]

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
    
output_dir = "/data/user/wkammeem/CORSIKA/TotalEdepScint"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_path = f"{output_dir}/{primary}_{energy}_{sin2theta}.npz"
np.savez_compressed(output_path, **all_data)
print(f"Saved files to {output_path}")

        
