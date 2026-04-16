#!/usr/bin/env python3
from corsikaio import CorsikaFile
from corsikaio import CorsikaParticleFile
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--particle", type=str, help="Primary particle type")
parser.add_argument("-e", "--energy", type=str, help="Energy of primary particle ex. lgE_16.0")
parser.add_argument("-z", "--zenith", type=str, help="sin2_zenith ex. 0.0")
args = parser.parse_args()

primary = args.particle
energy = args.energy
sin2theta = args.zenith

base_dir = "/data/user/wkammeem/CORSIKA"
relative_path = f'GroundTotalParticles/{primary}_{energy}_{sin2theta}.npz'
output_path = os.path.join(base_dir, relative_path)

all_data = {
    "run": [],
    "zenith": [],
    "azimuth": [],
    "nMu": [],
    "nEP": [],
    "nElectron": []
}
print('run zenith azimuth TotalMuGround TotalEPGround TotalElectronGround')
for i in range(0,200):
    runnums = i
    try:
        file_path = (f'/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energy}/sin2_{sin2theta}/{runnums:06d}/DAT{runnums:06d}')
        # Open the CORSIKA DAT file
        with CorsikaParticleFile(file_path, thinning=True) as f:

            run_number = f.run_header[1]

            for event in f: #loop for 1 shower

                nMupm = 0
                nEpm = 0
                nElectron = 0

                header = event.header
                zenith = header[10]   # Index 11 in CORSIKA manual unit in radian
                azimuth = header[11]  # Index 12 in OSRSIKA manual unit in radian
                
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
                
                nMupm = weight[idx_mu].sum() # number of total +-muon at ground
                nEpm = weight[idx_epm].sum() # number of total +-e at ground
                nElectron = weight[idx_electron].sum() # number of total electron at ground
            print(runnums, zenith, azimuth, nMupm, nEpm, nElectron)

                
    except FileNotFoundError:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue

    except OSError:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue
    
    
    all_data["run"].append(runnums)
    all_data["zenith"].append(zenith)
    all_data["azimuth"].append(azimuth)
    all_data["nMu"].append(nMupm)
    all_data["nEP"]. append(nEpm)
    all_data["nElectron"].append(nElectron)


output_dir = os.path.dirname(output_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# save everything in a compressed npz file
np.savez_compressed(output_path, **all_data)
print(f"Saved all data to {output_path}")