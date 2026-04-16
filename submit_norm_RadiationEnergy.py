#!/usr/bin/env python3

import numpy as np
import sys
import argparse
import os
#sys.path.append('../')
from functions import *

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--particle", type=str, help="Primary particle type")
parser.add_argument("-e", "--energy", type=str, help="Energy of primary particle ex. lgE_16.0")
parser.add_argument("-z", "--zenith", type=str, help="sin2_zenith ex. 0.0")
args = parser.parse_args()

p = args.particle
e = args.energy
sin2theta = args.zenith

rawdata = "/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern" 
fp_RadE_norm = lambda primary, energy, sin2theta: f'radEnergy/norm_sintheta2/{primary}_{energy}_{sin2theta}.npz'
fp_RadE = lambda primary, energy, sin2theta: f'radEnergy/raw/{primary}_{energy}_{sin2theta}.npz'
fp_inp = lambda primary, energy, sin2theta, runnum: f'{rawdata}/{primary}/{energy}/sin2_{sin2theta}/{runnum:06d}/SIM{runnum:06d}.inp'

base_dir = "/data/user/wkammeem/CORSIKA"
subout_path = fp_RadE_norm(p, e, sin2theta)
output_path = os.path.join(base_dir, subout_path)

subin_path1 = fp_RadE(p, e, sin2theta)
input_path1 = os.path.join(base_dir, subin_path1)

all_data = {
    "run": [],
    "alpha": [],
    "radE(eV)": [],
    "radE_filtered(eV)": [],
    "radE_vxB(eV)": [],
    "radE_filtered_vxB(eV)": [],
    
}

fileRadE = np.load(input_path1, allow_pickle= True)

RadE = fileRadE['radE(eV)']
RadE_f = fileRadE['radE_filtered(eV)']

RadE_vxB = fileRadE['radE_vxB(eV)']
RadE_f_vxB = fileRadE['radE_filtered_vxB(eV)']

# Normalization with sin(alpha)
for i in range(200):
    r = i
    subin_path2 = fp_inp(p, e, sin2theta,r)
    input_path2 = os.path.join(base_dir, subin_path2)
    finp = input_path2
    
    try:
        with open(finp) as f:
            for line in f:
                parts = line.split()
                if parts[0] == "THETAP":
                    thetap = float(parts[1])
                if parts[0] == "PHIP":
                    phip = float(parts[1])
                elif parts[0] == "MAGNET":
                    Bx, Bz = float(parts[1]), -float(parts[2]) # input Bz possitive downwards
                    By = 0

    except FileNotFoundError:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue
                
    theta = np.deg2rad(thetap)
    phi   = np.deg2rad(phip)
    
    rad  = 1 # basis v
    vx, vy, vz =  rad* np.sin(theta)*np.cos(phi), rad*np.sin(theta)*np.sin(phi), -rad*np.cos(theta)
    v = [vx, vy, vz]
    B = [Bx, By, Bz]
    
    magV, magB = np.linalg.norm(v), np.linalg.norm(B)
    
    alpha = np.arccos(np.dot(B,v)/(magV*magB))

    RadEnew = RadE[i]/(np.sin(alpha)**2)
    RadEnew_f = RadE_f[i]/(np.sin(alpha)**2)

    RadEnew_vxB = RadE_vxB[i]/(np.sin(alpha)**2)
    RadEnew_f_vxB = RadE_f_vxB[i]/(np.sin(alpha)**2)

    # store data in dictionary
    all_data["run"].append(r)
    all_data["radE(eV)"].append(RadEnew)
    all_data["radE_filtered(eV)"].append(RadEnew_f)
    all_data["alpha"].append(alpha)
    all_data["radE_vxB(eV)"].append(RadEnew_vxB)
    all_data["radE_filtered_vxB(eV)"].append(RadEnew_f_vxB)

output_dir = os.path.dirname(output_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# save everything in a compressed npz file
np.savez_compressed(output_path, **all_data)
print(f"Saved all data to {output_path}")
print(f'{p} {e} {sin2theta} normalization done.')

            