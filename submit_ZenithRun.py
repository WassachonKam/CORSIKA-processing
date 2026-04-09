#!/usr/bin/env python3
from corsikaio import CorsikaFile
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
relative_path = f'ZenithAngle/{primary}_{energy}_{sin2theta}.npz'
output_path = os.path.join(base_dir, relative_path)

all_data = {
    "run": [],
    "zenith": [],
    "azimuth": []
}

for i in range(0,200):
    runnums = i

    try:
        file_path = (f'/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energy}/sin2_{sin2theta}/{runnums:06d}/DAT{runnums:06d}')
        # Open the CORSIKA DAT file
        with CorsikaFile(file_path, thinning=True) as f:

            run_number = f.run_header[1]

            for event in f:

                header = event.header
                
                zenith = header[10]   # Index 11 in CORSIKA manual unit in radian
                azimuth = header[11]  # Index 12 in OSRSIKA manual unit in radian
                
                print(runnums, zenith, azimuth)
    except FileNotFoundError:
        continue
    
    all_data["run"].append(runnums)
    all_data["zenith"].append(zenith)
    all_data["azimuth"].append(azimuth)


output_dir = os.path.dirname(output_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# save everything in a compressed npz file
np.savez_compressed(output_path, **all_data)
print(f"Saved all data to {output_path}")