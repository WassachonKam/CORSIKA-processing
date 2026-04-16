#!/usr/bin/env python3

from scipy import signal, fft, constants, optimize
import argparse
import os
from functions import *


# Radiation energy calculation in both total radiation energy and separate components of electric field
# Usage: python3 RadiationEnergy.py -p proton -e lgE_16.2 -z 0.0 

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--particle", type=str, help="Primary particle type")
parser.add_argument("-e", "--energy", type=str, help="Energy of primary particle ex. lgE_16.0")
parser.add_argument("-z", "--zenith", type=str, help="sin2_zenith ex. 0.0")
args = parser.parse_args()

p = args.particle
e = args.energy
sin2theta = args.zenith

# antenna positions
ant_x = lambda xpos: xpos.astype(float)/1e2 # x position in m
ant_y =  lambda ypos: ypos.astype(float)/1e2 # y position in m

b, a = signal.butter(5, [70e6, 350e6] , 'bp', fs=5e9)

base_dir = "/data/user/wkammeem/CORSIKA"
relative_path = fp_RadE(p, e, sin2theta)
output_path = os.path.join(base_dir, relative_path)

runmin, runmax = 0, 199

# create dictionaries to store all data
all_data = {
    "run": [],
    "radE(eV)": [],
    "radE_filtered(eV)": [],
    "antx": [],
    "anty": [],
    "ef": [],
    "eff": [],
    "radE_vxB(eV)": [],
    "radE_filtered_vxB(eV)": [],
    "radE_vxvxB(eV)": [],
    "radE_filtered_vxvxB(eV)": []
}

print(f"processing radiation energy {p} {e} {sin2theta} ")
for run in range(runmin, runmax+1):

    # x-y antenna position
    try:
        flist = fp_list(p, e, sin2theta, run)
        xypos = pd.read_csv(flist, sep=r"\s+", header=None)
    except FileNotFoundError:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue   # skip this run entirely
    
    antx = ant_x(xypos[2])
    anty = ant_y(xypos[3])

    try:
        finp = fp_inp(p, e, sin2theta, run)
    except FileNotFoundError:
        for key in all_data.keys():
            all_data[key].append(np.nan)
        continue
    
    with open(finp) as f:
        for line in f:
            parts = line.split()
            if parts[0] == "THETAP":
                thetap = float(parts[1])
            if parts[0] == "PHIP":
                phip = float(parts[1])
            elif parts[0] == "MAGNET":
                Bx, Bz = float(parts[1]), -float(parts[2]) # input Bz possitive downward
                By = 0
                
    # shower direction
    theta = np.deg2rad(thetap)
    phi   = np.deg2rad(phip)
    
    # magnetic field 
    B = np.array([Bx, By, Bz])

    Nant = 160 #number of total antennas
    # create empty arrays of energy fluence, filtered energy fluence, vxB, and vxvxB
    ef, eff, ef_vxB, eff_vxB, ef_vxvxB, eff_vxvxB,  = (np.zeros(Nant) for _ in range(6))

    # 160 antenna (Nant = 160) loop
    for i, ant_no in enumerate(range (1,Nant+1)):

        try: radio_dat = np.loadtxt(fp_radio(p, e, sin2theta, run, ant_no)) #load raw_ant_{ant}.dat
        except FileNotFoundError: continue
            
        ### Ex, Ey, Ez components
        # extract data
        time = radio_dat[:,0]
        E_SI = 29979 # electric field converter from cgs unit to SI unit 1 statV/cm≈29,979 V/m
        Ex, Ey, Ez = radio_dat[:,1]*E_SI, radio_dat[:,2]*E_SI, radio_dat[:,3]*E_SI # north, west, vertical electric field in V/m

        #filter electric field using signal.butter
        # b, a = signal.butter(5, [70e6, 350e6] , 'bp', fs=5e9)
        Ex_f, Ey_f, Ez_f = (signal.filtfilt(b, a, i) for i in (Ex, Ey, Ez))
        Emag2, Emag2_f = magE2(Ex, Ey, Ez), magE2(Ex_f, Ey_f, Ez_f)  # |E|^2 raw and filtered
        E_sum, E_sum_f = sum(Emag2), sum(Emag2_f) # sum(|E|^2)
        
        # energy fluence 
        ef[i], eff[i] = energyfluence(E_sum), energyfluence(E_sum_f)

        ### EvxB, EvxvxB, Ev components
        # Create object
        trans = GroundtoShowerCoordinates(antx[i], anty[i], (theta, phi), B)

        E_vxB, E_vxvxB, E_v = [], [], []
        for ei in range (len(Ex)):
            EvxB, EvxvxB, Ev = trans.EtoShower(Ex[ei], Ey[ei], Ez[ei])
            E_vxB.append(EvxB)
            E_vxvxB.append(EvxvxB)
            E_v.append(Ev)
        E_vxB = np.array(E_vxB)
        E_vxvxB = np.array(E_vxvxB)
        E_v = np.array(E_v)

        #filter electric field using signal.butter
        E_vxB_f, E_vxvxB_f, E_v_f = (signal.filtfilt(b, a, i) for i in (E_vxB, E_vxvxB, E_v))
        #vxB component
        E_vxB_2, E_vxB_2_f = magE2(E_vxB, 0, 0), magE2(E_vxB_f, 0, 0)  # |E|^2 raw and filtered
        E_vxB_sum, E_vxB_sum_f = sum(E_vxB_2), sum(E_vxB_2_f) # sum(|E|^2)
        #vxvxB component
        E_vxvxB_2, E_vxvxB_2_f = magE2(E_vxvxB, 0, 0), magE2(E_vxvxB_f, 0, 0)  # |E|^2 raw and filtered
        E_vxvxB_sum, E_vxvxB_sum_f = sum(E_vxvxB_2), sum(E_vxvxB_2_f) # sum(|E|^2)

        # energy fluence 
        ef_vxB[i], eff_vxB[i] = energyfluence(E_vxB_sum), energyfluence(E_vxB_sum_f)
        ef_vxvxB[i], eff_vxvxB[i] = energyfluence(E_vxvxB_sum), energyfluence(E_vxvxB_sum_f)
    
    radE_vals = RadEnergy(finp, Nant, ef, eff, antx, anty)
    radE_vals_vxB = RadEnergy(finp, Nant, ef_vxB, eff_vxB, antx, anty)
    radE_vals_vxvxB = RadEnergy(finp, Nant, ef_vxvxB, eff_vxvxB, antx, anty)
    
    # store data in dictionary
    all_data["run"].append(run)
    all_data["radE(eV)"].append(radE_vals[0])
    all_data["radE_filtered(eV)"].append(radE_vals[1])
    all_data["ef"].append(ef)
    all_data["eff"].append(eff)
    all_data["antx"].append(antx)
    all_data["anty"].append(anty)
    all_data["radE_vxB(eV)"].append(radE_vals_vxB[0])
    all_data["radE_filtered_vxB(eV)"].append(radE_vals_vxB[1])
    all_data["radE_vxvxB(eV)"].append(radE_vals_vxvxB[0])
    all_data["radE_filtered_vxvxB(eV)"].append(radE_vals_vxvxB[1])


final_output = {}
for key, value in all_data.items():
    final_output[key] = np.array(value, dtype=object)

output_dir = os.path.dirname(output_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# save everything in a compressed npz file
np.savez_compressed(output_path, **final_output)
print(f"Saved all data to {output_path}")
