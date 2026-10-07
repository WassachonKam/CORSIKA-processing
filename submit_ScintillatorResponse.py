#!/usr/bin/env python3
# Usage: python3 submit_ScintillatorResponse.py -p proton -e lgE_16.2 -z 0.0 

import os
import logging
if not 'I3_BUILD' in os.environ:
    raise Exception('To run this script start an IceTray environment')

from icecube import dataio, dataclasses, icetray, recclasses
from icecube.icetray import I3Units
from icecube import icetray
import argparse
import numpy as np
import matplotlib.pyplot as plt
from corsikaio import CorsikaParticleFile
import pandas as pd
import random
import os
from scipy.interpolate import interp1d
 

############################# FUNCTIONS AND CONSTANTS #############################

# kinetic energy 
def Ekin(px,py,pz,m):
    return np.sqrt(px**2 + py**2 + pz**2 + m**2) -m 

def R(x,y):
    return np.sqrt(x**2 + y**2)

# constants
mumass = 0.105658 #muon mass in GeV/c^2
emass = 5.11e-4 #electron mass in GeV/c^2
gammamass = .0

class ScintillatorResponse:
    def __init__(self, log_ekin_geant4, log_edep_geant4, k_neighbors=50, p_interact=0.03):
        self.k = k_neighbors
        self.p_interact = p_interact

        sort_idx = np.argsort(log_ekin_geant4)
        self.log_ekin_sorted = log_ekin_geant4[sort_idx]
        self.log_edep_sorted = log_edep_geant4[sort_idx]

    def sample(self, log_ekin_corsika, rng=None):
        if rng is None:
            rng = np.random.default_rng()

        scalar = np.ndim(log_ekin_corsika) == 0
        log_ekin_corsika = np.atleast_1d(np.asarray(log_ekin_corsika, dtype=float))
        result = np.zeros(len(log_ekin_corsika))

        log_ekin_clamped = np.clip(log_ekin_corsika,
                                    self.log_ekin_sorted[0],
                                    self.log_ekin_sorted[-1])
        indices = np.searchsorted(self.log_ekin_sorted, log_ekin_clamped)

        for i, idx in enumerate(indices):
            # Step 1: energy-dependent or flat interaction probability
            p = float(np.clip(self.p_interact(log_ekin_clamped[i]), 0.0, 1.0)) \
                if callable(self.p_interact) else self.p_interact
            if rng.random() > p:
                continue

            # Step 2: k-NN sampling
            lo = max(0, idx - self.k // 2)
            hi = min(len(self.log_ekin_sorted), lo + self.k)
            lo = max(0, hi - self.k)
            j = rng.integers(lo, hi)
            result[i] = 10 ** self.log_edep_sorted[j]

        return float(result[0]) if scalar else result

def compute_interaction_probability(log_ekin_all, edep_all, bin_width=0.1):
    """
    log_ekin_all : log10(E_kin/GeV) for ALL simulated gammas (including non-interacting)
    edep_all     : E_dep/MeV for all gammas (0 for non-interacting)
    bin_width    : bin width in log10(E_kin) — Stef says 0.1
    """
    ekin_min = np.floor(log_ekin_all.min() / bin_width) * bin_width
    ekin_max = np.ceil(log_ekin_all.max()  / bin_width) * bin_width
    edges = np.arange(ekin_min, ekin_max + bin_width, bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])

    p_interact = np.zeros(len(centers))
    p_err      = np.zeros(len(centers))

    for i in range(len(centers)):
        mask = (log_ekin_all >= edges[i]) & (log_ekin_all < edges[i + 1])
        n_total    = mask.sum()
        n_interact = (edep_all[mask] > 0).sum()

        if n_total > 0:
            p = n_interact / n_total
            p_interact[i] = p
            p_err[i] = np.sqrt(p * (1 - p) / n_total)  # binomial uncertainty

    return centers, p_interact, p_err


############################# DATA PREPARATION #############################
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
    # print(event.particles.dtype.names)
    particle_id = event.particles['particle_description'] // 1000 # corsika particle ID
    x = event.particles['x'] # x coordinate
    y = event.particles['y'] # y coordinate
    px = event.particles['px'] # momentum component in x direction in GeV/c
    py = event.particles['py'] # momentum component in y direction in GeV/c
    pz = event.particles['pz'] # momentum component in z direction in GeV/c
    weight = event.particles['thinning_weight'] # particle weight 


    # get indices of particles
    particles = {
        'gamma':  (1, gammamass),
        'e+': (2, emass),
        'e-':   (3, emass),
        'mu+': (5, mumass),
        'mu-': (6, mumass)
    }

    Ek_data, weight_data = {}, {}

    for name, (pid, mass) in particles.items():
        idx = particle_id == pid    # index which particle == particle ID
        mask = R(x[idx], y[idx]) > minR * 1e3  # minimum R in cm
        Ek_data[name]     = Ekin(px[idx], py[idx], pz[idx], mass)[mask]
        weight_data[name] = weight[idx][mask]

    Ek_gamma, Ek_e_plus, Ek_e_min, Ek_mu_plus, Ek_mu_min = Ek_data['gamma'], Ek_data['e+'], Ek_data['e-'], Ek_data['mu+'], Ek_data['mu-']
    weight_gamma, weight_e_plus, weight_e_min, weight_mu_plus, weight_mu_min = weight_data['gamma'], weight_data['e+'], weight_data['e-'], weight_data['mu+'], weight_data['mu-']

    # zenith angle for normalization
    theta = np.deg2rad(thetap)
    
    ############################# SCINTILLATOR RESPONSE #############################   

    total_Edep_all = 0 

    energy_range = ['0.0001-0.001', '0.001-0.01', '0.01-0.1', '0.1-1.0', '1.0-10.0']
    particles = ['e-', 'gamma', 'mu+', 'mumin']

    particle_map = {
        'e-':    (dataclasses.I3Particle.EMinus,   r'$e^-$', 1.0, Ek_e_min, weight_e_min),
        'e+':    (dataclasses.I3Particle.EPlus,   r'$e^+$', 1.0, Ek_e_plus, weight_e_plus),
        'mu+':   (dataclasses.I3Particle.MuPlus,   r'$\mu^+$', 1.0, Ek_mu_plus, weight_mu_plus),
        'gamma': (dataclasses.I3Particle.Gamma,    r'$\gamma$', None, Ek_gamma, weight_gamma),
        'mumin': (dataclasses.I3Particle.MuMinus,  r'$\mu^-$', 1.0, Ek_mu_min, weight_mu_min),
    }


    for particle in particles:

        particleclass, label, p_interact, Ek_corsika, weight_corsika = particle_map[particle]

        total_q_frames = 0
        q_frames_no_edep = 0

        mctree = 'IceTopMCTree'
        energy_key = 'IceScintHitSeriesMap'

        Edep_all = []
        Ek_all = []

        # For energy-dependent P_interact (used for gamma)
        erange_centers     = []
        erange_p_interact  = []
        erange_p_err       = []

        for e in energy_range:
            infile = dataio.I3File("/data/user/wkammeem/CORSIKA/detector-response/" + f'scint_response_{particle}_{e}GeV.i3')

            n_total_file   = 0
            n_interact_file = 0

            for f in infile:
                if f.Stop == icetray.I3Frame.DAQ:
                    total_q_frames += 1
                    n_total_file   += 1

                    if energy_key in f:
                        energy_scint = f[energy_key]
                        if len(energy_scint) == 0:
                            q_frames_no_edep += 1
                            # Zero-dep: just count, no MCTree needed
                        else:
                            n_interact_file += 1
                            MCTree = f[mctree]
                            primary_Ek = None
                            for p in MCTree:
                                if p.type == particleclass:
                                    primary_Ek = p.energy
                                    break
                            if primary_Ek is not None:
                                Nscint_photon = energy_scint[list(energy_scint.keys())[0]][0].charge
                                # Edep = Nscint_photon / 8960
                                Edep = Nscint_photon 
                                Edep_all.append(Edep)
                                Ek_all.append(primary_Ek)
                    else:
                        logging.warning('Pulse list was empty')

            # P_interact for this energy range file
            if n_total_file > 0:
                emin, emax = [float(x) for x in e.split('-')]
                e_center = np.log10(np.sqrt(emin * emax))  # geometric mean in log space
                p = n_interact_file / n_total_file
                erange_centers.append(e_center)
                erange_p_interact.append(p)
                erange_p_err.append(np.sqrt(p * (1 - p) / n_total_file))

        E_kin_GeV = np.array(Ek_all, dtype=float)
        E_dep_MeV = np.array(Edep_all, dtype=float)
        log_ekin  = np.log10(E_kin_GeV)
        log_edep  = np.log10(E_dep_MeV) / np.cos(theta)

        if particle == 'gamma':
            erange_centers    = np.array(erange_centers)
            erange_p_interact = np.array(erange_p_interact)
            erange_p_err      = np.array(erange_p_err)

            p_interact_func = interp1d(erange_centers, erange_p_interact, kind='linear',
                                    bounds_error=False,
                                    fill_value=(erange_p_interact[0], erange_p_interact[-1]))
            p_interact_arg = p_interact_func

        else:
            p_interact_arg = p_interact  # 1.0 for charged particles

        response = ScintillatorResponse(log_ekin, log_edep,
                                        p_interact=p_interact_arg,
                                        k_neighbors=10)

        # Sample particle with Ek_corsika
        log_ekin_sampled = np.log10(Ek_corsika)
        edep_sampled_MeV = response.sample(log_ekin_sampled)

        # Weight deposited energy and calculate the total amount
        weightedEdep = edep_sampled_MeV * weight_corsika # apply weight factor 
        total_Edep = sum(weightedEdep) # total deposited energy of every particles shared the same type
        total_Edep_all += total_Edep # total deposited energy of every particles types

        if particle == 'e-': all_data["Edep_e"].append(total_Edep) 
        elif particle == 'mu+': all_data["Edep_mu"].append(total_Edep)

    all_data["Edep_tot"].append(total_Edep_all)
    print("run: ", run, "Total Nscint: ", total_Edep_all)

output_dir = f"/data/user/wkammeem/CORSIKA/TotalEdepScint/{minR}m"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_path = f"{output_dir}/{primary}_{energy}_{sin2theta}.npz"
np.savez_compressed(output_path, **all_data)
print(f"Saved files to {output_path}")


        
