#!/us/bin/env python
# ------------------------------------------------------------------------
import os
import logging
if not 'I3_BUILD' in os.environ:
    raise Exception('To run this script start an IceTray environment')

from icecube import dataio, dataclasses, icetray, recclasses
from icecube.icetray import I3Units
from icecube import icetray
import numpy as np
import pylab, math
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
    
def build_response_table(log_ekin, log_edep, n_bins=50): # find median to interpolate
    edges = np.linspace(log_ekin.min(), log_ekin.max(), n_bins + 1)
    centers, medians, p16, p84 = [], [], [], []

    for i in range(len(edges) - 1):
        mask = (log_ekin >= edges[i]) & (log_ekin < edges[i+1])
        if mask.sum() < 5:
            continue
        vals = log_edep[mask]
        centers.append(0.5 * (edges[i] + edges[i+1]))
        medians.append(np.median(vals))
        p16.append(np.percentile(vals, 16))
        p84.append(np.percentile(vals, 84))

    return np.array(centers), np.array(medians), np.array(p16), np.array(p84)

parser = argparse.ArgumentParser()

mctree  = 'IceTopMCTree'
energy  = 'IceScintHitSeriesMap'

plt.rcParams["font.size"] = 14
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
})

Edep_all = []
energy_all = []

parser.add_argument("--particle", type=str, default="e-", help='particle: e-/gamma/mu+/mumin')
parser.add_argument("--plot", action='store_true', help='--plot to plot')
args   = parser.parse_args()
#infile = dataio.I3File("/data/user/wkammeem/CORSIKA/detector-response/" + args.results)

energy_range = ['0.0001-0.001', '0.001-0.01', '0.01-0.1', '0.1-1.0', '1.0-10.0']

particle_map = {
    'e-':    (dataclasses.I3Particle.EMinus,   r'$e^-$'),
    'mu+':   (dataclasses.I3Particle.MuPlus,   r'$\mu^+$'),
    'gamma': (dataclasses.I3Particle.Gamma,    r'$\gamma$'),
    'mumin': (dataclasses.I3Particle.MuMinus,  r'$\mu^-$'),
}
if args.particle not in particle_map:
    raise ValueError(f"Unknown particle: {args.particle}")
particleclass, label = particle_map[args.particle]

print(rf"processing {args.particle} deposition")

total_q_frames = 0
q_frames_no_edep = 0

for e in energy_range:

    infile = dataio.I3File("/data/user/wkammeem/CORSIKA/detector-response/" + f'scint_response_{args.particle}_{e}GeV.i3')
    for f in infile:
       
        if f.Stop == icetray.I3Frame.DAQ:

            total_q_frames += 1

            if (energy in f):
                energy_scint = f[energy]
                if len(energy_scint) == 0: # if no deposited energy
                    q_frames_no_edep += 1
                else: # if there's deposited energy
                    energy_deposit = [energy_scint[key][0].charge for key in energy_scint] 
                    MCTree = f[mctree]
                    for primary in MCTree:
                        if primary.type == particleclass: 
                            Nscint_photon = energy_deposit[0]
                            Edep = Nscint_photon/(8960) # in MeV
                            Edep_all.append(Edep) # deposited energy
                            energy_all.append(np.log10(primary.energy)) #kinetic energy (log GeV)
                
            else:
                logging.warning('Pulse list was empty')

print(rf'Total Q-frame: {total_q_frames}')
print(rf'No deposited energy Q-Frame: {q_frames_no_edep}')

# Filter: only keep entries where both energy and deposition are positive
energy_all = np.array(energy_all, dtype=float)
Edep_all   = np.array(Edep_all,   dtype=float)

# find mean value
centers, medians, p16, p84 = build_response_table(energy_all, Edep_all)

# interpolate median value
median_interp = interp1d(centers, medians, kind='linear', fill_value='extrapolate')
width_interp  = interp1d(centers, (p84 - p16) / 2, kind='linear', fill_value='extrapolate')

x_smooth = np.linspace(centers.min(), energy_all.max(), 300)
if args.plot:
    plt.scatter(energy_all, Edep_all, s = 1, alpha= 0.3, color = 'darkcyan')
    # plt.plot(x_smooth, median_interp(x_smooth), color='orange', lw=2, label='Median (interpolated)')
    # plt.scatter(centers, medians, facecolors='none', edgecolors='orange', s = 20)
    plt.yscale('log')
    plt.xlabel(rf'$\mathrm{{log_{{10}}(E_{{kin}}/GeV)}}$')
    plt.ylabel(rf'$\mathrm{{N_{{scint \ photons}}}}$')
    plt.minorticks_on()
    plt.title(f'Scintillator Response to GEANT4 {label}')
    # plt.legend()
    plt.savefig(f'figures/ScintResponse/{args.particle}.png', bbox_inches="tight", dpi = 300)