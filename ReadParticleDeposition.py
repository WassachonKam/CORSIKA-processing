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
args   = parser.parse_args()
#infile = dataio.I3File("/data/user/wkammeem/CORSIKA/detector-response/" + args.results)

energy_range = ['0.0001-0.001', '0.001-0.01', '0.01-0.1', '0.1-1.0', '1.0-10.0']

if args.particle == 'e-': 
    particleclass = dataclasses.I3Particle.EMinus
    label = rf'$e^-$'
if args.particle == 'mu+': 
    particleclass = dataclasses.I3Particle.MuPlus
    label = rf'$\mu^+$'
if args.particle == 'gamma': 
    particleclass = dataclasses.I3Particle.Gamma
    label = rf'$\gamma$'
if args.particle == 'mumin': 
    particleclass = dataclasses.I3Particle.MuMinus
    label = rf'$\mu^-$'

for i in energy_range:
    infile = dataio.I3File("/data/user/wkammeem/CORSIKA/detector-response/" + f'scint_response_{args.particle}_{i}GeV.i3')
    for f in infile:
        if (energy in f):
            energy_scint = f[energy]
            if len(energy_scint)>0:
                energy_deposit = [energy_scint[key][0].charge for key in energy_scint] 
                MCTree = f[mctree]
                for i in MCTree:
                    if i.type==particleclass: 
                        Edep_all.append(energy_deposit[0])
                        energy_all.append(np.log10(i.energy))
            
        else:
            logging.warning('Pulse list was empty')

plt.scatter(energy_all, Edep_all, s = 1, alpha= 0.3, color = 'darkcyan')
plt.yscale('log')
plt.xlabel(rf'$\mathrm{{log_{{10}}(E_{{kin}}/GeV)}}$')
plt.ylabel(rf'$\mathrm{{E_{{deposited}}/MeV}}$')
plt.minorticks_on()
plt.title(f'Scintillator Response of {label}')
plt.savefig(f'figures/ScintResponse/{args.particle}.png', bbox_inches="tight", dpi = 300)