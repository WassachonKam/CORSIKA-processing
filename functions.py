import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import integrate
from scipy import signal, fft, constants, optimize
from corsikaio import CorsikaParticleFile
import random

#=====================================
# set primary particle parameters
#=====================================

primary = "proton"
energy = "lgE_16.0"
sin2theta = "0.0"
runnum = 0

#=====================================
# file paths 
#=====================================
# need to set path to simulation raw data
rawdata = "/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern" 

fp_long = lambda primary, energy, sin2theta, runnum: f'{rawdata}/{primary}/{energy}/sin2_{sin2theta}/{runnum:06d}/DAT{runnum:06d}.long'
fp_list = lambda primary, energy, sin2theta, runnum: f'{rawdata}/{primary}/{energy}/sin2_{sin2theta}/{runnum:06d}/SIM{runnum:06d}.list'
fp_inp = lambda primary, energy, sin2theta, runnum: f'{rawdata}/{primary}/{energy}/sin2_{sin2theta}/{runnum:06d}/SIM{runnum:06d}.inp'
fp_radio = lambda primary, energy, sin2theta, runnum, ant: f'{rawdata}/{primary}/{energy}/sin2_{sin2theta}/{runnum:06d}/SIM{runnum:06d}_coreas/raw_ant_{ant}.dat'
fp_Nmu = lambda primary, energy, sin2theta, runnum: f'Particles/{primary}/{energy}/sin2_{sin2theta}/DAT{runnum:06d}_mupm.npz'
fp_Ne = lambda primary, energy, sin2theta, runnum: f'Particles/{primary}/{energy}/sin2_{sin2theta}/DAT{runnum:06d}_epm.npz'
fp_Nelectron = lambda primary, energy, sin2theta, runnum: f'Particles/{primary}/{energy}/sin2_{sin2theta}/DAT{runnum:06d}_electron.npz'
fp_groundTot = lambda primary, energy, sin2theta: f'GroundTotalParticles/{primary}_{energy}_{sin2theta}.npz'
# fp_Nmu_tot = lambda primary, energy, sin2theta: f'Particles/{primary}/{energy}/sin2_{sin2theta}/TOTAL_mupm.dat'
# fp_Ne_tot = lambda primary, energy, sin2theta: f'Particles/{primary}/{energy}/sin2_{sin2theta}/TOTAL_epm.dat'
fp_RadE = lambda primary, energy, sin2theta: f'radEnergy/raw/{primary}_{energy}_{sin2theta}.npz'
# fp_RadE_norm1 = lambda primary, energy, sin2theta: f'radEnergy/norm_sintheta/{primary}_{energy}_{sin2theta}.dat'
fp_RadE_norm2 = lambda primary, energy, sin2theta: f'radEnergy/norm_sintheta2/{primary}_{energy}_{sin2theta}.npz'
fp_Xmax = lambda primary, energy, sin2theta: f'Xmax/{primary}_{energy}_sin2_{sin2theta}.dat'
fp_DAT = lambda primary, energy, sin2theta, runnum: f'{rawdata}/{primary}/{energy}/sin2_{sin2theta}/{runnum:06d}/DAT{runnum:06d}'
#=====================================
# constants / math functions
#=====================================

# constants
e0 = 8.854e-12 # vacuum permittivity constant (F/m)
c = 2.99e8 # speed of light in vacuum (m/s)
bin_width = 2e-10 #Time resolution
mumass = 0.105658 #muon mass in GeV/c^2
emass = 5.11e-4 #electron mass in GeV/c^2

# functions

# kinetic energy 
def Ekin(px,py,pz,m):
    return np.sqrt(px**2 + py**2 + pz**2 + m**2) -m 

# energy fluence
def energyfluence(sumE2):
    return e0 * c * bin_width* sumE2 * 6.242e+18 # energy fluence in eV 
    
# |E|^2
def magE2(Ex, Ey, Ez):
    return np.abs(np.sqrt(Ex**2 + Ey**2 + Ez**2))**2

#list multiplication
def mullist(r,f):
    return [a * b for a, b in zip(r, f)]

# sum(|Et|^2) = (1/N) sum(|Ef|^2)
def fftsum(N, e_fft2):
    return (1 / N) * np.sum(e_fft2) 

# coordinate conversion
class GroundtoShowerCoordinates:

    def __init__(self, x, y, angle , B):
        self.x = x
        self.y = y
        self.theta, self.phi = angle
        self.Bx, self.By, self.Bz = B

        # transformation 
        self.v = np.array([
            +np.sin(self.theta)*np.cos(self.phi),
            +np.sin(self.theta)*np.sin(self.phi),
            -np.cos(self.theta) ])
        
        # normalized B field
        self.B = self.normalized(np.array([self.Bx, self.By, self.Bz]))

        # basis vectors
        self.e1 = self.normalized(np.cross(self.v, self.B))
        self.e2 = self.normalized(np.cross(self.v, self.e1))
        
    def normalized(self, v):
        return v / np.linalg.norm(v)

    def vxB(self):
        dr = np.array([self.x, self.y, 0])
        return np.dot(dr, self.e1)

    def vxvxB(self):
        dr = np.array([self.x, self.y, 0])
        return np.dot(dr, self.e2)
    
    def EtoShower(self, Ex, Ey, Ez): # convert Ex, Ey, Ez to EvxB, Evx(vxB)
        E = np.array([Ex, Ey, Ez])
        E_vxB = np.dot(E, self.e1)
        E_vxvxB = np.dot(E, self.e2)
        E_v = np.dot(E, self.v)  # optional longitudinal component

        return E_vxB, E_vxvxB, E_v

# Radiation energy
def RadEnergy(inpfile, Nant, ef, eff, ant_x, ant_y):

    vxB, vxvxB,  = (np.empty(Nant) for _ in range(2))
    ### shower coordinate ###
    with open(inpfile) as f:
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
    
    # calculation loop
    for i in range(Nant):
        g2s = GroundtoShowerCoordinates(ant_x[i], ant_y[i], [theta, phi], [Bx, By, Bz])
        new_x = g2s.vxB()
        new_y = g2s.vxvxB()
        vxB[i] = new_x
        vxvxB[i] = new_y
    
    
    vxB = np.asarray(vxB)
    vxvxB = np.asarray(vxvxB)
    ef = np.asarray(ef)
    eff = np.asarray(eff)

    
    # select data on phi = 90 
    angle = np.arctan2(vxvxB, vxB)
    mask  = (np.round(angle,2) == np.round(np.pi/2,2)) & (vxvxB >= 0) 
    # mask = (np.abs(vxB) < 1e-1) & (vxvxB >= 0) #select x and y position only phi = 90
    r = vxvxB[mask]
    f1 = ef[mask]
    f2 = eff[mask]
    
    # r*f
    rf1 = r * f1
    rf2 = r * f2
    
    # trapezoidal integration
    int_rf1 = integrate.trapezoid(rf1,r)
    int_rf2 = integrate.trapezoid(rf2,r)
    
    # radiation energy calculation
    rad_E1 = 2*np.pi*int_rf1
    rad_E2 = 2*np.pi*int_rf2

    return rad_E1, rad_E2, f1, f2 #raw and filtered
    
#=====================================
# plotting functions 
#=====================================

# text boxes for plots 
class textboxes:

    def __init__(self, fig, E_sum, E_sum_f):
        self.E_sum = E_sum
        self.E_sum_f = E_sum_f
        self.fig = fig

    def allbands(self, time0, time1):
        self.fig.text(
            0.95, 0.80,
            rf'full frequency band' + '\n' 
            + rf'bin width = {time1-time0:.2e}' + '\n'
            + rf'$\Sigma |\mathrm{{E_t}}|^2$ = {self.E_sum:.2e} '
              r'$\mathrm{V^2\,m^{-2}}$' + '\n'
            + f'energy fluence = {energyfluence(self.E_sum):.2e} '
              r'$\mathrm{eV\,m^{-2}}$',
            ha='left',
            va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=11
        )

    def filtered(self):
        self.fig.text(
            0.95, 0.60,
            rf'70-350 MHz' + '\n' 
            + rf'$\Sigma |\mathrm{{E_t}}|^2$ = {self.E_sum_f:.2e} '
              r'$\mathrm{V^2\,m^{-2}}$' + '\n'
            + f'energy fluence = {energyfluence(self.E_sum_f):.2e} '
              r'$\mathrm{eV\,m^{-2}}$',
            ha='left',
            va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=11
        )
        
# plot particle number
def pltNpar(parx, pary, parw, ptype): # number of particles on x, and y axes with weight and label ('muon', 'electron')

    if ptype == 'muon': label = '\mu'
    elif ptype == 'electron': label = 'e'
        
    hist = 'step'
    fig = plt.figure(figsize=(7,7))
    gs = gridspec.GridSpec(2, 2, width_ratios=[2.5,1], height_ratios=[1,2.5],
                           wspace=0.05, hspace=0.05)
    
    ax1 = fig.add_subplot(gs[1, 0])
    h = ax1.hist2d(parx, pary, bins=50, weights = parw, norm=LogNorm(), label = ptype)
    # ax1.set_aspect('equal')
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.set_xlim(-1000, 1000)
    ax1.set_ylim(-1000, 1000)
    
    ax2 = fig.add_subplot(gs[0, 0], sharex=ax1)
    counts_par, bins_par, _ = ax2.hist(parx, bins=25,  weights = parw, histtype = hist, label = rf'${ptype}^{{\pm}}$')
    ax2.set_ylabel('particle number')
    ax2.set_yscale('log')
    # plt.legend()
    
    ax3 = fig.add_subplot(gs[1, 1], sharey=ax1)
    ax3.hist(pary, bins=25, orientation='horizontal',  weights = parw, histtype = hist, label = rf'${ptype}^{{\pm}}$')
    ax3.set_xlabel('particle number')
    ax3.set_xscale('log')
    # plt.legend()
    
    # Optional: remove ticks on shared axes for cleanliness
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    
    cbar_ax = fig.add_axes([0.12, -0.02, 0.5, 0.02])  # adjust numbers to move/resize
    cbar = fig.colorbar(h[3], cax=cbar_ax, orientation='horizontal')
    cbar.set_label(r'$\mu^{\pm}$ counts')
    
    
    # Text box
    fig.text(
        0.79, 0.79,                   
        rf'Total ${label}^{{\pm}}$' + '\n' + f'{counts_par.sum():.2e}',
        ha='center',
        va='center',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )
    
    plt.show()
    
    
# plot radius vs energy fluence
def pltef(vxB, vxvxB, ef, eff, method, xmin,xmax, scale, primary, energy, sin2theta):
    vxB = np.asarray(vxB)
    vxvxB = np.asarray(vxvxB)
    ef = np.asarray(ef)
    eff = np.asarray(eff)

    
    # select data on phi = 90 
    theta = np.arctan2(vxvxB, vxB)
    mask  = (np.round(theta,2) == np.round(np.pi/2,2)) & (vxvxB >= 0) 

    r = vxvxB[mask]
    f1 = ef[mask]
    f2 = eff[mask]
    
    # r*f
    rf1 = r * f1
    rf2 = r * f2
    
    # trapezoidal integration
    int_rf1 = integrate.trapezoid(rf1,r)
    int_rf2 = integrate.trapezoid(rf2,r)
    
    # radiation energy calculation
    rad_E1 = 2*np.pi*int_rf1
    rad_E2 = 2*np.pi*int_rf2

    ############# select data on phi = 90 and 270
    mask2  = (np.round(theta,2) == np.round(np.pi/2,2)) | (np.round(theta,2) == np.round(-np.pi/2,2))
    rr = vxvxB[mask2]
    ff1 = ef[mask2]
    ff2 = eff[mask2]
    
    # r*f
    rrf1 = rr * ff1
    rrf2 = rr * ff2
    
    # trapezoidal integration
    int_rrf1 = integrate.trapezoid(rrf1,rr)
    int_rrf2 = integrate.trapezoid(rrf2,rr)
    
    # radiation energy calculation
    rrad_E1 = np.pi*int_rrf1
    rrad_E2 = np.pi*int_rrf2

    fig = plt.figure()
    plt.rcParams["font.size"] = 12
    plt.scatter(r,f2, label = '70-350 Mz ($\phi = 90$)', color = 'red')
    plt.scatter(r,f1, label = r'full frequency band ($\phi = 90$)', color = 'blue')
    plt.xlabel("distance to shower axis (m)")
    plt.ylabel(rf"energy fluence  ($\mathrm{{eV}} \cdot \mathrm{{m}}^{{-2}}$)")
    if method == 'bp':
        plt.title(rf"Bandpass Filter $\phi = 90$ ({primary}, {energy}, sin2_{sin2theta})")
    if method == 'fft':
        plt.title(rf"FFT $\phi = 90$ ({primary}, {energy}, sin2_{sin2theta})")
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.text(
        1.05, 0.75,
        rf'$E_{{\rm rad, full \ bands }}= {rad_E1:.2e}$ eV' + '\n'
        rf'$E_{{\rm rad, 70-350 MHz }} = {rad_E2:.2e}$ eV',
        transform=plt.gca().transAxes,
        va='top'
    )
    plt.xscale(scale)
    # plt.xlim(xmin,xmax)
    plt.show()

    fig = plt.figure()
    plt.rcParams["font.size"] = 12
    plt.scatter(rr,ff1,label = 'full frequency band ($\phi = 90, 270$)', color = 'blue')
    plt.scatter(rr,ff2,label = '70-350 Mz ($\phi = 90, 270$)', color = 'red')
    plt.xlabel("distance to shower axis (m)")
    plt.ylabel(rf"energy fluence  ($\mathrm{{eV}} \cdot \mathrm{{m}}^{{-2}}$)")
    if method == 'bp':
        plt.title(rf"Bandpass Filter $\phi = 90,270$  ({primary}, {energy}, sin2_{sin2theta})")
    if method == 'fft':
        plt.title(rf"FFT $\phi = 90,270$  ({primary}, {energy}, sin2_{sin2theta})")
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.text(
        1.05, 0.75,
        rf'$E_{{\rm rad, full \ bands }} = {rrad_E1:.2e}$ eV' + '\n'
        rf'$E_{{\rm rad, 70-350 MHz }} = {rrad_E2:.2e}$ eV',
        transform=plt.gca().transAxes,
        va='top'
    )
    plt.xscale(scale)
    # plt.xlim(xmin,xmax)
    plt.show()

    return len(ff1), len(ff2)


# plot energy fluence color map
def pltefmap(finp, Nant, vxB, vxvxB, ant_x, ant_y, colors):
    
    ### ground coordinate ###
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.3)
    plot = axes[0].scatter(ant_x, ant_y, c=colors, s = 25, cmap = 'jet')
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    axes[0].set_title("Ground Coordinates")
    
    
    ### shower coordinate ###
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
    
    # calculation loop
    for i in range(Nant):
        g2s = GroundtoShowerCoordinates(ant_x[i], ant_y[i], [theta, phi], [Bx, By, Bz])
        new_x = g2s.vxB()
        new_y = g2s.vxvxB()
        vxB[i] = new_x
        vxvxB[i] = new_y
    
    # plot
    plot = axes[1].scatter(vxB, vxvxB, c=colors, s = 25, cmap = 'jet')
    axes[1].set_xlabel(rf"$\hat{{v}}\times\hat{{B}} \ \mathrm{{direction \ (m)}}$")
    axes[1].set_ylabel(rf"$\hat{{v}} \times (\hat{{v}}\times\hat{{B}}) \ \mathrm{{direction \ (m)}}$")
    axes[1].set_title("Shower Coordinates")
    # axes[1].set_xlim(-0.05,0.05)
    # axes[1].set_ylim(-0.05,0.05)
    

    cbar_ax = inset_axes(axes[1],
                         width="5%", # width = 5% of parent axes width
                         height="100%", # height = 100% of parent axes height
                         loc='right', # fixed location
                         borderpad=-3 # padding between the axes and colorbar
                        )
    cbar = fig.colorbar(plot, cax=cbar_ax)
    cbar.set_label(rf'energy fluence $(\mathrm{{eV}} \cdot \mathrm{{m}}^{{-2}})$', rotation=-90, labelpad=30)
    # plt.suptitle(rf'{primary}, {energy}, sin2_{sin2theta}, run {runnum}')
    plt.show()


# plot |E|^2 histogram
def pltEmag2(ant_no, timer, timef, Emag2, Emag2_f, E_sum, E_sum_f, style):

    
    fig = plt.figure(figsize=(8,5))
    plt.title(f'ant {ant_no}')
    
    if style == 'bp':
        plt.plot(timer, Emag2, label='full frequency band', color = 'blue')
        plt.plot(timef, Emag2_f, label='70–350 MHz', color = 'red')
        plt.xlim(0,50)
        plt.xlabel("time (ns)")
    elif style == 'fft':
        plt.scatter(timer, Emag2, label='full frequency band', color = 'blue', s = 5)
        plt.scatter(timef, Emag2_f, label='70–350 MHz', color = 'red', s = 5)
        plt.xlabel("Frequency (Hz)")
    
    
    plt.ylabel(r"$|\mathrm{E_t}|^2\ (\mathrm{V^2\,m^{-2}})$")
    plt.title(rf'{primary}, {energy}, sin2_{sin2theta}, run {runnum}')
    plt.legend()

    tb = textboxes(fig, E_sum, E_sum_f)
    tb.allbands(timer[0], timer[1])
    tb.filtered()

    plt.show()

# plt E
def pltE(ant_no, time, Ex,Ey, Ez):
        plt.figure(figsize=(8,5))
        plt.title(f'ant {ant_no}')
        plt.plot(time,Ex, label = rf'$E_x$')
        plt.plot(time,Ey, label = rf'$E_y$')
        plt.plot(time,Ez, label = rf'$E_z$')
        plt.xlim(0,50)
        plt.xlabel("time (ns)")
        plt.ylabel("E (V/m)")
        plt.title(rf'{primary}, {energy}, sin2_{sin2theta}, run {runnum}')
        plt.legend()
    
# plot longitudinal profile
def pltlp(atmdepth, positron, electron, muplus, muminus, tot_e, tot_mu, Xmax, Ne_Xmax, RadE, runnum, Xv):
    color_mu =  'steelblue'
    color_e = 'firebrick'
    
    fig, axes = plt.subplots(1, 3, sharey=True)
    fig.subplots_adjust(wspace=0.1)  

    Xmax = Xmax
    Ne_Xmax = Ne_Xmax

    
    axes[0].plot(positron, atmdepth , label = r'$e^+$', color = color_e, ls ='--')
    axes[0].plot(electron, atmdepth,  label = r'$e^-$', color = color_e, ls =':')
    axes[0].legend()
    axes[0].set_ylabel(r"atmospheric depth (g/$\mathrm{cm}^2$)")
    
    axes[1].plot(muplus, atmdepth, label = r'$\mu^+$', color = color_mu , ls ='--')
    axes[1].plot(muminus, atmdepth, label = r'$\mu^-$',  color = color_mu, ls = ':')
    axes[1].legend()
    axes[1].set_xlabel("particle number")
    # axes[1].set_title(f'{primary}, {energy}, sin2theta = {sin2theta}, run {runnum}')
    
    axes[2].plot(tot_mu, atmdepth, label = r'$\mu^{\pm} (\times 50)$', color = color_mu)
    axes[2].plot(tot_e, atmdepth, label = r'$e^{\pm}$', color = color_e)
    axes[2].hlines(Xmax, min(tot_e), max(tot_e), color = 'black')
    axes[2].hlines(Xv, min(tot_e), max(tot_e), color='black', linestyle='--') 
    axes[2].text(0, Xmax, 'Xmax', ha='left', va='bottom')
    axes[2].text(0, Xv, 'ground', ha='left', va='bottom')
    axes[2].legend()

    plt.gca().invert_yaxis()

    # Text box
    fig.text(
        0.95, 0.8, 
        rf'primary: {primary}' + '\n' +
        rf'energy: {energy}' + '\n' +
        rf'$sin^2 \theta = {sin2theta}$' + '\n' +
        rf'run {runnum}' + '\n' +        
        rf'$N_{{e,Xmax}}$ = {Ne_Xmax:.2e} ' + '\n' +
        rf'radiation energy = {RadE:.2e} eV/$sin^2 \alpha$',
        ha='left',
        va='top',
        #bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )

    
   

    save_dir = rf'figure/LongitudinalProfile/{primary}_{energy}_{sin2theta}'
    os.makedirs(save_dir, exist_ok=True)
    # plt.savefig(f'{save_dir}/{runnum}.jpg', bbox_inches="tight", dpi = 300)

    plt.show()

# plot correlation of Nmu and Ne
def pltmuecorr(energydir, sin2theta, energylabel):
    protoncolor = 'red'
    ironcolor = 'blue'
    scattersize = 5
    
    protonpath = f'GroundTotalParticles/proton_lgE_{energydir}_{sin2theta}.npz'
    protondata = np.load(protonpath)
    ptotmu = protondata['nMu'] # total +-muon at ground
    ptote = protondata['nEP']  # total +- e at ground
    
    plt.scatter(ptote, ptotmu, color = protoncolor, s = scattersize, label = 'proton')
    
    ironpath = f'GroundTotalParticles/iron_lgE_{energydir}_{sin2theta}.npz'
    irondata = np.load(ironpath)
    fetotmu = irondata['nMu'] # total +-muon at ground
    fetote = irondata['nEP'] # total +- e at ground
    
    plt.scatter(fetote, fetotmu, color = ironcolor, s = scattersize, label = 'iron')
    plt.text(
        np.mean(fetote), max(fetotmu) * 1.2,                   
        fr'{energylabel}',
        fontsize=14,
        ha='center',
        va='center',
    )
    if energylabel == "1 PeV": plt.legend()
    plt.xscale("log")
    plt.ylim(2e4,3e7)
    plt.yscale("log")
    
    plt.xlabel("number of electrons")
    plt.ylabel("number of muon")


# plot correlation between Ne_Xmax and RadE. RadE normalized by sin2theta, norm = True 
def pltRadE_NeXmax(sin2thetas, primary, energy, labels, norm, filtering): 
    colors = plt.cm.viridis(np.linspace(0, 1, len(sin2thetas)))
    for p in primary:
        for i in range (len(energy)):
            e = energy[i]
            for sin2theta, c in zip(sin2thetas, colors):
                    
                    fNe_tot = fp_Xmax(p, e, sin2theta)

                    if norm == True:
                        fRadE = fp_RadE_norm2(p, e, sin2theta)
                        plt.ylabel(rf'radiation energy (eV) / $sin^2 \alpha$')
                    if norm == False:
                        fRadE = fp_RadE(p, e, sin2theta)
                        plt.ylabel(f'radiation energy (eV)')
                    fileNe = np.loadtxt(fNe_tot)
                    Ne = fileNe[:,3]
                    fileRadE = np.load(fRadE, allow_pickle=True)
  
                    if filtering == True: RadE = fileRadE['radE_filtered(eV)']
                    if filtering == False: RadE = fileRadE['radE(eV)']
                    
                    
                    plt.scatter(Ne, RadE, s = 7, color = c, label = rf"$sin^2\theta$ = {sin2theta}")
            
                    if sin2theta == "0.7":
                        if norm == False: 
                            txty = max(RadE) * 1.5
                            txtx = np.mean(Ne)
                        if norm == True: 
                            txty = max(RadE) * 3
                            txtx = np.mean(Ne) 
                        plt.text(
                            txtx, txty,                   
                            fr'{labels[i]}',
                            fontsize=14,
                            ha='center',
                            va='center',
                            color = 'red'
                            )
                    if e == "lgE_16.0":
                        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
            
                            
                
                    plt.xlabel(rf'$N_e$ at $X_{{max}}$')
                    
                    plt.xscale('log')
                    plt.yscale('log')
                    # plt.ylim(2e4,2e12)
                    plt.title(rf'primary particle: {p}, filtering = {filtering}')
        plt.show()

# plot correlation between Ne ratio vs Xmax, style = 'horiz' or style = 'verti'
def pltNeRatio(energy, style):
    
    sin2thetas = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7",  "0.8", "0.9"]
    primary = ["proton", "iron"]
    X0 = 697.6 #g/cm^2

    if style == 'verti': fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(8, 16), sharex = True, sharey = True)
    elif style == 'horiz': fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(16, 8), sharex = True, sharey = True)
    fig.suptitle(f'{energy[0]}', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    for i in range(len(sin2thetas)):
        sin2theta = sin2thetas[i]

        #determine row and column for plotting
        if style == 'verti':
            if i % 2 == 0:
                r = int(i/2)
                c = 0
            else:
                r = int((i-1)/2)
                c = 1
        elif style == 'horiz':
            if i <= 4:
                r = 0
                c = int(i)
            else:
                r = 1
                c = int(i-5)

        for j in range(len(energy)):
            e = energy[j]
            for k in range(len(primary)):
                p = primary[k]
                fileXmax = np.loadtxt(fp_Xmax(p, e, sin2theta))
                fileNeground = np.load(fp_groundTot(p, e, sin2theta), allow_pickle= True)
                Ne_Xmax = fileXmax[:,5]
                Xmax = fileXmax[:,1]
                Ne_ground = fileNeground['nEP'] # total number of +-e at ground
                
                ratio = Ne_ground/Ne_Xmax

                thetafloat = float(sin2theta)
                costheta = np.cos(np.arcsin(np.sqrt(thetafloat)))

                Xv = X0/ costheta

                if p == "proton": 
                    color = 'red'
                if p == "iron": color = 'blue'

                ax = axs[r, c]
                ax.scatter(Xmax,ratio, c = color, s = 5, label = f"{p}")
                ax.set_title(rf"$sin^2 \theta = ${sin2theta}")
                ax.axvline(x=Xv, color='black', linestyle='--') 
                
                if Xv < 1000: ax.text(Xv +8 , 0.5, 'ground', ha='left', va='bottom')
                if e == "lgE_17.0": plt.legend()

                if style == 'verti':
                    if c == 0: ax.set_ylabel(rf"$N_{{e_{{ground}}}}/N_{{e_{{Xmax}}}}$")
                    if r == 4: ax.set_xlabel(rf"$X_{{max}}$ (g/cm$^2$)")

                if style == 'horiz': 
                    if c == 0: ax.set_ylabel(rf"$N_{{e_{{ground}}}}/N_{{e_{{Xmax}}}}$")
                    if r == 1: ax.set_xlabel(rf"$X_{{max}}$ (g/cm$^2$)")
                
                ax.set_xlim(500,1100)
                ax.set_ylim(- 0.1,1.2)

    plt.subplots_adjust(wspace=0.1)
    plt.legend()
    plt.show()


# plot correlation between Ne vs Xmax vs radiation energy with all zenith angle 
# filtering = True, 70-350 MHz
# filtering = False, full band
# style = 'horiz' or style = 'verti'

def pltNeXmaxRadEAllzenith( primary, energy, filtering, style, corr):
    sin2thetas = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7",  "0.8", "0.9"]
    X0 = 697.6 #g/cm^2

    if style == 'verti': 
        fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(10, 16), sharex = True, sharey = True)
        fig.subplots_adjust(wspace=0.15, hspace = 0.15, top = 0.95)
    elif style == 'horiz': 
        fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 12), sharex = False, sharey = True)
        # fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 12), sharex = False, sharey = False)
        fig.subplots_adjust(wspace=0.05, hspace = 0.1, top = 0.93)


    fig.suptitle(f'{primary[0]}, {energy[0]}, filtering = {filtering}', fontsize=16)

    for i in range(len(sin2thetas)):
        sin2theta = sin2thetas[i]

        #determine row and column for plotting
        if style == 'verti':
            if i % 2 == 0:
                r = int(i/2)
                c = 0
            else:
                r = int((i-1)/2)
                c = 1
        elif style == 'horiz':
            if i <= 4:
                r = 0
                c = int(i)
            else:
                r = 1
                c = int(i-5)

        for j in range(len(energy)):
            e = energy[j]
            for k in range(len(primary)):
                p = primary[k]
                
                fXmax = fp_Xmax(p, e, sin2theta)
                fRadE = fp_RadE_norm2(p, e, sin2theta)
                fRadE_unnorm = fp_RadE(p, e, sin2theta)

                fileXmax = np.loadtxt(fXmax)
                Ne = fileXmax[:,5] # Ne at Xmax
                Xmax = fileXmax[:,6]
                runnum = fileXmax[:,0]
                fileRadE = np.load(fRadE, allow_pickle=True)
                fileRadE_unnorm = np.load(fRadE_unnorm, allow_pickle=True)

        
                if filtering == True: 
                    RadE = fileRadE['radE_filtered(eV)']
                    RadE_unnorm = fileRadE_unnorm['radE_filtered(eV)']
                    RadE_vxB = fileRadE_unnorm['radE_filtered_vxB(eV)']
                    RadE_vxvxB = fileRadE_unnorm['radE_filtered_vxvxB(eV)']
                    RadE_vxB_norm = fileRadE['radE_filtered_vxB(eV)']
                    
                if filtering == False: 
                    RadE = fileRadE['radE(eV)']
                    RadE_unnorm = fileRadE_unnorm['radE(eV)']
                    RadE_vxB = fileRadE_unnorm['radE_vxB(eV)']
                    RadE_vxvxB = fileRadE_unnorm['radE_vxvxB(eV)']
                    RadE_vxB_norm = fileRadE['radE_vxB(eV)']
                

                # angle value
                alpha = fileRadE['alpha']
                sin2alpha = np.sin(alpha)**2
                thetafloat = float(sin2theta)
                theta = np.arcsin(np.sqrt(thetafloat))
                costheta = np.cos(theta)

                Xv = X0/ costheta
                
                # print(len(sin2alpha), len(runnum), len(RadE), len(RadE_unnorm))
                # Exclude data if Xmax is negative value 
                mask = Xmax >= 0
                Xmax = Xmax[mask]
                Ne = Ne[mask]
                runnum = runnum[mask]
                sin2alpha = sin2alpha[mask]

                RadE = RadE[mask]
                RadE_unnorm = RadE_unnorm[mask]
                RadE_vxB = RadE_vxB[mask]
                RadE_vxvxB = RadE_vxvxB[mask]
                RadE_vxB_norm = RadE_vxB_norm[mask]

                ax = axs[r, c]

                if corr == 'NeErad_norm':
                    sc = ax.scatter(Ne, RadE, s=7, c=Xmax, cmap="viridis", vmin = None, vmax = None)
                    # sc = ax.scatter(Ne, RadE, s=7, c=Xmax, cmap="viridis", vmin = 400, vmax = 650) # 400, 600 for iron, 400, 1100 for proton
                    ax.set_yscale('log')
                    # ax.set_xlim(0.3e7,1e7)
                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r" $\mathrm{{E_{{rad}}}}$ (eV) / $\sin^2\alpha$")
                        if r == 4: ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        plt.colorbar(sc, label=rf"$X_{{max}}$ (g/cm$^2$)")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad}}}}$ (eV)/ $\sin^2\alpha$")
                        ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        plt.colorbar(sc, label=rf"$X_{{max}}$ (g/cm$^2$)", orientation='horizontal')
                
                elif corr == 'NeErad':
                    sc = ax.scatter(Ne, RadE_unnorm, s=7, c=Xmax, cmap="viridis", vmin = None, vmax = None)
                    # sc = ax.scatter(Ne, RadE_unnorm, s=7, c=Xmax, cmap="viridis", vmin = 400, vmax = 650)
                    ax.set_yscale('log')
                    # ax.set_xlim(0.3e7,1e7)
                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad}}}}$ (eV)")
                        if r == 4: ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        plt.colorbar(sc, label=rf"$X_{{max}}$ (g/cm$^2$)")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad}}}}$ (eV)")
                        ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        plt.colorbar(sc, label=rf"$X_{{max}}$ (g/cm$^2$)", orientation='horizontal')
                    
                elif corr == 'NeXmax':
                    norm = LogNorm(vmin=RadE.min(), vmax=RadE.max())
                    # norm = LogNorm(vmin=10,vmax = 1e7)
                    sc = ax.scatter(Ne, Xmax, s=7, c=RadE, cmap="viridis", norm=norm)
                    ax.axhline(y=Xv, color='black', linestyle='--') 
                    # ax.set_xlim(0.3e7,1e7)
                    ax.set_ylim(400,1100)
                    if Xv < 1000: ax.text(np.mean(Ne) , Xv +8, 'ground', ha='left', va='bottom')

                    if style == 'verti':
                        if c == 0: ax.set_ylabel(rf"$X_{{max}}$ (g/cm$^2$)")
                        if r == 4: ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        plt.colorbar(sc, label=r"$\mathrm{{E_{{rad}}}}$ (eV)/ $\sin^2\alpha$")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(rf"$X_{{max}}$ (g/cm$^2$)")
                        # if r == 1: ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        ax.set_xlabel(rf"$N_{{e,Xmax}}$")
                        plt.colorbar(sc, label=r"$\mathrm{{E_{{rad}}}}$ (eV)/ $\sin^2\alpha$", orientation='horizontal')

                elif corr == 'EradSina':
                    sc = ax.scatter(sin2alpha, RadE_unnorm, s=7, c=Ne, cmap="viridis")
                    ax.set_yscale('log')

                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad}}}}$ (eV)")
                        if r == 4: ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad}}}}$ (eV)")
                        ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$", orientation='horizontal')

                elif corr == 'EvxBsina':
                    sc = ax.scatter(sin2alpha, RadE_vxB, s=7, c=Ne, cmap="viridis")
                    ax.set_yscale('log')

                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxB}}}}$ (eV)")
                        if r == 4: ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxB}}}}$ (eV)")
                        ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$", orientation='horizontal')

                elif corr == 'EvxvxBNe':
                    sc = ax.scatter(Ne, RadE_vxvxB, s=7, c=sin2alpha, cmap="viridis")
                    ax.set_yscale('log')

                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxvxB}}}}$ (eV)")
                        if r == 4: ax.set_xlabel(rf'$N_{{e,Xmax}}$')
                        plt.colorbar(sc, label=rf"$sin^2 \alpha$")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxvxB}}}}$ (eV)")
                        ax.set_xlabel(rf'$N_{{e,Xmax}}$')
                        plt.colorbar(sc, label=rf"$sin^2 \alpha$", orientation='horizontal')

                    # ax.set_xlim(0.3e7,1e7)
                
                elif corr == 'EvxvxBsina':
                    sc = ax.scatter(sin2alpha, RadE_vxvxB, s=7, c=Ne, cmap="viridis")
                    ax.set_yscale('log')

                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxvxB}}}}$ (eV)")
                        if r == 4: ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxvxB}}}}$ (eV)")
                        ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$", orientation='horizontal')

                elif corr == 'EvxBsina_norm':
                    # print(RadE_vxB_norm)
                    sc = ax.scatter(sin2alpha, RadE_vxB_norm, s=7, c=Ne, cmap="viridis")
                    ax.set_yscale('log')

                    if style == 'verti':
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxB}}}}/ sin^2 \alpha$ (eV)")
                        if r == 4: ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$")

                    elif style == 'horiz': 
                        if c == 0: ax.set_ylabel(r"$\mathrm{{E_{{rad,vxB}}}}/ sin^2 \alpha$ (eV)")
                        ax.set_xlabel(rf'$sin^2 \alpha$')
                        plt.colorbar(sc, label=rf"$N_{{e,Xmax}}$", orientation='horizontal')

                ax.set_title(rf"$sin^2 \theta = ${sin2theta} ({np.rad2deg(theta):.2f}$^\circ$)")
                
    plt.show()
    
# plot correlation between Ne_ground and RadE. RadE normalized by sin2theta, norm = True 
def pltRadE_NeGround(sin2thetas, primary, energy, labels, norm, filtering):
    colors = plt.cm.viridis(np.linspace(0, 1, len(sin2thetas)))
    for p in primary:
        for i in range (len(energy)):
            e = energy[i]
            for sin2theta, c in zip(sin2thetas, colors):

                fNe_tot = fp_groundTot(p, e, sin2theta)
                fRadE = fp_RadE(p, e, sin2theta)
                if norm == True:
                    fRadE = fp_RadE_norm2(p, e, sin2theta)
                    plt.ylabel(rf'radiation energy (eV) / $sin^2 \alpha$')
                if norm == False:
                    fRadE = fp_RadE(p, e, sin2theta)
                    plt.ylabel(f'radiation energy (eV)')
                fileNe = np.load(fNe_tot, allow_pickle=True)
                fileRadE = np.load(fRadE, allow_pickle=True)
                Ne = fileNe["nEP"]
                if filtering == True: RadE = fileRadE['radE_filtered(eV)']
                if filtering == False: RadE = fileRadE['radE(eV)']

                plt.scatter(Ne, RadE, s = 7, color = c, label = rf"$sin^2\theta$ = {sin2theta}")
                
                if sin2theta == "0.7":
                    plt.text(
                        np.mean(Ne), max(RadE) * 1.2,                   
                        fr'{labels[i]}',
                        fontsize=14,
                        ha='center',
                        va='center',
                        )
                if e == "lgE_16.0":
                    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

        
        plt.xlabel(rf'$N_e$ at ground level')
        plt.xscale('log')
        plt.yscale('log')
        plt.title(rf'primary particle: {p}, filtering = {filtering}')
        plt.show()



def pltEdepEkin(primary, energy, sin2theta, run):

    file_path = (f'/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energy}/sin2_{sin2theta}/{run:06d}/DAT{run:06d}')
    file_input = (f"/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energy}/sin2_{sin2theta}/{run:06d}/SIM{run:06d}.inp")


    with open(file_input) as f:
        for line in f:
            parts = line.split()
            if parts[0] == "THETAP":
                thetap = float(parts[1])


    with CorsikaParticleFile(file_path, thinning= True) as file:
        # we only have one event per file, we can grab it like this
        event = next(file)

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
    Ek_mu = Ekin(px[idx_mu], py[idx_mu], pz[idx_mu], mumass) * weight[idx_mu]
    Ek_epm = Ekin(px[idx_epm], py[idx_epm], pz[idx_epm], emass) * weight[idx_epm]
    Ek_electron = Ekin(px[idx_electron], py[idx_electron], pz[idx_electron], emass) * weight[idx_electron]

    # zenith angle for normalization
    theta = np.deg2rad(thetap)

    ############################# SCINTILLATOR RESPONSE #############################   

    # digitized data from Agnieszka's thesis
    dfe = pd.read_csv('ScintillatorResponse/electron_0.0.csv')
    dmu = pd.read_csv('ScintillatorResponse/muon_0.0.csv')
    digitized_files = [dfe, dmu]
    Ek_array = [Ek_electron, Ek_mu]
    plt_title = [rf'$e^{{-}}$', rf'$\mu^{{\pm}}$']
    savefile = ['electron', 'muon']

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

        total_Edep = sum(corsika_y)


        # make a scatter plot
        fig1, ax1 = plt.subplots() 
        fig2, ax2 = plt.subplots()

        # generate randon Gaussian fluctuation
        Edep = interpolate_y
        for ft in range (len(Edep)):
            yfluc = []
            xfluc = []
            for j in range(20):
                fluct = random.gauss(1, 0.5)
                fEdep = Edep[ft]*fluct
                yfluc.append(fEdep)
                xfluc.append(interpolate_x[ft])
            if ft == 0:
                ax1.scatter(xfluc, yfluc, s = 2, alpha = 0.3, color = 'lightskyblue', label = 'Gaussian fluctuation')
            else: ax1.scatter(xfluc, yfluc, s = 2, alpha = 0.3, color = 'lightskyblue')


        interpolate_xcorr = np.linspace(min(corsika_x), max(corsika_x), num = 400)
        interpolate_ycorr = np.interp(interpolate_xcorr, df_x, df_y)
        Edep = interpolate_ycorr

        for ft in range (len(Edep)):
            cor_yfluc = []
            cor_xfluc = []
            for j in range(20):
                fluct = random.gauss(1, 0.5)
                fEdep = Edep[ft]*fluct
                cor_yfluc.append(fEdep)
                cor_xfluc.append(interpolate_xcorr[ft])
            if ft == 0:
                ax2.scatter(cor_xfluc, cor_yfluc, s = 2, alpha = 0.3, color = 'wheat', label = 'Gaussian fluctuation')
            else: ax2.scatter(cor_xfluc, cor_yfluc, s = 2, alpha = 0.3, color = 'wheat')

        fig1.suptitle(f'Digitized Data {plt_title[i]} ')
        ax1.set_title(f'primary: {primary}, {energy}, {sin2theta}, run {run} ')
        ax1.scatter(df_x, df_y, s  = 20, color = 'mediumblue', label = 'digitized data')
        ax1.plot(interpolate_x, interpolate_y, color = 'black', label = '1-D interpolation')
        ax1.legend()
        ax1.set_yscale('log')
        ax1.set_xlabel(rf'$\mathrm{{log_{{10}}(E_{{kin}}/GeV)}}$')
        ax1.set_ylabel(rf'$\mathrm{{E_{{deposited}}/MeV}}$')
        ax1.set_ylim(1e-4, 1e2)
        ax1.set_xlim(-3, 1)
        ax1.grid('-')

        fig2.suptitle(f'CORSIKA Data {plt_title[i]} ')
        ax2.set_title(f'primary: {primary}, {energy}, {sin2theta}, run {run} ')
        ax2.scatter(corsika_x, corsika_y, s = 1, color = 'red', label = 'CORSIKA with 1-D interp')
        ax2.legend()
        ax2.set_yscale('log')
        ax2.set_xlabel(rf'$\mathrm{{log_{{10}}(E_{{kin}}/GeV)}}$')
        ax2.set_ylabel(rf'$\mathrm{{E_{{deposited}}/MeV}}$')
        ax2.set_ylim(1e-4, 1e2)
        ax2.set_xlim(-3, None)
        ax2.grid('-')

        fig1.savefig(f'ScintillatorResponse/{savefile[i]}_digitized', bbox_inches='tight', dpi=400)
        fig2.savefig(f'ScintillatorResponse/{savefile[i]}_CORSIKA', bbox_inches='tight', dpi=400)