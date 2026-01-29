import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import integrate
from scipy import signal, fft, constants, optimize


#=====================================
# set primary particle parameters
#=====================================

primary = "proton"
energy = "lgE_17.0"
sin2theta = "0.4"
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
fp_Nmu = lambda primary, energy, sin2theta, runnum: f'Particles/{primary}/{energy}/sin2_{sin2theta}/DAT{runnum:06d}_mupm.dat'
fp_Ne = lambda primary, energy, sin2theta, runnum: f'Particles/{primary}/{energy}/sin2_{sin2theta}/DAT{runnum:06d}_epm.dat'
fp_Nmu_tot = lambda primary, energy, sin2theta: f'Particles/{primary}/{energy}/sin2_{sin2theta}/TOTAL_mupm.dat'
fp_Ne_tot = lambda primary, energy, sin2theta: f'Particles/{primary}/{energy}/sin2_{sin2theta}/TOTAL_epm.dat'
fp_RadE = lambda primary, energy, sin2theta: f'radEnergy/{primary}_{energy}_{sin2theta}.dat'
fp_radE_op = lambda primary, energy, sin2theta: f'radEnergy/{primary}_{energy}_{sin2theta}.dat'


#=====================================
# constants / math functions
#=====================================

# constants
e0 = 8.854e-12 # vacuum permittivity constant (F/m)
c = 2.99e8 # speed of light in vacuum (m/s)
bin_width = 2e-10 #Time resolution

# functions

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
    mask = (np.abs(vxB) < 1e-1) & (vxvxB >= 0) #select x and y position only phi = 90
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

    def __init__(self, E_sum, E_sum_f):
        self.E_sum = E_sum
        self.E_sum_f = E_sum_f

    def allbands(self, time0, time1):
        plt.gcf().text(
            0.54, 0.55,
            rf'all frequency bands' + '\n' 
            + rf'bin width = {time1-time0:.2e}' + '\n'
            + rf'$\Sigma |\mathrm{{E_t}}|^2$ = {self.E_sum:.2e} '
              r'$\mathrm{V^2\,m^{-2}}$' + '\n'
            + f'energy fluence = {energyfluence(self.E_sum):.2e} '
              r'$\mathrm{eV\,m^{-2}}$',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=11
        )

    def filtered(self):
        plt.gcf().text(
            0.54, 0.38,
            rf'70-350 MHz' + '\n' 
            + rf'$\Sigma |\mathrm{{E_t}}|^2$ = {self.E_sum_f:.2e} '
              r'$\mathrm{V^2\,m^{-2}}$' + '\n'
            + f'energy fluence = {energyfluence(self.E_sum_f):.2e} '
              r'$\mathrm{eV\,m^{-2}}$',
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
    # ax1.set_xlim(-1000, 1000)
    # ax1.set_ylim(-1000, 1000)
    
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
def pltef(vxB, vxvxB, ef, eff, method):
    vxB = np.asarray(vxB)
    vxvxB = np.asarray(vxvxB)
    ef = np.asarray(ef)
    eff = np.asarray(eff)

    
    # select data on phi = 90 
    mask = (np.abs(vxB) < 1e-1) & (vxvxB >= 0) 
    # mask = (np.round(vxB, 3) == 0) & (vxvxB >= 0) #select x and y position only phi = 90
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

    plt.scatter(r,f1, label = 'all frequency bands', color = 'blue')
    plt.scatter(r,f2, label = '70-350 Mz', color = 'red')
    plt.xlabel("distance to shower axis (m)")
    plt.ylabel(rf"energy fluence  ($\mathrm{{eV}} \cdot \mathrm{{m}}^{{-2}}$)")
    if method == 'bp':
        plt.title(rf"Bandpass Filter")
    if method == 'fft':
        plt.title(rf"FFT")
    plt.legend()
    plt.text(
        0.45, 0.75,
        rf'$E_{{\rm rad_{{all \ bands}} }}= {rad_E1:.2e}$ eV' + '\n'
        rf'$E_{{\rm rad_{{70-350 MHz}} }} = {rad_E2:.2e}$ eV',
        transform=plt.gca().transAxes,
        va='top'
    )
    plt.show()
    return len(f1), len(f2)


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
    axes[1].set_xlabel(rf"$\hat{{v}}\times\hat{{B}} \ direction (m)$")
    axes[1].set_ylabel(rf"$\hat{{v}} \times (\hat{{v}}\times\hat{{B}}) \ direction \ (m)$")
    axes[1].set_title("Shower Coordinates")
    
    cbar_ax = inset_axes(axes[1],
                         width="5%", # width = 5% of parent axes width
                         height="100%", # height = 100% of parent axes height
                         loc='right', # fixed location
                         borderpad=-3 # padding between the axes and colorbar
                        )
    cbar = fig.colorbar(plot, cax=cbar_ax)
    cbar.set_label(rf'energy fluence $(\mathrm{{eV}} \cdot \mathrm{{m}}^{{-2}})$', rotation=-90, labelpad=30)
    plt.show()


# plot |E|^2 histogram
def pltEmag2(timer, timef, Emag2, Emag2_f, E_sum, E_sum_f, style):

    plt.figure(figsize=(8,5))
    
    if style == 'plot':
        plt.plot(timer, Emag2, label='all frequency bands', color = 'blue')
        plt.plot(timef, Emag2_f, label='70–350 MHz', color = 'red')
    elif style == 'scatter':
        plt.scatter(timer, Emag2, label='all frequency bands', color = 'blue', s = 5)
        plt.scatter(timef, Emag2_f, label='70–350 MHz', color = 'red', s = 5)
    
    plt.xlabel("Time (s)")
    plt.ylabel(r"$|\mathrm{E_t}|^2\ (\mathrm{V^2\,m^{-2}})$")
    plt.legend()

    tb = textboxes(E_sum, E_sum_f)
    tb.allbands(timer[0], timer[1])
    tb.filtered()

    plt.show()
    
# plot longitudinal profile
def pltlp(atmdepth, positron, electron, muplus, muminus, tot_e, tot_mu):
    color_mu =  'steelblue'
    color_e = 'firebrick'
    
    fig, axes = plt.subplots(1, 3, sharey=True)
    fig.subplots_adjust(wspace=0.1)  
    
    axes[0].plot(positron, atmdepth , label = r'$e^+$', color = color_e, ls ='--')
    axes[0].plot(electron, atmdepth,  label = r'$e^-$', color = color_e, ls =':')
    axes[0].legend()
    axes[0].set_ylabel(r"atmospheric depth (g/$\mathrm{cm}^2$)")
    
    axes[1].plot(muplus, atmdepth, label = r'$\mu^+$', color = color_mu , ls ='--')
    axes[1].plot(muminus, atmdepth, label = r'$\mu^-$',  color = color_mu, ls = ':')
    axes[1].legend()
    axes[1].set_xlabel("particle number")
    axes[1].set_title(f'{primary}, {energy}, sin2theta = {sin2theta}')
    
    axes[2].plot(tot_mu, atmdepth, label = r'$\mu^{\pm} (\times 50)$', color = color_mu)
    axes[2].plot(tot_e, atmdepth, label = r'$e^{\pm}$', color = color_e)
    axes[2].legend()

    plt.gca().invert_yaxis()
    plt.show()

# plot correlation of Nmu and Ne
def pltmuecorr(energydir, sin2theta, energylabel):
    protoncolor = 'red'
    ironcolor = 'blue'
    scattersize = 5
    
    protonpath = f'Particles/proton/{energydir}/sin2_{sin2theta}'
    mudata = np.loadtxt(protonpath + '/TOTAL_mupm.dat')
    ptotmu = mudata[:,5]
    edata = np.loadtxt(protonpath + '/TOTAL_epm.dat')
    ptote = edata[:,5]
    
    plt.scatter(ptote, ptotmu, color = protoncolor, s = scattersize, label = 'proton')
    
    ironpath = f'Particles/iron/{energydir}/sin2_{sin2theta}'
    mudata = np.loadtxt(ironpath + '/TOTAL_mupm.dat')
    fetotmu = mudata[:,5]
    edata = np.loadtxt(protonpath + '/TOTAL_epm.dat')
    fetote = edata[:,5]
    
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