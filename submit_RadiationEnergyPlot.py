#!/usr/bin/env python3

from matplotlib.backends.backend_pdf import PdfPages  
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
from matplotlib.colors import LogNorm


plt.rcParams["font.size"] = 14
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
})

def pltNeXmaxRadEAllzenith( p, e, filtering, style, corr, pdf_pages=None):
    sin2thetas = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7",  "0.8", "0.9"]
    X0 = 697.6 #g/cm^2

    if style == 'verti': 
        fig, axs = plt.subplots(nrows=5, ncols=2, figsize=(10, 16), sharex = True, sharey = True)
        fig.subplots_adjust(wspace=0.15, hspace = 0.15, top = 0.95)
    elif style == 'horiz': 
        fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 12), sharex = False, sharey = True)
        # fig, axs = plt.subplots(nrows=2, ncols=5, figsize=(20, 12), sharex = False, sharey = False)
        fig.subplots_adjust(wspace=0.05, hspace = 0.1, top = 0.93)


    fig.suptitle(f'{p}, {e}, filtering = {filtering}, R > 0', fontsize=16)

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
                
        fXmax = f'Xmax/{p}_{e}_sin2_{sin2theta}.dat'
        fRadE = f'radEnergy/norm_sintheta2/{p}_{e}_{sin2theta}.npz'
        fRadE_unnorm = f'radEnergy/raw/{p}_{e}_{sin2theta}.npz'
        fMuon = f'GroundTotalParticles/{p}_{e}_{sin2theta}.npz'
        fEdep = f'TotalEdepScint/0m/{p}_{e}_{sin2theta}.npz'

        fileXmax = np.loadtxt(fXmax)
        Ne = fileXmax[:,5] # Ne at Xmax
        Xmax = fileXmax[:,6]
        runnum = fileXmax[:,0]
        fileRadE = np.load(fRadE, allow_pickle=True)
        fileRadE_unnorm = np.load(fRadE_unnorm, allow_pickle=True)
        fileMuon = np.load(fMuon, allow_pickle= True)
        fileEdep = np.load(fEdep, allow_pickle=True)
        Nmu_ground = fileMuon['nMu'] # number of total +-mu at ground
        Edep = fileEdep["Edep_tot"]


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
        

        alpha = fileRadE['alpha']
        sin2alpha = np.sin(alpha)**2
        thetafloat = float(sin2theta)
        theta = np.arcsin(np.sqrt(thetafloat))
        costheta = np.cos(theta)

        Xv = X0/ costheta

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
        Edep = Edep[mask]
        Nmu_ground = Nmu_ground[mask]

        ax = axs[r, c]

        if corr == 'NeErad_norm':
            sc = ax.scatter(Ne, RadE, s=7, c=Xmax, cmap="viridis", vmin = None, vmax = None)
            ax.set_yscale('log')
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
            ax.set_yscale('log')
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
        
        elif corr == 'EradEdep':
            sc = ax.scatter(Edep, RadE, s=7, c= Nmu_ground, cmap="viridis", vmin = None, vmax = None)
            ax.set_yscale('log')

            if style == 'verti':
                if c == 0: ax.set_ylabel(rf" $\mathrm{{E_{{rad}}}}$ (eV) / $\sin^2\alpha$")
                if r == 4: ax.set_xlabel(rf"$\mathrm{{E_{{deposited}}}}$ (eV)")
                plt.colorbar(sc, label=rf"$N_{{\mu,ground}}$")

            elif style == 'horiz': 
                if c == 0: ax.set_ylabel(rf"$\mathrm{{E_{{rad}}}}$ (eV)/ $\sin^2\alpha$")
                ax.set_xlabel(rf"$\mathrm{{E_{{deposited}}}}$ (eV)")
                plt.colorbar(sc, label=rf"$N_{{\mu,ground}}$", orientation='horizontal')

        

        ax.set_title(rf"$sin^2 \theta = ${sin2theta} ({np.rad2deg(theta):.2f}$^\circ$)")
                
         
    # save_dir = rf'figure/RadiationEnergy/{p}/{e}/{p}_{e}_{corr}'
    save_filename = f'{p}_{e}_{corr}.png'
    if pdf_pages:
        pdf_pages.savefig(fig, bbox_inches="tight")

    png_path = os.path.join(output_dir, save_filename)
    plt.savefig(png_path, bbox_inches="tight", dpi = 300)
    plt.close(fig)

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--particle", type=str, help="Primary particle type")
parser.add_argument("-e", "--energy", type=str, help="Energy of primary particle ex. lgE_16.0")
args = parser.parse_args()
p = args.particle
e = args.energy

style = 'horiz'
corr =['NeErad', 'NeErad_norm', 'NeXmax', 'EradSina', 'EvxBsina', 'EvxBsina_norm', 'EvxvxBsina', 'EvxvxBNe', 'EradEdep']
# corr = ['EradEdep']
minR = 0

output_dir = f'/data/user/wkammeem/CORSIKA/figures/RadiationEnergy/{p}/{e}/'
if not os.path.exists(output_dir):
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)

output_path = os.path.join(output_dir, f'{p}_{e}_AllPlots.pdf')

# pltNeXmaxRadEAllzenith(p, e, filtering = True, style=style, corr=corr[0], pdf_pages=None)

with PdfPages(output_path) as pdf:
    for i in corr:
        print(f"Plotting: {i}")
        pltNeXmaxRadEAllzenith(p, e, filtering = True, style=style, corr=i, pdf_pages=pdf)
        