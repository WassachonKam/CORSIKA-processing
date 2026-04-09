import argparse
import numpy as np
import matplotlib.pyplot as plt
from corsikaio import CorsikaParticleFile
import pandas as pd
import random

############################# FUNCTIONS AND CONSTANTS #############################

# kinetic energy 
def Ekin(px,py,pz,m):
    return np.sqrt(px**2 + py**2 + pz**2 + m**2) -m 

# constants
mumass = 0.105658 #muon mass in GeV/c^2
emass = 5.11e-4 #electron mass in GeV/c^2

############################# DATA PREPARATION #############################

primary = 'proton'
energies = [f"lgE_{x/10:.1f}" for x in range(160, 181)] # energy bin lgE_16.0, lgE_16.1, ..., lgE_18.0
sin2thetas = [f"sin2_{x/10:.1f}" for x in range(0, 10)] # zenith bin sin2_0.0, sin2_0.1, ..., sin2_0.9
# sin2thetas = ["sin2_0.0"]
# energies = ['lgE_16.0']
runnums = 0

for t in range(len(sin2thetas)):
    print(f'processing {sin2thetas[t]}')

    # dictionary for storing data
    all_data = {
                "energy": np.arange(16.0, 18.05, 0.1), # energy array 16.0, 16.1, ..., 18.0 
                "Edep_e": [],
                "Edep_mu": [],
                "Edep_epm": []

            }

    for e in range (len(energies)):

        file_path = (f'/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energies[e]}/{sin2thetas[t]}/{runnums:06d}/DAT{runnums:06d}')
        file_input = (f"/data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern/{primary}/{energies[e]}/{sin2thetas[t]}/{runnums:06d}/SIM{runnums:06d}.inp")
        
        try:
            with open(file_input) as f:
                for line in f:
                    parts = line.split()
                    if parts[0] == "THETAP":
                        thetap = float(parts[1])
        except FileNotFoundError:
            all_data["energy"] = np.delete(all_data["energy"], e)
            continue 

        try:
            with CorsikaParticleFile(file_path, thinning= True) as file:
                # we only have one event per file, we can grab it like this
                event = next(file)

        except (OSError, IOError) as err:
            print(f"Skipping {energies[e]}: {err}")
            all_data["energy"] = np.delete(all_data["energy"], e)
            continue
        except StopIteration:
            print(f"Skipping {energies[e]}: File is empty.")
            all_data["energy"] = np.delete(all_data["energy"], e)
            continue

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
        dfe = pd.read_csv('electron_0.0.csv')
        dmu = pd.read_csv('muon_0.0.csv')
        digitized_files = [dfe, dmu]
        Ek_array = [Ek_electron, Ek_mu]
        plt_title = [rf'$e^{{-}}$', rf'$\mu^{{\pm}}$']
        savefile = ['electron', 'muon']

        for i in range(len(digitized_files)):

            df = digitized_files[i]
            df_x = df['x']
            df_y = df[' y'] / np.cos(theta)

            # interpolate digitized data
            interpolate_x = np.linspace(min(df_x), max(df_x), num = 400)
            interpolate_y = np.interp(interpolate_x, df_x, df_y)

            # Deposited energy as a function of Ek from CORSIKA file
            logEk = np.log10(Ek_array[i])
            corsika_x = logEk
            corsika_y = np.interp(corsika_x, df_x, df_y)

            total_Edep = sum(corsika_y)

            if i == 0: all_data["Edep_e"].append(total_Edep)
            elif i == 1: all_data["Edep_mu"].append(total_Edep)

    output_path = f"TotalEdep/{sin2thetas[t]}"
    np.savez_compressed(output_path, **all_data)


            # # make a scatter plot
            # fig1, ax1 = plt.subplots() 
            # fig2, ax2 = plt.subplots()

            # # generate randon Gaussian fluctuation
            # Edep = interpolate_y
            # for ft in range (len(Edep)):
            #     yfluc = []
            #     xfluc = []
            #     for j in range(20):
            #         fluct = random.gauss(1, 0.5)
            #         fEdep = Edep[ft]*fluct
            #         yfluc.append(fEdep)
            #         xfluc.append(interpolate_x[ft])
            #     if ft == 0:
            #         ax1.scatter(xfluc, yfluc, s = 2, alpha = 0.3, color = 'lightskyblue', label = 'Gaussian fluctuation')
            #     else: ax1.scatter(xfluc, yfluc, s = 2, alpha = 0.3, color = 'lightskyblue')
            

            # interpolate_xcorr = np.linspace(min(corsika_x), max(corsika_x), num = 400)
            # interpolate_ycorr = np.interp(interpolate_xcorr, df_x, df_y)
            # Edep = interpolate_ycorr
            
            # for ft in range (len(Edep)):
            #     cor_yfluc = []
            #     cor_xfluc = []
            #     for j in range(20):
            #         fluct = random.gauss(1, 0.5)
            #         fEdep = Edep[ft]*fluct
            #         cor_yfluc.append(fEdep)
            #         cor_xfluc.append(interpolate_xcorr[ft])
            #     if ft == 0:
            #         ax2.scatter(cor_xfluc, cor_yfluc, s = 2, alpha = 0.3, color = 'wheat', label = 'Gaussian fluctuation')
            #     else: ax2.scatter(cor_xfluc, cor_yfluc, s = 2, alpha = 0.3, color = 'wheat')

            # fig1.suptitle(f'Digitized Data {plt_title[i]} ')
            # ax1.set_title(f'primary: {primary}, {energies[e]}, {sin2thetas[t]}, run {runnums} ')
            # ax1.scatter(df_x, df_y, s  = 20, color = 'mediumblue', label = 'digitized data')
            # ax1.plot(interpolate_x, interpolate_y, color = 'black', label = '1-D interpolation')
            # ax1.legend()
            # ax1.set_yscale('log')
            # ax1.set_xlabel(rf'$\mathrm{{log_{{10}}(E_{{kin}}/GeV)}}$')
            # ax1.set_ylabel(rf'$\mathrm{{E_{{deposited}}/MeV}}$')
            # ax1.set_ylim(1e-4, 1e2)
            # ax1.set_xlim(-3, 1)
            # ax1.grid('-')

            # fig2.suptitle(f'CORSIKA Data {plt_title[i]} ')
            # ax2.set_title(f'primary: {primary}, {energies[e]}, {sin2thetas[t]}, run {runnums} ')
            # ax2.scatter(corsika_x, corsika_y, s = 1, color = 'red', label = 'CORSIKA with 1-D interp')
            # ax2.legend()
            # ax2.set_yscale('log')
            # ax2.set_xlabel(rf'$\mathrm{{log_{{10}}(E_{{kin}}/GeV)}}$')
            # ax2.set_ylabel(rf'$\mathrm{{E_{{deposited}}/MeV}}$')
            # ax2.set_ylim(1e-4, 1e2)
            # ax2.set_xlim(-3, None)
            # ax2.grid('-')

            # fig1.savefig(f'{savefile[i]}_digitized', bbox_inches='tight', dpi=400)
            # fig2.savefig(f'{savefile[i]}_CORSIKA', bbox_inches='tight', dpi=400)
