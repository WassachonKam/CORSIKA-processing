from scipy import signal, fft, constants, optimize
from functions import *


output_path = fp_RadE(primary, energy, sin2theta)


# antenna positions
ant_x = lambda xpos: xpos.astype(float)/1e2 # x position in m
ant_y =  lambda ypos: ypos.astype(float)/1e2 # y position in m

runmin, runmax = 0, 199

# create dictionaries to store all data
all_data = {
    "run": [],
    "radE(eV)": [],
    "radE_filtered(eV)": [],
    "antx": [],
    "anty": [],
    "ef": [],
    "eff": []
}

print(f"processing radiation energy {primary} {energy} {sin2theta} ")
for run in range(runmin, runmax+1):

    Nant = 160 #number of total antennas
    # create empty arrays of energy fluence, filtered energy fluence, vxB, and vxvxB
    ef, eff,  = (np.zeros(Nant) for _ in range(2))
    # 160 antenna (Nant = 160) loop
    for i, ant_no in enumerate(range (1,Nant+1)):

        try: radio_dat = np.loadtxt(fp_radio(primary, energy, sin2theta, run, ant_no))
        except FileNotFoundError: continue
            
    
        # extract data
        time = radio_dat[:,0]
        E_SI = 29979 # electric field converter from cgs unit to SI unit 1 statV/cm≈29,979 V/m
        Ex, Ey, Ez = radio_dat[:,1]*E_SI, radio_dat[:,2]*E_SI, radio_dat[:,3]*E_SI # north, west, vertical electric field in V/m
    
        #filter electric field using signal.butter
        b, a = signal.butter(5, [70e6, 350e6] , 'bp', fs=5e9)
        Ex_f, Ey_f, Ez_f = (signal.filtfilt(b, a, i) for i in (Ex, Ey, Ez))
        Emag2, Emag2_f = magE2(Ex, Ey, Ez), magE2(Ex_f, Ey_f, Ez_f)  # |E|^2 raw and filtered
        E_sum, E_sum_f = sum(Emag2), sum(Emag2_f) # sum(|E|^2)
         
        # energy fluence 
        ef[i], eff[i] = energyfluence(E_sum), energyfluence(E_sum_f)
    
    # x-y antenna position
    try:
        flist = fp_list(primary, energy, sin2theta, run)
        xypos = pd.read_csv(flist, sep=r"\s+", header=None)
    except FileNotFoundError:
        continue   # skip this run entirely
    
    antx = ant_x(xypos[2])
    anty = ant_y(xypos[3])

    try:
        finp = fp_inp(primary, energy, sin2theta, run)
    except FileNotFoundError:
        continue
    
    radE_vals = RadEnergy(finp, Nant, ef, eff, antx, anty)
    
    # store data in dictionary
    all_data["run"].append(run)
    all_data["radE(eV)"].append(radE_vals[0])
    all_data["radE_filtered(eV)"].append(radE_vals[1])
    all_data["ef"].append(ef)
    all_data["eff"].append(eff)
    all_data["antx"].append(antx)
    all_data["anty"].append(anty)

    # print(run, len(radE_vals[2]), len(radE_vals[3]))
    print(f"run {run}")

# save everything in a compressed npz file
np.savez_compressed(output_path, **all_data)
print(f"Saved all data to {output_path}")