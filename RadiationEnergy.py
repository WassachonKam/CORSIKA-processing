from scipy import signal, fft, constants, optimize
from functions import *

output_path = fp_radE_op(primary, energy, sin2theta)


# antenna positions
ant_x = lambda xpos: xpos.astype(float)/1e2 # x position in m
ant_y =  lambda ypos: ypos.astype(float)/1e2 # y position in m

runmin, runmax = 0, 199

open(output_path, "w").close() #clean file
with open(output_path, 'a') as file:
    file.write('#run radE(eV) radE_filtered(eV)' + '\n') #header

for run in range(runmin, runmax+1):

    Nant = 160 #number of total antennas
    # create empty arrays of energy fluence, filtered energy fluence, vxB, and vxvxB
    ef, eff,  = (np.zeros(Nant) for _ in range(2))
    # 160 antenna (Nant = 160) loop
    for i, ant_no in enumerate(range (1,Nant+1)):
        
        radio_dat = np.loadtxt(fp_radio(primary, energy, sin2theta, run, ant_no))
    
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
    flist = fp_list(primary, energy, sin2theta, run)
    xypos = pd.read_csv(flist, sep=r"\s+",header=None)
    antx = ant_x(xypos[2])
    anty = ant_y(xypos[3])

    finp = fp_inp(primary, energy, sin2theta, run)
    
    # append output on file
    with open(output_path, 'a') as file:
        file.write( f'{run}' + ' '
                   + f'{RadEnergy(finp, Nant, ef, eff, antx, anty)[0]}' 
                   + ' '
                   + f'{RadEnergy(finp, Nant, ef, eff, antx, anty)[1]}' 
                   + '\n')
    print(run, len(RadEnergy(finp, Nant, ef, eff, antx, anty)[2]), len(RadEnergy(finp, Nant, ef, eff, antx, anty)[3]))