#%% downnload modules
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
from matplotlib.backends.backend_pdf import PdfPages  
from sklearn.metrics import mean_squared_error, r2_score
from scipy.optimize import curve_fit
from scipy.interpolate import griddata

#%% parameter setting
#=====================================
# set primary particle parameters
#=====================================

primary = "proton"
energy = "lgE_17.0"
sin2theta = "0.4"
runnum = 0

#%% file path and constants
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
fp_RadE = lambda primary, energy, sin2theta: f'radEnergy/raw/{primary}_{energy}_{sin2theta}.npz'
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
X0 = 697.6

#%% radiation energy

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

#%% regrssion function
def linearEdep(X, alpha, beta):
	ne, nmu = X
	Edep = alpha *ne + beta *nmu
	return Edep

def linearErad_NeXmax(NeXmax, gamma, b):
	Erad = gamma*NeXmax + b
	return Erad

def ExpoNeRatio(dXmax, delta,A):
	return A*np.exp(-dXmax/ delta)

def Pred_Nmu(alpha, beta, delta, gamma, A, b, Edep, Erad, dXmax):
	return (1/beta)*(Edep - ((alpha*A/gamma)*np.exp(-dXmax/delta)*(Erad-b)))

def RegressionFixedEnergyBin(ThreshR, filteringXmax, energy_all, lgE_GeV_bin, NmuNorm):
	df = pd.read_parquet(f'RandomForestRegression/data_for_regression_{ThreshR}_Xmax_{filteringXmax}_NmuNorm_{NmuNorm}.parquet')
	titlelabel =  f'removing underground Xmax = {filteringXmax}, R > {ThreshR}'

	# binning energy and zenith angle 
	costheta = df['costheta']
	energy = df['energy']
	particle = df['particle']
	sin2theta = np.sin(np.arccos(costheta))**2
	zenith_bins = np.linspace(0.0, 0.9, 10)
	energy_bins = np.linspace(7.0, 9.0, 21)
	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)



	particle_bins = ['proton', 'helium', 'oxygen', 'iron']

	colors = ['red','gold', 'green', 'blue']
	lgE_bin = lgE_GeV_bin - 9

	cols, rows = 1, 7
	fig_ax, axes = plt.subplots(rows, cols, figsize=(7, 2 * rows), constrained_layout=True, sharex= True)
	axes = axes.flatten() # one sequence list
	fig, ax = plt.subplots(figsize=(8,5), constrained_layout=True)

	# remove showers in sin2_0.9 bin
	zenith_bins = zenith_bins[:-1]
	handle_res = []
	handle_bias = []
	for i in range (len(particle_bins)):

		particle_bin = particle_bins[i]
		# fit parameters
		params = {'alpha': [], 'alpha_sd': [],
				'beta': [], 'beta_sd': [],
				'delta': [], 'delta_sd': [],
				'gamma': [], 'gamma_sd': [],
				'A': [], 'A_sd': [],
				'b': [], 'b_sd': [],
				'R2': []}

		bias_array = []
		reso_array = []


		for sin2_bin in zenith_bins:
			# masking
			sin2_index = np.digitize([sin2_bin], zenith_bins)[0]
			lgE_index = np.digitize([lgE_bin], energy_bins)[0]
			mask_particle = (particle == particle_bin)
			mask_energy = (energy_indices == lgE_index)
			mask_zenith = (zenith_indices == sin2_index)
			energy_label = rf'$\mathrm{{lgE = {lgE_GeV_bin}}}$'

			combined_mask = mask_particle  & mask_zenith & mask_energy

			# df_mask = df
			if energy_all == True:
				combined_mask = mask_particle & mask_zenith
				energy_label = rf'$\mathrm{{lgE = all \ bins}}$'
			df_mask = df[combined_mask]

			# define parameters
			Edep = df_mask['Edep']
			Erad = df_mask['Erad']
			Xmax = df_mask['Xmax']
			Ne_ground = df_mask['Ne_ground']
			Ne_Xmax = df_mask['Ne']
			Nmu_ground = df_mask['Nmu']
			costheta = df_mask['costheta']
			Ne_ground_Xmax = Ne_ground/Ne_Xmax
			Xv = X0/ costheta
			dXmax = Xv - Xmax

			# fitting 
			X_Nmu_Ne = (Ne_ground, Nmu_ground)
			popt, pcov = curve_fit(linearEdep, X_Nmu_Ne, Edep)
			perr = np.sqrt(np.diag(pcov))
			alpha, beta = popt
			alpha_sd, beta_sd = perr[0], perr[1] 

			popt2, pcov2 = curve_fit(linearErad_NeXmax, Ne_Xmax, Erad)
			perr2 = np.sqrt(np.diag(pcov2))
			gamma , b, gamma_sd , b_sd = popt2[0], popt2[1], perr2[0], perr2[1]

			popt3, pcov3 = curve_fit(ExpoNeRatio, dXmax, Ne_ground_Xmax, p0 = [250,1])
			perr3 = np.sqrt(np.diag(pcov3))
			delta, A, delta_sd, A_sd = popt3[0], popt3[1], perr3[0], perr3[1]

			Nmu_pred  = Pred_Nmu(alpha, beta, delta, gamma, A, b, Edep, Erad, dXmax)

			mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
			Nmu_pred = Nmu_pred[mask_Nmu_pred]
			Nmu_ground = Nmu_ground[mask_Nmu_pred]

			#print(Nmu_ground, Nmu_pred)
			R2 = r2_score(Nmu_ground, Nmu_pred)
			
			residual = Nmu_pred - Nmu_ground
			bias = np.mean(residual)
			reso = np.std(residual)

			bias_array.append(bias)
			reso_array.append(reso)
			

			params['alpha'].append(alpha)
			params['beta'].append(beta)
			params['delta'].append(delta)
			params['gamma'].append(gamma)
			params['A'].append(A)
			params['b'].append(b)
			params['alpha_sd'].append(alpha_sd)
			params['beta_sd'].append(beta_sd)
			params['delta_sd'].append(delta_sd)
			params['gamma_sd'].append(gamma_sd)
			params['A_sd'].append(A_sd)
			params['b_sd'].append(b_sd)
			params['R2'].append(R2)
			
			# plotting
			color = 'steelblue'
			# if sin2_bin == 0.0:
				# fig, ax = plt.subplots(figsize= (6,5), dpi = 150)
				# plt.scatter(Nmu_ground, Nmu_pred, s = 5, label = rf'$N_{{\mu}}$', color = color)
				# plt.xlabel(rf' $N_{{\mu}}^{{true}}$')
				# plt.ylabel(rf' $N_{{\mu}}^{{predicted}}$')
				# plt.text(0.7, 0.15, rf'$R^2$ = {R2:.3f}', transform=ax.transAxes)
				# x = np.linspace(min(Nmu_ground), max(Nmu_ground))
				# plt.plot(x,x, color = 'red', label = rf'$N_{{\mu}}^{{predicted}}$ = $N_{{\mu}}^{{true}}$')
				# plt.legend()
				# fig_ax.suptitle(rf'{particle_bin} lgE_{lgE_GeV_bin}')

		size = 5

		bias = ax.scatter(zenith_bins, bias_array, color = colors[i], marker = 'x', label = particle_bin )
		res = ax.scatter(zenith_bins, reso_array, color = colors[i], marker = '^', label = particle_bin)

		handle_res.append(res)
		handle_bias.append(bias)

		axes[0].errorbar(zenith_bins, params['alpha'], params['alpha_sd'],fmt='o', markersize = size, color = colors[i], label = particle_bin )
		axes[1].errorbar(zenith_bins, params['beta'], params['beta_sd'], fmt='o', markersize = size, color = colors[i])
		axes[2].errorbar(zenith_bins, params['delta'], params['delta_sd'], fmt='o', markersize = size, color = colors[i])
		axes[3].errorbar(zenith_bins, params['gamma'], params['gamma_sd'], fmt='o', markersize = size, color = colors[i])
		axes[4].errorbar(zenith_bins, params['A'], params['A_sd'], fmt='o', markersize = size, color = colors[i])
		axes[5].errorbar(zenith_bins, params['b'], params['b_sd'], fmt='o', markersize = size, color = colors[i])
		axes[6].scatter(zenith_bins, params['R2'], s = 15, color = colors[i])

		axes[0].set_ylabel(rf'$\alpha$')
		axes[1].set_ylabel(rf'$\beta$')
		axes[2].set_ylabel(rf'$\delta \ (\mathrm{{g/cm^2}})$')
		axes[3].set_ylabel(rf'$\gamma$')
		axes[4].set_ylabel(rf'A')
		axes[5].set_ylabel(rf'b')
		axes[6].set_ylabel(rf'$R^2$')
		axes[6].set_xlabel(rf'$\mathrm{{sin^2 \theta}}$')

	fig_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig_ax.suptitle(rf'lgE_{lgE_GeV_bin} eV (filtered Xmax)')

	fig_ax.suptitle(titlelabel + '\n' + energy_label)
	fig.suptitle(titlelabel + '\n' + energy_label)

	leg1 = ax.legend(bbox_to_anchor=(1.01, 1), handles= handle_res, loc='upper left', title="Resolution")
	ax.add_artist(leg1) 
	ax.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_bias, loc='upper left', title="Bias")
	ax.set_xlabel(rf'$\mathrm{{sin^2 \theta}}$')
	ax.set_ylabel(rf'$\mathrm{{log(N_{{\mu^{{\pm}}}}^{{true}}) - log(N_{{\mu^{{\pm}}}}^{{predicted}})}}$')
	ax.set_title(f'')
	ax.set_ylim(-0.16,0.16)
	ax.hlines(y= 0, xmin = min(zenith_bins), xmax= max(zenith_bins), colors = 'k', linestyles =  '--')

def RegressionFixedZenithBin(ThreshR, filteringXmax, zenith_all, sin2_bin, NmuNorm):
	df = pd.read_parquet(f'RandomForestRegression/data_for_regression_{ThreshR}_Xmax_{filteringXmax}_NmuNorm_{NmuNorm}.parquet')
	titlelabel =  f'removing underground Xmax = {filteringXmax}, R > {ThreshR}'

	# binning energy and zenith angle 
	costheta = df['costheta']
	energy = df['energy']
	particle = df['particle']
	sin2theta = np.sin(np.arccos(costheta))**2
	zenith_bins = np.linspace(0.0, 0.9, 10)
	energy_bins = np.linspace(7.0, 9.0, 21)

	# remove showers in sin2_0.9 bin
	zenith_bins = zenith_bins[:-1]

	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)

	particle_bins = ['proton', 'helium', 'oxygen', 'iron']
	colors = ['red','gold', 'green', 'blue']
	
	sin2_index = np.digitize([sin2_bin], zenith_bins)[0]

	cols, rows = 1, 7
	fig_ax, axes = plt.subplots(rows, cols, figsize=(7, 2 * rows), constrained_layout=True, sharex= True)
	axes = axes.flatten()
	fig, ax = plt.subplots(figsize=(8,5), constrained_layout=True)
	fig2, ax2 = plt.subplots(figsize=(8,5), constrained_layout=True)

	handle_bias = []
	handle_res = []
	handle_true = []
	handle_recon = []

	for i, particle_bin in enumerate(particle_bins):

		energy_array = []
		bias_array =[]
		reso_array = []
		Nmu_pred_array = []
		Nmu_ground_array =[]


		params = {'alpha': [], 'alpha_sd': [],
				'beta': [], 'beta_sd': [],
				'delta': [], 'delta_sd': [],
				'gamma': [], 'gamma_sd': [],
				'A': [], 'A_sd': [],
				'b': [], 'b_sd': [],
				'R2': []}
		
		# for sin2_bin in zenith_bins:
		for lgE_idx in range(len(energy_bins)):

			# masking
			mask_particle = (particle == particle_bin)
			mask_energy = (energy_indices == lgE_idx)
			mask_zenith = (zenith_indices == sin2_index)
			combined_mask = mask_particle & mask_energy & mask_zenith
			zenith_label = rf'$\mathrm{{sin^2 \theta = {sin2_bin}}}$'

			if zenith_all == True:
				combined_mask = mask_particle & mask_energy
				zenith_label = rf'$\mathrm{{sin^2 \theta = all \ bins}}$'
			df_mask = df[combined_mask]

			if df_mask.empty: 
				continue

			# define parameters
			Edep = df_mask['Edep']
			Erad = df_mask['Erad']
			Xmax = df_mask['Xmax']
			Ne_ground = df_mask['Ne_ground']
			Ne_Xmax = df_mask['Ne']
			Nmu_ground = df_mask['Nmu']
			costheta = df_mask['costheta']
			Ne_ground_Xmax = Ne_ground/Ne_Xmax
			Xv = X0/ costheta
			dXmax = Xv - Xmax

			# fitting 
			X_Nmu_Ne = (Ne_ground, Nmu_ground)
			popt, pcov = curve_fit(linearEdep, X_Nmu_Ne, Edep)
			perr = np.sqrt(np.diag(pcov))
			alpha, beta = popt
			alpha_sd, beta_sd = perr[0], perr[1] 

			popt2, pcov2 = curve_fit(linearErad_NeXmax, Ne_Xmax, Erad)
			perr2 = np.sqrt(np.diag(pcov2))
			gamma , b, gamma_sd , b_sd = popt2[0], popt2[1], perr2[0], perr2[1]

			popt3, pcov3 = curve_fit(ExpoNeRatio, dXmax, Ne_ground_Xmax, p0 = [250,1])
			perr3 = np.sqrt(np.diag(pcov3))
			delta, A, delta_sd, A_sd = popt3[0], popt3[1], perr3[0], perr3[1]

			Nmu_pred  = Pred_Nmu(alpha, beta, delta, gamma, A, b, Edep, Erad, dXmax)

			mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
			Nmu_pred = Nmu_pred[mask_Nmu_pred]
			Nmu_ground = Nmu_ground[mask_Nmu_pred]

			R2 = r2_score(Nmu_ground, Nmu_pred)

			residual = Nmu_pred - Nmu_ground
			bias = np.mean(residual)
			reso = np.std(residual)

			bias_array.append(bias)
			reso_array.append(reso)
			energy_array.append(energy_bins[lgE_idx])
			Nmu_ground_array.append(np.mean(Nmu_ground))
			Nmu_pred_array.append(np.mean(Nmu_pred))

			params['alpha'].append(alpha)
			params['beta'].append(beta)
			params['delta'].append(delta)
			params['gamma'].append(gamma)
			params['A'].append(A)
			params['b'].append(b)
			params['alpha_sd'].append(alpha_sd)
			params['beta_sd'].append(beta_sd)
			params['delta_sd'].append(delta_sd)
			params['gamma_sd'].append(gamma_sd)
			params['A_sd'].append(A_sd)
			params['b_sd'].append(b_sd)
			params['R2'].append(R2)

		bias = ax.scatter(energy_array, bias_array, color = colors[i], marker = 'x', label = particle_bin )
		res = ax.scatter(energy_array, reso_array, color = colors[i], marker = '^', label = particle_bin)
		
		handle_res.append(res)
		handle_bias.append(bias)
		true = ax2.scatter(energy_array, Nmu_ground_array, color = colors[i], marker = 'x', label = particle_bin, s = 50)
		recon = ax2.scatter(energy_array, Nmu_pred_array, color = colors[i], marker = 'o', label = particle_bin, s = 20,
						edgecolors = 'black', linewidths=0.5)

		handle_true.append(true)
		handle_recon.append(recon)

		axes[0].errorbar(energy_array, params['alpha'], params['alpha_sd'],fmt='o', markersize = 5, color = colors[i], label = particle_bin )
		axes[1].errorbar(energy_array, params['beta'], params['beta_sd'], fmt='o', markersize = 5, color = colors[i])
		axes[2].errorbar(energy_array, params['delta'], params['delta_sd'], fmt='o', markersize = 5, color = colors[i])
		axes[3].errorbar(energy_array, params['gamma'], params['gamma_sd'], fmt='o', markersize = 5, color = colors[i])
		axes[4].errorbar(energy_array, params['A'], params['A_sd'], fmt='o', markersize = 5, color = colors[i])
		axes[5].errorbar(energy_array, params['b'], params['b_sd'], fmt='o', markersize = 5, color = colors[i])
		axes[6].scatter(energy_array, params['R2'], s = 15, color = colors[i])

		axes[0].set_ylabel(rf'$\alpha$')
		axes[1].set_ylabel(rf'$\beta$')
		axes[2].set_ylabel(rf'$\delta \ (\mathrm{{g/cm^2}})$')
		axes[3].set_ylabel(rf'$\gamma$')
		axes[4].set_ylabel(rf'A')
		axes[5].set_ylabel(rf'b')
		axes[6].set_ylabel(rf'$R^2$')
		axes[6].set_xlabel(rf'$\mathrm{{log(E (GeV))}}$')

	fig_ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig_ax.suptitle(titlelabel + '\n' + zenith_label)
	

	fig.suptitle(titlelabel + '\n' + zenith_label)
	leg1 = ax.legend(bbox_to_anchor=(1.01, 1), handles= handle_res, loc='upper left', title="Resolution")
	ax.add_artist(leg1) 
	ax.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_bias, loc='upper left', title="Bias")
	ax.set_xlabel(rf'$\mathrm{{log(E (GeV))}}$')
	ax.set_ylabel(rf'$\mathrm{{log(N_{{\mu^{{\pm}}}}^{{true}}) - log(N_{{\mu^{{\pm}}}}^{{predicted}})}}$')
	ax.set_title(f'')
	# ax.set_ylim(-0.16,0.16)
	ax.hlines(y= 0, xmin = min(energy_array), xmax= max(energy_array), colors = 'k', linestyles =  '--')

	fig2.suptitle(titlelabel + '\n' + zenith_label)
	leg2 = ax2.legend(bbox_to_anchor=(1.01, 1), handles= handle_true, loc='upper left', title="True Value")
	ax2.add_artist(leg2) 
	ax2.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_recon, loc='upper left', title="Reconstructed Value")
	ax2.set_xlabel(rf'$\mathrm{{log(E (GeV))}}$')
	ax2.set_ylabel(rf'$\mathrm{{log(N_{{\mu^{{\pm}}}})}}$')

# fitting and return constant array by averaging fit constants across primaries
def GetConstFromRegression(ThreshR, filteringXmax, NmuNorm):
	df = pd.read_parquet(f'RandomForestRegression/data_for_regression_{ThreshR}_Xmax_{filteringXmax}_NmuNorm_{NmuNorm}.parquet')

	# binning energy and zenith angle 
	costheta = df['costheta']
	energy = df['energy']
	particle = df['particle']
	sin2theta = np.sin(np.arccos(costheta))**2
	zenith_bins = np.linspace(0.0, 0.9, 10)
	energy_bins = np.logspace(7.0, 9.0, 21)
	particle_bins = ['proton', 'helium', 'oxygen', 'iron']

	# remove showers in sin2_0.9 bin
	zenith_bins = zenith_bins[:-1]

	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)

	params_list = []
	Nmu_list = []


	for lgE_idx in range(len(energy_bins)):
		for sin2_index in range(len(zenith_bins)):

			params = {'alpha': [], 
				'beta': [], 
				'delta': [], 
				'gamma': [], 
				'A': [], 
				'b': [], 	
				'R2': [],
				'bias': [],
				'reso': []}
			

			for i, particle_bin in enumerate(particle_bins):
				# masking

				mask_particle = (particle == particle_bin)
				mask_energy = (energy_indices == lgE_idx +1) # bin 0 is less than  lgE_16.0
				mask_zenith = (zenith_indices == sin2_index +1) # bin 0 is less than sin2_0.0 
				combined_mask = mask_particle & mask_energy & mask_zenith
				
				df_mask = df[combined_mask]
				
				if df_mask.empty: 
					continue
				
				# define parameters
				Edep = df_mask['Edep']
				Erad = df_mask['Erad']
				Xmax = df_mask['Xmax']
				Ne_ground = df_mask['Ne_ground']
				Ne_Xmax = df_mask['Ne']
				Nmu_ground = df_mask['Nmu']
				costheta = df_mask['costheta']
				Ne_ground_Xmax = Ne_ground/Ne_Xmax
				Xv = X0/ costheta
				dXmax = Xv - Xmax
		
				# fitting 
				X_Nmu_Ne = (Ne_ground, Nmu_ground)
				popt, pcov = curve_fit(linearEdep, X_Nmu_Ne, Edep)
				perr = np.sqrt(np.diag(pcov))
				alpha, beta = popt
				alpha_sd, beta_sd = perr[0], perr[1] 

				popt2, pcov2 = curve_fit(linearErad_NeXmax, Ne_Xmax, Erad)
				perr2 = np.sqrt(np.diag(pcov2))
				gamma , b, gamma_sd , b_sd = popt2[0], popt2[1], perr2[0], perr2[1]

				popt3, pcov3 = curve_fit(ExpoNeRatio, dXmax, Ne_ground_Xmax, p0 = [250,1])
				perr3 = np.sqrt(np.diag(pcov3))
				delta, A, delta_sd, A_sd = popt3[0], popt3[1], perr3[0], perr3[1]

				Nmu_pred  = Pred_Nmu(alpha, beta, delta, gamma, A, b, Edep, Erad, dXmax)

				mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
				Nmu_pred = Nmu_pred[mask_Nmu_pred]
				Nmu_ground = Nmu_ground[mask_Nmu_pred]

				R2 = r2_score(Nmu_ground, Nmu_pred)

				residual = Nmu_pred - Nmu_ground
				bias = np.mean(residual)
				reso = np.std(residual)

				params['alpha'].append(alpha)
				params['beta'].append(beta)
				params['delta'].append(delta)
				params['gamma'].append(gamma)
				params['A'].append(A)
				params['b'].append(b)
				params['R2'].append(R2)
				params['bias'].append(bias)
				params['reso'].append(reso)
			
			# average constants of four primary particles
			if len(params['alpha']) > 0: # check if arrays don't empty 

				params_all = {'sin2theta': zenith_bins[sin2_index],
					'energy': energy_bins[lgE_idx],
					'alpha': np.mean(params['alpha']), 
					'beta': np.mean(params['beta']), 
					'delta': np.mean(params['delta']), 
					'gamma': np.mean(params['gamma']), 
					'A': np.mean(params['A']), 
					'b': np.mean(params['b']), 
					'R2': np.mean(params['R2']),
					'bias':np.mean(params['bias']),
					'reso': np.mean(params['reso'])}
				params_list.append(params_all)

				for i, particle_bin in enumerate(particle_bins):
					# masking

					mask_particle = (particle == particle_bin)
					mask_energy = (energy_indices == lgE_idx +1) # bin 0 is less than  lgE_16.0
					mask_zenith = (zenith_indices == sin2_index +1) # bin 0 is less than sin2_0.0 
					combined_mask = mask_particle & mask_energy & mask_zenith
					
					df_mask = df[combined_mask]
					
					if df_mask.empty: 
						continue
					
					# define parameters
					Edep = df_mask['Edep']
					Erad = df_mask['Erad']
					Xmax = df_mask['Xmax']
					Nmu_ground = df_mask['Nmu']
					costheta = df_mask['costheta']
					Xv = X0/ costheta
					dXmax = Xv - Xmax
					Nmu_pred = Pred_Nmu(np.mean(params['alpha']), np.mean(params['beta']), np.mean(params['delta'])
						, np.mean(params['gamma']), np.mean(params['A']), np.mean(params['b']), Edep, Erad, dXmax)
					
					residual = Nmu_ground - Nmu_pred
					bias = np.mean(residual)
					reso = np.std(residual)

					Nmu_all = {'particle': particle_bin, 
			   					'sin2theta': zenith_bins[sin2_index],
								'energy': energy_bins[lgE_idx],
								'Nmu_true': np.mean(Nmu_ground),
			   					'Nmu_pred': np.mean(Nmu_pred),
								'bias': bias,
								'reso': reso}
					Nmu_list.append(Nmu_all)

	df_params = pd.DataFrame(params_list)
	df_Nmu = pd.DataFrame(Nmu_list)
	return df_params, df_Nmu

# fitting and return constant array by fitting all data across primaries
def GetConstFromRegression2(ThreshR, filteringXmax, NmuNorm):
	df = pd.read_parquet(f'RandomForestRegression/data_for_regression_{ThreshR}_Xmax_{filteringXmax}_NmuNorm_{NmuNorm}.parquet')

	# binning energy and zenith angle 
	costheta = df['costheta']
	energy = df['energy']
	particle = df['particle']
	sin2theta = np.sin(np.arccos(costheta))**2
	zenith_bins = np.linspace(0.0, 0.9, 10)
	energy_bins = np.logspace(7.0, 9.0, 21)
	particle_bins = ['proton', 'helium', 'oxygen', 'iron']

	# remove showers in sin2_0.9 bin
	zenith_bins = zenith_bins[:-1]

	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)

	params_list = []
	Nmu_list = []


	for lgE_idx in range(len(energy_bins)):
		for sin2_index in range(len(zenith_bins)):

			
			mask_energy = (energy_indices == lgE_idx +1) # bin 0 is less than  lgE_16.0
			mask_zenith = (zenith_indices == sin2_index +1) # bin 0 is less than sin2_0.0 
			combined_mask = mask_energy & mask_zenith
			
			df_mask = df[combined_mask]
			
			if df_mask.empty: 
				continue
			
			# define parameters
			Edep = df_mask['Edep']
			Erad = df_mask['Erad']
			Xmax = df_mask['Xmax']
			Ne_ground = df_mask['Ne_ground']
			Ne_Xmax = df_mask['Ne']
			Nmu_ground = df_mask['Nmu']
			costheta = df_mask['costheta']
			Ne_ground_Xmax = Ne_ground/Ne_Xmax
			Xv = X0/ costheta
			dXmax = Xv - Xmax
	
			# fitting 
			X_Nmu_Ne = (Ne_ground, Nmu_ground)
			popt, pcov = curve_fit(linearEdep, X_Nmu_Ne, Edep)
			perr = np.sqrt(np.diag(pcov))
			alpha, beta = popt
			alpha_sd, beta_sd = perr[0], perr[1] 

			popt2, pcov2 = curve_fit(linearErad_NeXmax, Ne_Xmax, Erad)
			perr2 = np.sqrt(np.diag(pcov2))
			gamma , b, gamma_sd , b_sd = popt2[0], popt2[1], perr2[0], perr2[1]

			popt3, pcov3 = curve_fit(ExpoNeRatio, dXmax, Ne_ground_Xmax, p0 = [250,1])
			perr3 = np.sqrt(np.diag(pcov3))
			delta, A, delta_sd, A_sd = popt3[0], popt3[1], perr3[0], perr3[1]

			Nmu_pred  = Pred_Nmu(alpha, beta, delta, gamma, A, b, Edep, Erad, dXmax)

			# mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
			# Nmu_pred = Nmu_pred[mask_Nmu_pred]
			# Nmu_ground = Nmu_ground[mask_Nmu_pred]

			R2 = r2_score(Nmu_ground, Nmu_pred)

			residual = Nmu_pred - Nmu_ground
			bias = np.mean(residual)
			reso = np.std(residual)
		

			params_all = {'sin2theta': zenith_bins[sin2_index],
				'energy': energy_bins[lgE_idx],
				'alpha': alpha, 
				'beta': beta, 
				'delta': delta, 
				'gamma': gamma, 
				'A': A, 
				'b': b, 
				'R2': R2,
				'bias': bias,
				'reso': reso}
			params_list.append(params_all)

			for i, particle_bin in enumerate(particle_bins):
				# masking

				mask_particle = (particle == particle_bin)
				mask_energy = (energy_indices == lgE_idx +1) # bin 0 is less than  lgE_16.0
				mask_zenith = (zenith_indices == sin2_index +1) # bin 0 is less than sin2_0.0 
				combined_mask = mask_particle & mask_energy & mask_zenith
				
				df_mask = df[combined_mask]
				
				if df_mask.empty: 
					continue
				
				# define parameters
				Edep = df_mask['Edep']
				Erad = df_mask['Erad']
				Xmax = df_mask['Xmax']
				Nmu_ground = df_mask['Nmu']
				costheta = df_mask['costheta']
				Xv = X0/ costheta
				dXmax = Xv - Xmax
				Nmu_pred = Pred_Nmu(alpha, beta, delta
				, gamma , A, b , Edep, Erad, dXmax)
				
				Nmu_ground = np.log10(Nmu_ground)
				Nmu_pred = np.log10(Nmu_pred)
				
				residual = Nmu_ground - Nmu_pred
				bias = np.mean(residual)
				reso = np.std(residual)

				Nmu_all = {'particle': particle_bin, 
							'sin2theta': zenith_bins[sin2_index],
							'energy': energy_bins[lgE_idx],
							'Nmu_true': np.mean(Nmu_ground), #log10(Nmu_true)
							'Nmu_pred': np.mean(Nmu_pred), #log10(Nmu_pred)
							'SD_Nmu_true': np.std(Nmu_ground),
							'bias': bias,
							'reso': reso}
				Nmu_list.append(Nmu_all)

	df_params = pd.DataFrame(params_list)
	df_Nmu = pd.DataFrame(Nmu_list)
	return df_params, df_Nmu

# constant interpolation
# df_params is from GetConstFromRegression 
def const_interp(df_params, primEnergy, sin2, Edep, Erad,Xmax, plotting, printing):
	# preparing input parameters
	costheta = np.cos(np.arcsin(np.sqrt(sin2)))
	X0 = 697.6
	lgE = np.log10(primEnergy)
	lgEdep = np.log10(Edep)
	lgErad = np.log10(Erad)
	Xv = X0/ costheta
	dXmax = Xv - Xmax

	# preparing data from DataFrame
	const_list = ['alpha', 'beta', 'delta', 'gamma', 'A', 'b', 'bias', 'reso', 'R2']
	const_label = [rf'$\alpha$', rf'$\beta$' , rf'$\delta$' , rf'$\gamma$' , 'A', 'b', 'bias', 'resolution', rf'$\mathrm{{R^2}}$']
	interp_dict ={'alpha':[],
					'beta': [],
					'delta': [],
					'gamma': [],
					'A': [],
					'b': [],
					'bias': [],
					'reso': [],
					'R2':[]}

	for i, const in enumerate(const_list):
		points = df_params[['sin2theta', 'energy']].values
		values = df_params[const].values

		# create grid space
		yi = np.linspace(7.0, 9, 1000)    #  energy_bin range  
		xi = np.linspace(0, 0.8, 1000)    #  sin2_theta range
		X, Y = np.meshgrid(xi, yi)	

		# interpolate grid space
		Z = griddata(points, values, (X, Y), method='cubic')

		if plotting == True:
			plt.figure(figsize=(8, 6))
			plt.imshow(Z, 
					extent=(min(xi), max(xi), min(yi), max(yi)), 
					origin='lower', 
					aspect='auto', 
					cmap='viridis')

			plt.colorbar(label=const_label[i])
			plt.xlabel(rf'$\mathrm{{sin^2 \theta}}$')
			plt.ylabel(rf'log(E(GeV))')
			plt.title(rf'2D Interpolation of {const_label[i]}')
			plt.show()


		interp_const = griddata(points, values, (sin2, lgE), method='cubic')
		interp_dict[const].append(interp_const)
		
	alpha = interp_dict['alpha'][0] # [0] just to evaluate the number 
	beta = interp_dict['beta'][0]
	delta = interp_dict['delta'][0]
	gamma = interp_dict['gamma'][0]
	A = interp_dict['A'][0]
	b = interp_dict['b'][0]


	Nmu_pred  = Pred_Nmu(alpha, beta, delta, gamma, A, b, lgEdep, lgErad, dXmax)
	if printing == True:
		print('Input:')
		print(rf'primary energy {primEnergy:.2e} GeV')
		print(rf'sin2theta {sin2} ')
		print(rf'Edep {Edep} GeV')
		print(rf'Erad {Erad:.2e} GeV')
		print(rf'Xmax {Xmax} g/cm^2')
		print(rf'Predicted logNmu is {Nmu_pred}')
	# return 
	return alpha, beta, delta, gamma, A, b
#%% plotting function

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

def pltRadEAllzenith( primary, energy, filtering, style, corr):
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
				
				fMuon = f'GroundTotalParticles/{p}_{e}_{sin2theta}.npz'
				fXmax = fp_Xmax(p, e, sin2theta)
				fRadE = fp_RadE_norm2(p, e, sin2theta)
				fRadE_unnorm = fp_RadE(p, e, sin2theta)
				fEdep = f'TotalEdepScint/200m/{p}_{e}_{sin2theta}.npz'

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
				Edep = Edep[mask]
				Nmu_ground = Nmu_ground[mask]

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

# Plot deposited energy vs kinetic energy
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


# Plot deposited energy vs primary energy 
def pltEdepEprim(primary, threshR):
	energies = [f"lgE_{x/10:.1f}" for x in range(160,181)]
	sin2theta = [x/10 for x in range(0, 10)]
	s = 5
	alpha = 0.3
	output_dir = f'ScintillatorResponse/'
	output_path = os.path.join(output_dir, f'{primary}.pdf')
	with PdfPages(output_path) as pdf:
		for theta in range(len(sin2theta)):
			fig = plt.figure()
			for i in energies:
				fileScint = f'TotalEdepScint/{threshR}m/{primary}_{i}_{str(sin2theta[theta])}.npz'
				dataScint = np.load(fileScint)

				x = dataScint["primaryE"]
				y1 = dataScint["Edep_e"]
				y2 = dataScint["Edep_mu"]
				y3 = dataScint["Edep_tot"]

				
				sc2 = plt.scatter(x,y2, color = 'yellowgreen', s = s, alpha = alpha)
				sc1 = plt.scatter(x,y1, color= 'orange', s = s, alpha = alpha)
				sc3 = plt.scatter(x,y3, color = 'lightseagreen', s = s, alpha = alpha )

			zenith = np.arcsin(np.sqrt(sin2theta[theta]))

			plt.yscale("log")
			plt.xscale("log")
			plt.xlabel("primary energy (GeV)")
			plt.ylabel(f"$\mathrm{{E_{{deposited}}}}$ (MeV)")
			plt.legend([sc1, sc2, sc3], [f'$e^-$', f'$\mu^{{\pm}}$', 'all particles'])
			plt.title(rf'primary: {primary}, $\mathrm{{sin^2 \theta}}$ = {sin2theta[theta]} ({np.rad2deg(zenith):.2f}$^\circ$)')
			plt.show()

			pdf.savefig(fig, bbox_inches="tight")
	
			plt.close(fig)


def pltSpecificBinFit(sin2_bin, lgE_bin_GeV, particle,filteringXmax,ThreshR,NmuNorm,dataforregression):
	df = pd.read_parquet(f'RandomForestRegression/data_for_regression_{ThreshR}_Xmax_{filteringXmax}_NmuNorm_{NmuNorm}.parquet')
	lgE_bin = lgE_bin_GeV - 9
	colors = {'proton': 'red',
			'helium': 'gold',
			'oxygen': 'green',
			'iron': 'blue'}
	cmap = ['viridis', 'plasma']
	costheta = df['costheta']
	energy = df['energy']
	particle_all = df['particle']
	sin2theta = np.sin(np.arccos(costheta))**2
	zenith_bins = np.linspace(0.0, 0.9, 10)
	energy_bins = np.linspace(7.0, 9.0, 21)
	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)

	X0 = 697.6

	figsize = (7, 5)
	s = 15


	fig1, ax1 = plt.subplots(figsize= figsize, dpi = 150, sharey = True)
	fig2, ax2 = plt.subplots(figsize= figsize ,dpi = 150, sharey = False)
	fig3 = plt.figure(figsize=(8, 8)) 
	ax3 = fig3.add_subplot(projection='3d')
	ax3.set_box_aspect(None, zoom=0.88)
	ax3.zaxis.set_rotate_label(False)
	# fig4, (ax4) = plt.subplots(figsize= figsize ,dpi = 150, sharey = True)

	str_Erad = rf'$\mathrm{{E_{{rad}}}}$ (GeV)/ $\sin^2\alpha$'
	str_Ne_Xmax = rf'$\mathrm{{N_{{e,Xmax}}}}$'
	str_Xmax = rf'$\mathrm{{X_{{max}} (g/cm^2)}} $'
	str_dXmax = rf'$\mathrm{{dX_{{max}} (g/cm^2)}}$'
	str_Ne_ratio = rf'$\mathrm{{N_{{e,ground}} / N_{{e,Xmax}}}}$'   
	str_Ne_ground = rf'$\mathrm{{N_{{e, ground}}}}$'
	str_Nmu_ground = rf'$\mathrm{{N_{{\mu, ground}}}}$'

	for i, p in enumerate(particle):
		fXmax = fp_Xmax(p, f'lgE_{lgE_bin_GeV}', str(sin2_bin))
		fRadE = fp_RadE_norm2(p, f'lgE_{lgE_bin_GeV}', str(sin2_bin))
		fEdep = f'TotalEdepScint/{ThreshR}m/{p}_lgE_{lgE_bin_GeV}_{sin2_bin}.npz'

		fileXmax = np.loadtxt(fXmax)
		fileground = np.load(fp_groundTot(p, f'lgE_{lgE_bin_GeV}', str(sin2_bin)), allow_pickle= True)
		fileRadE = np.load(fRadE, allow_pickle=True)
		fileEdep = np.load(fEdep)

		if dataforregression == False:
			Ne = fileXmax[:,5]
			Xmax = fileXmax[:,6]
			Ne_ground = fileground['nEP']
			Nmu_ground = (fileground['nMu'])
			RadE = (fileRadE['radE_filtered(eV)'])
			Edep = (fileEdep["Edep_tot"])
			zenith = fileground['zenith']
			costheta = np.cos(zenith)
			Xv = X0/ costheta
			combined_mask = (Xmax <= Xv) & (Xmax >= 0)
		
		if dataforregression == True:
			Ne = 10**(df['Ne'])
			Xmax = df['Xmax']
			Ne_ground = 10**(df['Ne_ground'])
			Nmu_ground = 10**(df['Nmu'])
			RadE = 10**(df['Erad'])
			Edep = 10**(df['Edep'])
			costheta = df['costheta']

			sin2_index = np.digitize([sin2_bin], zenith_bins)[0]
			lgE_index = np.digitize([lgE_bin], energy_bins)[0]
			mask_particle = (particle_all == p)
			mask_energy = (energy_indices == lgE_index)
			mask_zenith = (zenith_indices == sin2_index)
			combined_mask = mask_particle  & mask_zenith & mask_energy

		Ne = Ne[combined_mask]
		Xmax = Xmax[combined_mask]
		Ne_ground = Ne_ground[combined_mask]
		Nmu_ground = Nmu_ground[combined_mask]
		RadE = RadE[combined_mask]
		Edep = Edep[combined_mask]
		costheta = costheta[combined_mask]

		Xv = X0/ costheta
		dXmax = Xv - Xmax
		Ne_ground_per_Xmax = Ne_ground/Ne


		##################### Fitting #####################

		X_Nmu_Ne = (Ne_ground, Nmu_ground)
		popt, pcov = curve_fit(linearEdep, X_Nmu_Ne, Edep)
		alpha, beta = popt

		popt2, pcov2 = curve_fit(linearErad_NeXmax, Ne, RadE)
		gamma , b = popt2[0], popt2[1]

		popt3, pcov3 = curve_fit(ExpoNeRatio, dXmax, Ne_ground_per_Xmax, p0 = [250,1])
		delta, A = popt3[0], popt3[1]
		# delta = popt3[0]

		print(rf"alpha = {alpha:.4f}, beta = {beta:.4f}, gamma = {gamma:.2e}, b = {b:.2e}, delta = {delta:.2e}")

		##################### Error Calculation #####################
		res1 = RadE - linearErad_NeXmax(Ne, *popt2)
		dof1 = len(Ne) - len(popt2)
		rmse = np.sqrt(np.sum(res1**2) / dof1)

		res2 = Ne_ground_per_Xmax - ExpoNeRatio(dXmax, *popt3)
		dof2 = len(dXmax) - len(popt3)
		rmse2 = np.sqrt(np.sum(res2**2) / dof2)

		res3 = Edep - linearEdep(X_Nmu_Ne, *popt)
		dof3 = len(X_Nmu_Ne*2) - len(popt)
		rmse3 = np.sqrt(np.sum(res3**2) / dof3)


		##################### Plotting #####################
		x1 = np.linspace(min(Ne), max(Ne), num = len(Ne))
		x2 = np.linspace(min(dXmax), max(dXmax), num = len(dXmax))

		ax1.errorbar(Ne, RadE, yerr=rmse, fmt='o', color = colors[p], label = p, zorder = 2, elinewidth=1, ms = s/5)
		ax1.plot(x1, linearErad_NeXmax(x1, *popt2), color = 'k', zorder = 3)
		# ax1.set_title(f'{p} {e} sin2_{sin2theta} all data')
		ax1.set_title(f'lgE_{lgE_bin_GeV} sin2_{sin2_bin}')
		# ax1.text(min(Ne), max(RadE)*0.9, rf'$E_{{rad}} = \gamma N_{{e,Xmax}} + b$' +'\n' + rf'$\gamma$ = {gamma:.2e}' +'\n' + rf'b = {b:.2e}')
		ax1.set_ylabel(str_Erad)
		ax1.set_xlabel(str_Ne_Xmax)

		ax2.errorbar(dXmax, Ne_ground_per_Xmax, yerr=rmse2, fmt='o', color = colors[p], label = p, zorder = 2, elinewidth=1, ms = s/5)
		# ax2.text(400, max(Ne_ground_per_Xmax)*0.9, rf'$N_{{e,ground}}/N_{{e,Xmax}} = A* e^{{-dXmax/ \delta}}$' +'\n' + rf'$\delta$ = {delta:.2e}' +'\n' + rf'A = {A:.2e}')
		ax2.plot(x2, ExpoNeRatio(x2, *popt3), color = 'k', zorder = 3)
		ax2.set_title(f'lgE_{lgE_bin_GeV} sin2_{sin2_bin}')
		ax2.set_ylabel(str_Ne_ratio)
		ax2.set_xlabel(str_dXmax)

		ax3.errorbar(Nmu_ground, Ne_ground, Edep, zerr=rmse3, fmt='o', color = colors[p], ms = s/5, label = p)
		ax3.set_title(f'lgE_{lgE_bin_GeV} sin2_{sin2_bin}')
		# ax3.text(min(Ne_ground), max(Nmu_ground)*0.9, rf'$E_{{dep}} = \alpha N_e + \beta N_{{\mu}}$' +'\n' + rf'$\alpha$ = {alpha:.2e}' +'\n' + rf'$\beta$ = {beta:.2e}')
		ax3.set_xlabel(str_Nmu_ground)
		ax3.set_ylabel(str_Ne_ground)
		ax3.set_zlabel('Edep (GeV)', labelpad=10, rotation=90)

		if p == 'proton':
			ax1.plot(x1, linearErad_NeXmax(x1, *popt2), color = 'k', zorder = 3, label = 'linear fit')
			ax2.plot(x2, ExpoNeRatio(x2, *popt3), color = 'k', zorder = 3, label = 'exp fit')
		# vmin2, vmax2 = min(Ne_ground), max(Ne_ground)
		# sc3 = ax4.scatter(Nmu_ground, Edep, c=Ne_ground, cmap=cmap[i], s = s, vmin=vmin2, vmax=vmax2, label = p)
		# cbar = fig4.colorbar(sc3, ax=[ax4])
		# cbar.set_label(str_Ne_ground)
		# ax4.set_title(f'lgE_{lgE_bin_GeV} sin2_{sin2_bin}')
		# ax4.set_ylabel('Edep (GeV)')
		# ax4.set_ylabel('Edep (GeV)')
		# ax4.set_yscale('log')
		# ax4.set_xscale('log')
		# ax4.set_xlabel(str_Nmu_ground)

	ax1.legend()
	ax2.legend()
	ax3.legend()
	# ax4.legend()

#%% plot the resulution og the reconstruction by regression across the primary energies
def plt_recon_res (df_Nmu, filteringXmax, ThreshR, sin2):
	energy_bins = np.logspace(7.0, 9.0, 21)
	sin2theta = df_Nmu['sin2theta']
	energy = df_Nmu['energy']
	particle = df_Nmu['particle']

	fig, ax = plt.subplots(figsize=(8,5), constrained_layout=True)
	fig2, ax2 = plt.subplots(figsize=(8,5), constrained_layout=True)
	fig3, ax3 = plt.subplots(figsize =(8,5), constrained_layout = True)

	handle_bias = []
	handle_res = []
	handle_true = []
	handle_recon = []
	colors = {'proton': 'red',
				'helium': 'gold',
				'oxygen': 'green',
				'iron': 'blue'}
	particle_bins = ['proton', 'helium', 'oxygen', 'iron']

	titlelabel =  f'removing underground Xmax = {filteringXmax}, R > {ThreshR} m'
	zenith_label = rf'$\mathrm{{sin^2 \theta = {sin2}}}$'
	x_label = rf'$\mathrm{{log}}(E_0 \mathrm{{(eV)}})$'


	for p_idx, p in enumerate(particle_bins):

		Nmu_true_all = []
		Nmu_pred_all = []
		SD_Nmu_true_all = []
		bias_all = []
		reso_all = []
		energy_all = []

		for e_idx, e in enumerate(energy_bins):

			mask = np.isclose(sin2theta, sin2) & (energy == e) & (particle == p)
			df_mask = df_Nmu[mask]
			Nmu_true = df_mask['Nmu_true'].values
			Nmu_pred = df_mask['Nmu_pred'].values
			SD_Nmu_true = df_mask['SD_Nmu_true'].values
			bias = df_mask['bias'].values
			reso = df_mask['reso'].values
			
			if not df_mask.empty:
				Nmu_true_all.append(Nmu_true)
				Nmu_pred_all.append(Nmu_pred)
				SD_Nmu_true_all.append(SD_Nmu_true)
				bias_all.append(bias)
				reso_all.append(reso)
				energy_all.append(e)
				
		energy_all = np.log10(energy_all) + 9
		bias = ax.scatter(energy_all, bias_all, color = colors[p], marker = 'x', label = p )
		res = ax.scatter(energy_all, reso_all, color = colors[p], marker = '^', label = p)
		
		true = ax2.scatter(energy_all, Nmu_true_all, color = colors[p], marker = 'x', label = p, s = 50, zorder = 2)
		recon = ax2.scatter(energy_all, Nmu_pred_all, color = colors[p], marker = 'o', label = p, s = 20, zorder = 3,
						edgecolors = 'black', linewidths=0.5, alpha = 0.7)

		SD_true = ax3.scatter(energy_all, SD_Nmu_true_all,color = colors[p], marker = 'o', label = p)

		handle_res.append(res)
		handle_bias.append(bias)
		handle_true.append(true)
		handle_recon.append(recon)

	fig.suptitle(titlelabel + '\n' + zenith_label)
	leg1 = ax.legend(bbox_to_anchor=(1.01, 1), handles= handle_res, loc='upper left', title="Resolution")
	ax.add_artist(leg1) 
	ax.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_bias, loc='upper left', title="Bias")
	ax.set_xlabel(x_label)
	ax.set_ylabel(rf'$\mathrm{{log(N_{{\mu^{{\pm}}}}^{{true}}) - log(N_{{\mu^{{\pm}}}}^{{recon}})}}$')
	ax.set_title(f'')
	ax.set_ylim(-0.23,0.23)
	ax.hlines(y= 0, xmin = min(energy_all), xmax= max(energy_all), colors = 'k', linestyles =  '--')
	ax.minorticks_on()

	fig2.suptitle(titlelabel + '\n' + zenith_label)
	leg2 = ax2.legend(bbox_to_anchor=(1.01, 1), handles= handle_true, loc='upper left', title="True Value")
	ax2.add_artist(leg2) 
	ax2.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_recon, loc='upper left', title="Reconstructed Value")
	ax2.set_xlabel(x_label)
	ax2.set_ylabel(rf'$\mathrm{{log(N_{{\mu^{{\pm}}}})}}$')
	ax2.minorticks_on()

	fig3.suptitle(titlelabel + '\n' + zenith_label)
	ax3.legend()
	ax3.set_xlabel(x_label)
	ax3.set_ylim(0,0.15)
	ax3.set_ylabel(rf'$\mathrm{{\sigma_{{log(N_{{\mu}}^{{true}})}}}}$')
	ax3.minorticks_on()