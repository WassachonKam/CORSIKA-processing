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
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import textwrap

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
	energy_bins = np.logspace(7.0, 9.0, 21)
	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)



	particle_bins = ['proton', 'helium', 'oxygen', 'iron']

	colors = ['red','gold', 'green', 'blue']
	lgE_bin = 10**(lgE_GeV_bin - 9)

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

			Nmu_ground = np.log10(Nmu_ground)
			Nmu_pred = np.log10(Nmu_pred)

			# mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
			# Nmu_pred = Nmu_pred[mask_Nmu_pred]
			# Nmu_ground = Nmu_ground[mask_Nmu_pred]

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
	energy_bins = np.logspace(7.0, 9.0, 21)

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

			Nmu_ground = np.log10(Nmu_ground)
			Nmu_pred = np.log10(Nmu_pred)

			# mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
			# Nmu_pred = Nmu_pred[mask_Nmu_pred]
			# Nmu_ground = Nmu_ground[mask_Nmu_pred]

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


				# mask_Nmu_pred = (Nmu_pred < 10) & (Nmu_pred > 0)
				# Nmu_pred = Nmu_pred[mask_Nmu_pred]
				# Nmu_ground = Nmu_ground[mask_Nmu_pred]

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
	Nmu_list_events = []


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

				for j in range(len(df_mask)):
					Nmu_all_events = {
						'particle': particle_bin,
						'sin2theta': zenith_bins[sin2_index],
						'energy': energy_bins[lgE_idx],
						'Nmu_true': Nmu_ground.iloc[j],
						'Nmu_pred': Nmu_pred.iloc[j],
						'residual': residual.iloc[j],
					}
					Nmu_list_events.append(Nmu_all_events)

	df_params = pd.DataFrame(params_list)
	df_Nmu = pd.DataFrame(Nmu_list)
	df_Nmu_all_events = pd.DataFrame(Nmu_list_events)
	return df_params, df_Nmu, df_Nmu_all_events

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
		
		df_params['log_energy'] = np.log10(df_params['energy'])
		points = df_params[['sin2theta', 'log_energy']].values
		values = df_params[const].values

		# create grid space
		yi = np.logspace(7.0, 9, 1000)    #  energy_bin range  
		xi = np.linspace(0, 0.8, 1000)    #  sin2_theta range

		yi = np.log10(yi)
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


def pltSpecificBinFit1(sin2_bin, lgE_bin_GeV, particle,filteringXmax,ThreshR,NmuNorm,dataforregression):
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

	str_Erad = rf'Normalized $\mathrm{{E_{{rad}}}}$ (GeV)'
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
		dof3 = len(Edep) - len(popt)
		rmse3 = np.sqrt(np.sum(res3**2) / dof3)


		##################### Plotting #####################
		x1 = np.linspace(min(Ne), max(Ne), num = len(Ne))
		x2 = np.linspace(min(dXmax), max(dXmax), num = len(dXmax))

		ax1.errorbar(Ne, RadE, yerr=rmse, fmt='o', color = colors[p], zorder = 2, elinewidth=1, ms = s/5)
		ax1.scatter(Ne, RadE, color = colors[p], label = p, s = s)
		# ax1.plot(x1, linearErad_NeXmax(x1, *popt2), color = 'k', zorder = 3) # regression line
		# ax1.set_title(f'{p} {e} sin2_{sin2theta} all data') 
		ax1.set_title(f'lgE_{lgE_bin_GeV} sin2_{sin2_bin}')
		# ax1.text(min(Ne), max(RadE)*0.9, rf'$E_{{rad}} = \gamma N_{{e,Xmax}} + b$' +'\n' + rf'$\gamma$ = {gamma:.2e}' +'\n' + rf'b = {b:.2e}')
		ax1.set_ylabel(str_Erad)
		ax1.set_xlabel(str_Ne_Xmax)

		ax2.errorbar(dXmax, Ne_ground_per_Xmax, yerr=rmse2, fmt='o', color = colors[p], zorder = 2, elinewidth=1, ms = s/5)
		ax2.scatter(dXmax, Ne_ground_per_Xmax, color = colors[p], label = p, s = s)
		# ax2.text(400, max(Ne_ground_per_Xmax)*0.9, rf'$N_{{e,ground}}/N_{{e,Xmax}} = A* e^{{-dXmax/ \delta}}$' +'\n' + rf'$\delta$ = {delta:.2e}' +'\n' + rf'A = {A:.2e}')
		# ax2.plot(x2, ExpoNeRatio(x2, *popt3), color = 'k', zorder = 3) # regression line
		ax2.set_title(f'lgE_{lgE_bin_GeV} sin2_{sin2_bin}')
		ax2.set_ylabel(str_Ne_ratio)
		ax2.set_xlabel(str_dXmax)

		ax3.errorbar(Nmu_ground, Ne_ground, Edep, zerr=rmse3, fmt='o', color = colors[p], ms = s/5, alpha = 0.2)
		ax3.scatter(Nmu_ground, Ne_ground, Edep, color = colors[p], label = p, s = s)
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
		# variance = [x ** 2 for x in reso_all]
		reso_avg = np.sqrt(np.mean(np.square(reso_all)))

		handle_res.append(res)
		handle_bias.append(bias)
		handle_true.append(true)
		handle_recon.append(recon)

		# print(p, np.mean(bias_all), np.sqrt(np.sum(variance))/len(reso_all))
		print(p, np.mean(bias_all), reso_avg)

	fig.suptitle(titlelabel + '\n' + zenith_label)
	leg1 = ax.legend(bbox_to_anchor=(1.01, 1), handles= handle_res, loc='upper left', title="Resolution")
	ax.add_artist(leg1) 
	ax.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_bias, loc='upper left', title="Bias")
	ax.set_xlabel(x_label)
	ax.set_ylabel(rf'$\mathrm{{log}}(N_{{\mu^{{\pm}}}}^{{true}}) - \mathrm{{log}}(N_{{\mu^{{\pm}}}}^{{recon}})$')
	ax.set_title(f'')
	ax.set_ylim(-0.23,0.23)
	ax.hlines(y= 0, xmin = min(energy_all), xmax= max(energy_all), colors = 'k', linestyles =  '--')
	ax.minorticks_on()

	fig2.suptitle(titlelabel + '\n' + zenith_label)
	leg2 = ax2.legend(bbox_to_anchor=(1.01, 1), handles= handle_true, loc='upper left', title="True Value")
	ax2.add_artist(leg2) 
	ax2.legend(bbox_to_anchor=(1.01, 0.6), handles=handle_recon, loc='upper left', title="Reconstructed Value")
	ax2.set_xlabel(x_label)
	ax2.set_ylabel(rf'$\mathrm{{log}}(N_{{\mu^{{\pm}}}})$')
	ax2.minorticks_on()

	fig3.suptitle(titlelabel + '\n' + zenith_label)
	ax3.legend()
	ax3.set_xlabel(x_label)
	ax3.set_ylim(0,0.15)
	ax3.set_ylabel(rf'$\mathrm{{\sigma_{{log(N_{{\mu}}^{{true}})}}}}$')
	ax3.minorticks_on()




#%% main plot

# ---------------------------------------------------------------- style tokens
SURF = '#fcfcfb'      # chart surface
INK = '#0b0b0b'       # primary text
INK2 = '#52514e'      # secondary text (ticks)
MUTED = '#8a8983'     # annotations
GRID = '#e8e7e3'      # gridlines
CONTEXT = '#d8d7d2'   # the other primaries, recessive
ACCENT = '#2a78d6'    # this panel's primary
BAND = '#9ec5f4'      # +/-1 RMSE ribbon
FITLINE = '#104281'   # fit curve
WARM = '#e34948'      # diverging warm pole (residual below the fit)

# Conventional mass-group colours. Each is used as ONE panel's accent against the
# grey context, never against the other three, so they are identity labels rather
# than a palette that has to self-separate.
PRIMARY_COLORS = {
	'proton': '#ff0000',   # red
	'helium': '#ffd700',   # gold
	'oxygen': '#008000',   # green
	'iron':   '#0000ff',   # blue
}


def _shade(c, amt):
	"""Blend a hex colour toward white (amt > 0) or black (amt < 0)."""
	c = c.lstrip('#')
	rgb = [int(c[i:i + 2], 16) for i in (0, 2, 4)]
	t = 255 if amt > 0 else 0
	f = abs(amt)
	return '#%02x%02x%02x' % tuple(round(v + (t - v) * f) for v in rgb)


BASE_FONT = 9.0   # every other size in the figure is a multiple of this


def _vizstyle(scale=1.0):
	"""rcParams for a recessive, print-safe scientific figure.

	`scale` multiplies every text size at once. The helpers below read their
	sizes back out of rcParams rather than hard-coding numbers, so one scale
	moves titles, labels, ticks and annotations together.
	"""
	b = BASE_FONT * scale
	return {
		'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
		'font.size': b, 'axes.titlesize': b * 10 / 9, 'axes.labelsize': b,
		'figure.titlesize': b * 12 / 9,
		'axes.edgecolor': '#c9c8c3', 'axes.linewidth': 0.8,
		'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
		'xtick.labelsize': b * 8 / 9, 'ytick.labelsize': b * 8 / 9,
		'legend.frameon': False, 'figure.dpi': 500,
	}


def _fs(k=1.0):
	"""A size relative to the current base font."""
	return plt.rcParams['font.size'] * k


def _pow10(arr):
	"""Exponent to divide out so the axis carries no corner offset box."""
	m = np.nanmax(np.abs(np.asarray(arr, dtype=float)))
	if not np.isfinite(m) or m == 0:
		return 0
	e = int(np.floor(np.log10(m)))
	return e if abs(e) >= 3 else 0


def _label(base, e):
	return base if e == 0 else rf'{base}  ($10^{{{e}}}$)'


def _subtitle(fig, text, y, width=96):
	"""Place the subtitle, wrapped, and record how much vertical room it took.

	Naming the pooled fit takes more characters than a one-line subtitle can
	hold at this figure width, so wrap rather than let it run off the canvas.
	"""
	# bigger type fits fewer characters per line and eats more vertical room,
	# so both the wrap width and the reserved height track the font size
	k = plt.rcParams['font.size'] / BASE_FONT
	lines = textwrap.wrap(text, width=max(24, int(width / k))) or ['']
	fig.text(0.5, y, '\n'.join(lines), ha='center', va='top', fontsize=_fs(),
	         color=MUTED, linespacing=1.4)
	fig._viz_top = y - 0.022 * k * len(lines)   # read back for tight_layout's rect
	return fig._viz_top


def _sub(what, scope, rmse, params, alt='dashed grey = that primary’s own fit'):
	"""Subtitle naming which fit is drawn -- a pooled fit must say so, or a
	reader assumes each panel was fitted on its own points."""
	if scope == 'per-primary':
		return f'{what} fit per primary; band is ' + u'±' + '1 RMSE about that fit'
	tail = '' if rmse is None else f'; band is ±1 RMSE = {rmse:.4g}'
	head = f'one {what} fit to all primaries pooled ({params}){tail}'
	if scope == 'both':
		head += f'; {alt}'
	return head


def _res_rmse(res, npar):
	"""RMSE with the dof taken from the residual array.

	Never from the input container: len() on an (Ne, Nmu) tuple returns 2 and
	inflates the result by sqrt(N/2) without raising.
	"""
	return np.sqrt(np.sum(res ** 2) / (res.size - npar))


def _facets(n, title, subtitle, figsize=(8.2, 6.6)):
	"""2-column grid of shared-scale panels with a title block."""
	nrow = int(np.ceil(n / 2))
	fig, axes = plt.subplots(nrow, 2, figsize=figsize, sharex=True, sharey=True,
	                         squeeze=False)
	# mathtext in the title is taller than its font size -- keep the two lines
	# well apart or the subscripts collide with the subtitle.
	fig.suptitle(title, fontsize=plt.rcParams['figure.titlesize'], y=0.985, va='top')
	_subtitle(fig, subtitle, 0.923)
	for ax in axes.flat:
		ax.grid(True, color=GRID, lw=0.7, zorder=0)
		ax.set_axisbelow(True)
		for side in ('top', 'right'):
			ax.spines[side].set_visible(False)
	for ax in axes.flat[n:]:          # blank any unused cell
		ax.set_visible(False)
	return fig, axes


def _facets3d(n, title, subtitle, figsize=(10.5, 9.0), elev=20, azim=-58):
	"""2-column grid of 3-D panels, styled to match the 2-D facets."""
	nrow = int(np.ceil(n / 2))
	fig = plt.figure(figsize=figsize)
	fig.suptitle(title, fontsize=plt.rcParams['figure.titlesize'], y=0.985, va='top')
	_subtitle(fig, subtitle, 0.945, width=104)
	axes = []
	for i in range(n):
		ax = fig.add_subplot(nrow, 2, i + 1, projection='3d')
		ax.view_init(elev=elev, azim=azim)
		ax.set_box_aspect(None, zoom=0.94)
		ax.zaxis.set_rotate_label(False)   # keep the z label upright
		# default 3-D panes are a heavy grey box; lighten them to the 2-D tokens
		for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
			axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
			axis.pane.set_edgecolor(GRID)
			axis._axinfo['grid'].update(color=GRID, linewidth=0.6)
			axis.line.set_color('#c9c8c3')
		ax.tick_params(labelsize=plt.rcParams['xtick.labelsize'], pad=1.5)
		axes.append(ax)
	fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=fig._viz_top,
	                    wspace=0.0, hspace=0.0)
	return fig, axes


def _note(ax, text, corner='lower right'):
	"""Fit parameters, on a surface-coloured plate so data never shows through."""
	va, ha = corner.split()
	ax.text(0.97 if ha == 'right' else 0.03, 0.05 if va == 'lower' else 0.95,
	        text, transform=ax.transAxes, ha=ha, va='bottom' if va == 'lower' else 'top',
	        fontsize=_fs(8 / 9), color=MUTED, linespacing=1.5, zorder=6,
	        bbox=dict(facecolor=SURF, edgecolor='none', alpha=0.78, pad=2.5))


def _panel(ax, x, y, ctx_x, ctx_y, model, popt, rmse, name, note,
           corner='lower right', alt_popt=None, accent=ACCENT):
	"""One primary in accent, the rest as grey context, fit + RMSE ribbon.

	`popt` is the fit the ribbon belongs to. `alt_popt`, when given, is drawn as
	a thin dashed line for comparison -- the per-primary fit next to the pooled
	one, so any mass dependence shows as a divergence between the two.
	"""
	ax.scatter(ctx_x, ctx_y, s=9, color=CONTEXT, lw=0, zorder=1, rasterized=True)

	# band and fit line are shades of the panel's own colour, so the three marks
	# read as one series. gold needs the darkening most -- at 1.4:1 on white an
	# undarkened fit line is invisible.
	band, fitline = _shade(accent, 0.62), _shade(accent, -0.45)

	# fit drawn only across the range this primary actually covers
	xf = np.linspace(np.min(x), np.max(x), 200)
	yf = model(xf, *popt)
	ax.fill_between(xf, yf - rmse, yf + rmse, color=band, alpha=0.5, lw=0, zorder=2)
	ax.plot(xf, yf, color=fitline, lw=2, zorder=4)
	if alt_popt is not None:
		# neutral, not WARM -- red is proton's identity colour here
		ax.plot(xf, model(xf, *alt_popt), color=INK2, lw=1.3, ls=(0, (4, 2)),
		        zorder=5)

	ax.scatter(x, y, s=14, color=accent, lw=0.4, edgecolor=_shade(accent, -0.35),
	           zorder=3, rasterized=True)

	ax.set_title(name, loc='left', fontweight='bold', color=INK, pad=6)
	_note(ax, note, corner)


# ------------------------------------------------------------------- the plot
def pltSpecificBinFit(sin2_bin, lgE_bin_GeV, particle, filteringXmax, ThreshR,
                      NmuNorm, dataforregression, show_3d=False,
                      fit_scope='per-primary', colors=None, font_scale=1.0,
                      figsize=None, figsize3d=None):
	"""fit_scope: 'per-primary' fits each primary separately (the original
	behaviour); 'global' fits one relation to all primaries pooled and draws
	that in every panel; 'both' draws the pooled fit solid with the primary's
	own fit dashed over it.

	colors: {primary: hex} overriding PRIMARY_COLORS for the panel accents.
	font_scale: multiplies every text size (1.3 is a good poster/slide value).
	figsize / figsize3d: enlarge the canvas if bigger type starts to crowd."""
	colors = dict(PRIMARY_COLORS, **(colors or {}))
	df = pd.read_parquet(
		f'RandomForestRegression/data_for_regression_{ThreshR}_Xmax_{filteringXmax}_NmuNorm_{NmuNorm}.parquet')
	lgE_bin = lgE_bin_GeV - 9
	costheta_all = df['costheta']
	energy = df['energy']
	particle_all = df['particle']
	sin2theta = np.sin(np.arccos(costheta_all)) ** 2
	zenith_bins = np.linspace(0.0, 0.9, 10)
	energy_bins = np.linspace(7.0, 9.0, 21)
	zenith_indices = np.digitize(sin2theta, zenith_bins)
	energy_indices = np.digitize(energy, energy_bins)

	X0 = 697.6

	# ---------------------------------------------------------------- pass 1
	# load, mask and fit every primary first, so each panel can draw the others
	# as context and so the fits are available before any axes exist.
	D = {}
	for p in particle:
		fXmax = fp_Xmax(p, f'lgE_{lgE_bin_GeV}', str(sin2_bin))
		fRadE = fp_RadE_norm2(p, f'lgE_{lgE_bin_GeV}', str(sin2_bin))
		fEdep = f'TotalEdepScint/{ThreshR}m/{p}_lgE_{lgE_bin_GeV}_{sin2_bin}.npz'

		fileXmax = np.loadtxt(fXmax)
		fileground = np.load(fp_groundTot(p, f'lgE_{lgE_bin_GeV}', str(sin2_bin)),
		                     allow_pickle=True)
		fileRadE = np.load(fRadE, allow_pickle=True)
		fileEdep = np.load(fEdep)

		if dataforregression:
			Ne = 10 ** (df['Ne'])
			Xmax = df['Xmax']
			Ne_ground = 10 ** (df['Ne_ground'])
			Nmu_ground = 10 ** (df['Nmu'])
			RadE = 10 ** (df['Erad'])
			Edep = 10 ** (df['Edep'])
			costheta = df['costheta']

			sin2_index = np.digitize([sin2_bin], zenith_bins)[0]
			lgE_index = np.digitize([lgE_bin], energy_bins)[0]
			combined_mask = ((particle_all == p) & (zenith_indices == sin2_index)
			                 & (energy_indices == lgE_index))
		else:
			Ne = fileXmax[:, 5]
			Xmax = fileXmax[:, 6]
			Ne_ground = fileground['nEP']
			Nmu_ground = fileground['nMu']
			RadE = fileRadE['radE_filtered(eV)']
			Edep = fileEdep['Edep_tot']
			costheta = np.cos(fileground['zenith'])
			combined_mask = (Xmax <= X0 / costheta) & (Xmax >= 0)

		# np.asarray keeps the two branches (numpy / pandas) interchangeable
		# downstream -- curve_fit and the plotting helpers both want arrays.
		sel = lambda a: np.asarray(a)[np.asarray(combined_mask)]
		Ne, Xmax, Ne_ground = sel(Ne), sel(Xmax), sel(Ne_ground)
		Nmu_ground, RadE, Edep, costheta = sel(Nmu_ground), sel(RadE), sel(Edep), sel(costheta)

		dXmax = X0 / costheta - Xmax
		Ne_ground_per_Xmax = Ne_ground / Ne

		# ---- fits
		X_Nmu_Ne = (Ne_ground, Nmu_ground)
		popt, _ = curve_fit(linearEdep, X_Nmu_Ne, Edep)
		popt2, _ = curve_fit(linearErad_NeXmax, Ne, RadE)
		popt3, _ = curve_fit(ExpoNeRatio, dXmax, Ne_ground_per_Xmax, p0=[250, 1])
		alpha, beta = popt
		gamma, b = popt2
		delta, A = popt3

		rmse = _res_rmse(RadE - linearErad_NeXmax(Ne, *popt2), len(popt2))
		rmse2 = _res_rmse(Ne_ground_per_Xmax - ExpoNeRatio(dXmax, *popt3), len(popt3))
		rmse3 = _res_rmse(Edep - linearEdep(X_Nmu_Ne, *popt), len(popt))

		print(f'{p:>7}  alpha = {alpha:.4f}, beta = {beta:.4f}, '
		      f'gamma = {gamma:.2e}, b = {b:.2e}, delta = {delta:.2e}')

		D[p] = dict(Ne=Ne, RadE=RadE, dXmax=dXmax, ratio=Ne_ground_per_Xmax,
		            Ne_ground=Ne_ground, Nmu_ground=Nmu_ground, Edep=Edep,
		            popt=popt, popt2=popt2, popt3=popt3,
		            rmse=rmse, rmse2=rmse2, rmse3=rmse3,
		            alpha=alpha, beta=beta, gamma=gamma, b=b, delta=delta, A=A)

	cat = lambda k: np.concatenate([D[p][k] for p in particle])
	bin_str = rf'lgE {lgE_bin_GeV},  $\sin^2\theta$ = {sin2_bin}'

	# ------------------------------------------------- pooled ("global") fits
	# One relation fitted to every primary at once. Comparing a primary's
	# scatter about THIS fit with its scatter about its own fit is the test of
	# whether the relation is mass-independent.
	G = {}
	G['popt'], _ = curve_fit(linearEdep, (cat('Ne_ground'), cat('Nmu_ground')),
	                         cat('Edep'))
	G['popt2'], _ = curve_fit(linearErad_NeXmax, cat('Ne'), cat('RadE'))
	G['popt3'], _ = curve_fit(ExpoNeRatio, cat('dXmax'), cat('ratio'), p0=[250, 1])
	G['rmse'] = _res_rmse(cat('RadE') - linearErad_NeXmax(cat('Ne'), *G['popt2']),
	                      len(G['popt2']))
	G['rmse2'] = _res_rmse(cat('ratio') - ExpoNeRatio(cat('dXmax'), *G['popt3']),
	                       len(G['popt3']))
	G['rmse3'] = _res_rmse(
		cat('Edep') - linearEdep((cat('Ne_ground'), cat('Nmu_ground')), *G['popt']),
		len(G['popt']))
	G['gamma'], G['b'] = G['popt2']
	G['delta'], G['A'] = G['popt3']
	G['alpha'], G['beta'] = G['popt']

	# Per-primary residuals ABOUT THE POOLED FIT. The mean is the interesting
	# number: a non-zero bias is exactly the mass dependence the pooled fit
	# cannot absorb, and it is invisible when every primary gets its own fit.
	for p in particle:
		d = D[p]
		d['bias'] = np.mean(d['RadE'] - linearErad_NeXmax(d['Ne'], *G['popt2']))
		d['bias2'] = np.mean(d['ratio'] - ExpoNeRatio(d['dXmax'], *G['popt3']))
		d['bias3'] = np.mean(d['Edep']
		                     - linearEdep((d['Ne_ground'], d['Nmu_ground']), *G['popt']))
		d['grmse'] = _res_rmse(d['RadE'] - linearErad_NeXmax(d['Ne'], *G['popt2']), 0)
		d['grmse2'] = _res_rmse(d['ratio'] - ExpoNeRatio(d['dXmax'], *G['popt3']), 0)
		d['grmse3'] = _res_rmse(
			d['Edep'] - linearEdep((d['Ne_ground'], d['Nmu_ground']), *G['popt']), 0)

	if fit_scope not in ('per-primary', 'global', 'both'):
		raise ValueError("fit_scope must be 'per-primary', 'global' or 'both'")
	glob = fit_scope in ('global', 'both')
	F = G if glob else None          # which fit the ribbon belongs to
	print(f'  pooled  alpha = {G["alpha"]:.4f}, beta = {G["beta"]:.4f}, '
	      f'gamma = {G["gamma"]:.2e}, b = {G["b"]:.2e}, delta = {G["delta"]:.2e}')

	str_Erad = r'Normalized $\mathrm{E_{rad}}$ (GeV)'
	str_Ne_Xmax = r'$\mathrm{N_{e,Xmax}}$'
	str_dXmax = r'$\mathrm{dX_{max}}$ (g/cm$^2$)'
	str_Ne_ratio = r'$\mathrm{N_{e,ground} / N_{e,Xmax}}$'
	str_Ne_ground = r'$\mathrm{N_{e, ground}}$'
	str_Nmu_ground = r'$\mathrm{N_{\mu, ground}}$'

	figs = {}
	with plt.rc_context(_vizstyle(font_scale)):
		# ------------------------------------------------ fig1: Erad vs Ne,Xmax
		ex, ey = _pow10(cat('Ne')), _pow10(cat('RadE'))
		sx, sy = 10.0 ** ex, 10.0 ** ey
		fig1, axes = _facets(len(particle),
		                     rf'Normalized $E_{{rad}}$ vs $N_{{e,Xmax}}$    {bin_str}',
		                     _sub('linear', fit_scope, G['rmse'] / sy,
		                          rf'$\gamma$ = {G["gamma"]:.2e}, b = {G["b"]:.2e}'),
		                     **({'figsize': figsize} if figsize else {}))
		for ax, p in zip(axes.flat, particle):
			d = D[p]
			# if glob:
			# 	note = (f'offset = {d["bias"] / sy:+.3f}' + '\n'
			# 	        + f'RMSE = {d["grmse"] / sy:.3f}   n = {d["Ne"].size}'
			# 	        + '\n' + rf'(own fit: $\gamma$ = {d["gamma"]:.2e})')
			# else:
			# 	note = (rf'$\gamma$ = {d["gamma"]:.2e}' + '\n' + rf'b = {d["b"]:.2e}'
			# 	        + '\n' + f'RMSE = {d["rmse"] / sy:.3f}   n = {d["Ne"].size}')
			note = None
			_panel(ax, d['Ne'] / sx, d['RadE'] / sy, cat('Ne') / sx, cat('RadE') / sy,
			       lambda t, *q: linearErad_NeXmax(t * sx, *q) / sy,
			       (F or d)['popt2'], (F or d)['rmse'] / sy, p, note,
			       alt_popt=d['popt2'] if fit_scope == 'both' else None,
			       accent=colors[p])
		fig1.supxlabel(_label(str_Ne_Xmax, ex), fontsize=_fs())
		fig1.supylabel(_label(str_Erad, ey), fontsize=_fs())
		fig1.tight_layout(rect=[0, 0, 1, fig1._viz_top])
		figs['erad'] = fig1

		# ------------------------------------------- fig2: Ne ratio vs dXmax
		fig2, axes = _facets(len(particle),
		                     rf'$N_{{e,ground}}/N_{{e,Xmax}}$ vs $dX_{{max}}$    {bin_str}',
		                     _sub(r'$A\,e^{-dX_{max}/\delta}$', fit_scope, G['rmse2'],
		                          rf'$\delta$ = {G["delta"]:.2e}, A = {G["A"]:.2e}'))
		for ax, p in zip(axes.flat, particle):
			d = D[p]
			# if glob:
			# 	note = (f'offset = {d["bias2"]:+.4f}' + '\n'
			# 	        + f'RMSE = {d["grmse2"]:.4f}   n = {d["dXmax"].size}'
			# 	        + '\n' + rf'(own fit: $\delta$ = {d["delta"]:.2e})')
			# else:
			# 	note = (rf'$\delta$ = {d["delta"]:.2e}' + '\n' + rf'A = {d["A"]:.2e}'
			# 	        + '\n' + f'RMSE = {d["rmse2"]:.4f}   n = {d["dXmax"].size}')
			note = None
			_panel(ax, d['dXmax'], d['ratio'], cat('dXmax'), cat('ratio'),
			       ExpoNeRatio, (F or d)['popt3'], (F or d)['rmse2'], p, note,
			       corner='upper right',
			       alt_popt=d['popt3'] if fit_scope == 'both' else None,
			       accent=colors[p])
		fig2.supxlabel(str_dXmax, fontsize=_fs())
		fig2.supylabel(str_Ne_ratio, fontsize=_fs())
		fig2.tight_layout(rect=[0, 0, 1, fig2._viz_top])
		figs['ratio'] = fig2

		# ------------------------ fig3: Edep vs (Ne, Nmu) -- 3-D, same style as 1/2
		# One panel per primary on shared limits: this primary in its own colour,
		# every other primary behind it in grey. The fit is reported numerically
		# in the corner rather than drawn -- a 2-variable linear fit is a plane,
		# and the plane is what fills the cube and hides the points.
		ex3, ey3, ez3 = _pow10(cat('Ne_ground')), _pow10(cat('Nmu_ground')), _pow10(cat('Edep'))
		sx3, sy3, sz3 = 10.0 ** ex3, 10.0 ** ey3, 10.0 ** ez3
		# The plane is no longer drawn, so the pooled parameters have nowhere to
		# live but the subtitle -- the panel notes carry each primary's own.
		sub3 = 'one primary per panel in colour, the others in grey; '
		if glob:
			sub3 += (rf'pooled fit: $\alpha$ = {G["alpha"]:.3e}, '
			         rf'$\beta$ = {G["beta"]:.3e}, '
			         rf'RMSE = {G["rmse3"] / sz3:.3f}')
		else:
			sub3 += ('each panel fitted on its own points; '
			         rf'all primaries pooled would give $\alpha$ = {G["alpha"]:.3e}, '
			         rf'$\beta$ = {G["beta"]:.3e}')
		fig3, axes3 = _facets3d(
			len(particle),
			rf'$E_{{dep}} = \alpha N_e + \beta N_\mu$    {bin_str}', sub3,
			**({'figsize': figsize3d} if figsize3d else {}))

		xlim = np.array([cat('Ne_ground').min(), cat('Ne_ground').max()]) / sx3
		ylim = np.array([cat('Nmu_ground').min(), cat('Nmu_ground').max()]) / sy3
		zlim = np.array([cat('Edep').min(), cat('Edep').max()]) / sz3

		for ax, p in zip(axes3, particle):
			d = D[p]
			x, y = d['Ne_ground'] / sx3, d['Nmu_ground'] / sy3
			z = d['Edep'] / sz3

			# context: the OTHER primaries only. In 2-D the accent points are drawn
			# over the context, but 3-D sorts by depth, so including this primary's
			# own points in the grey layer would let grey land on top of colour.
			other = [q for q in particle if q != p]
			if other:
				ax.scatter(np.concatenate([D[q]['Ne_ground'] for q in other]) / sx3,
				           np.concatenate([D[q]['Nmu_ground'] for q in other]) / sy3,
				           np.concatenate([D[q]['Edep'] for q in other]) / sz3,
				           s=7, color=CONTEXT, lw=0, depthshade=False, rasterized=True)

			ax.scatter(x, y, z, s=13, color=colors[p], lw=0.3,
			           edgecolor=_shade(colors[p], -0.35), depthshade=False,
			           rasterized=True)

			ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_zlim(*zlim)
			kpad = plt.rcParams['font.size'] / BASE_FONT
			ax.set_xlabel(_label(str_Ne_ground, ex3), labelpad=14 * kpad)
			ax.set_ylabel(_label(str_Nmu_ground, ey3), labelpad=14 * kpad)
			ax.set_zlabel(_label('Edep (GeV)', ez3), labelpad=10 * kpad, rotation=90)
			# a 3-D axes' title sits far above the cube; place it in figure space
			ax.text2D(0.02, 0.97, p, transform=ax.transAxes, ha='left', va='top',
			          fontweight='bold', color=INK,
			          fontsize=plt.rcParams['axes.titlesize'])
			# alpha is O(1e-3) here -- fixed-point rounds every primary to 0.001
			# if glob:
			# 	txt = (f'offset = {d["bias3"] / sz3:+.3f}' + '\n'
			# 	       + f'RMSE = {d["grmse3"] / sz3:.3f}   n = {z.size}' + '\n'
			# 	       + rf'(own fit: $\alpha$ = {d["alpha"]:.2e}, $\beta$ = {d["beta"]:.2e})')
			# else:
			# 	txt = (rf'$\alpha$ = {d["alpha"]:.2e}   $\beta$ = {d["beta"]:.2e}'
			# 	       + '\n' + f'RMSE = {d["rmse3"] / sz3:.3f}   n = {z.size}')
			# ax.text2D(0.02, 0.90, txt, transform=ax.transAxes, ha='left', va='top',
			#           fontsize=_fs(8 / 9), color=MUTED, linespacing=1.5)
		figs['edep'] = fig3

		# ------------------------------------------- optional: the 3-D scatter
		if show_3d:
			fig4 = plt.figure(figsize=(8, 8))
			ax4 = fig4.add_subplot(projection='3d')
			ax4.set_box_aspect(None, zoom=0.88)
			ax4.zaxis.set_rotate_label(False)
			shades = ['#86b6ef', '#3987e5', '#256abf', '#104281']  # ordinal by mass
			for p, c in zip(particle, shades):
				d = D[p]
				ax4.scatter(d['Nmu_ground'], d['Ne_ground'], d['Edep'],
				            color=c, s=14, lw=0.3, edgecolor=SURF, label=p,
				            depthshade=False)
			ax4.set_title(bin_str)
			ax4.set_xlabel(str_Nmu_ground)
			ax4.set_ylabel(str_Ne_ground)
			ax4.set_zlabel('Edep (GeV)', labelpad=10, rotation=90)
			ax4.legend(loc='upper left')
			figs['edep3d'] = fig4

	return figs