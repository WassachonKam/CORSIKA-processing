import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import integrate

zh_version = 123
co_version = 126
energy = 17 # in lg(eV)
zenith = 80
azimuth = 0.0 
Bx =  0   # µT 
By =   0.0  # µT  
Bz = -50.0 # µT  
primary = 'Proton'
#saveDir = '/work/icecube/users/pgalvezm/sims/CoREAS_vs_ZHAireS/plots/'


class _GroundToShowerCoords:
    """Convert ground-plane (x, y) positions to shower coordinates (vxB, vxvxB)."""

    def __init__(self, x, y, theta_rad, phi_rad, Bx, By, Bz):
        self.x = x
        self.y = y
        v = np.array([
            +np.sin(theta_rad) * np.cos(phi_rad),
            +np.sin(theta_rad) * np.sin(phi_rad),
            -np.cos(theta_rad),
        ])
        B = np.array([Bx, By, Bz])
        B = B / np.linalg.norm(B)
        e1 = np.cross(v, B);  e1 /= np.linalg.norm(e1)
        e2 = np.cross(v, e1); e2 /= np.linalg.norm(e2)
        self._e1 = e1
        self._e2 = e2

    def vxB(self):
        return np.dot(np.array([self.x, self.y, 0.0]), self._e1)

    def vxvxB(self):
        return np.dot(np.array([self.x, self.y, 0.0]), self._e2)


def compute_radiation_energy_2D(t, ant_x, ant_y, Ex, Ey, Ez,
                                 zenith_deg, azimuth_deg, Bx, By, Bz):
    """
    Compute full-band 2D radiation energy (eV) via azimuthal integration in
    shower coordinates, matching the energyfluence() method in functions.py:
        fluence = e0 * c * dt * sum(|Et|²)  [eV/m²]

    The integration follows:
        E_rad = 2π ∫ r · f(r) dr
    where f(r) is the unfiltered energy fluence along the φ=90° arm
    (vxvxB direction) and r is the perpendicular distance to the shower axis.

    Parameters
    ----------
    t : 2D array (n_ant, n_time)
        Time in ns.
    ant_x, ant_y : 1D arrays (n_ant,)
        Antenna positions in meters (ground coordinates).
    Ex, Ey, Ez : 2D arrays (n_ant, n_time)
        Electric field components in µV/m.
    zenith_deg : float
        Shower zenith angle in degrees.
    azimuth_deg : float
        Shower azimuth angle in degrees.
    Bx, By, Bz : float
        Geomagnetic field components in any consistent unit (only direction matters).

    Returns
    -------
    E_rad : float
        Full-band radiation energy in eV.
    fluence : 1D array
        Unfiltered energy fluence per antenna (eV/m²).
    """

    # Physical constants — match functions.py exactly
    e0       = 8.854e-12  # F/m  vacuum permittivity
    c        = 2.99e8     # m/s  speed of light
    eV_per_J = 6.242e18   # eV/J

    theta_rad = np.deg2rad(zenith_deg)
    phi_rad   = np.deg2rad(azimuth_deg)

    n_ant = len(ant_x)

    # dt uniform across antennas; inferred from first antenna's time array
    dt_s = (t[0][1] - t[0][0]) * 1e-9   # ns → s

    # --- Unfiltered energy fluence per antenna (eV/m²) ---
    # µV/m → V/m: |E|² in (µV/m)² needs ×1e-12
    fluence = np.array([
        e0 * c * np.sum((Ex[i]**2 + Ey[i]**2 + Ez[i]**2) * 1e-12) * dt_s * eV_per_J
        for i in range(n_ant)
    ])

    # --- Project onto shower coordinates ---
    vxB_arr   = np.empty(n_ant)
    vxvxB_arr = np.empty(n_ant)
    for i in range(n_ant):
        g2s           = _GroundToShowerCoords(ant_x[i], ant_y[i],
                                               theta_rad, phi_rad, Bx, By, Bz)
        vxB_arr[i]    = g2s.vxB()
        vxvxB_arr[i]  = g2s.vxvxB()

    # --- Select φ = 90° arm (vxvxB ≥ 0, |vxB| ≈ 0) ---
    angle    = np.arctan2(vxvxB_arr, vxB_arr)
    mask     = (np.round(angle, 2) == np.round(np.pi / 2, 2)) & (vxvxB_arr >= 0)
    sort_idx = np.argsort(vxvxB_arr[mask])
    r        = vxvxB_arr[mask][sort_idx]
    f_sel    = fluence[mask][sort_idx]

    # --- 2π ∫ r · f(r) dr ---
    E_rad = 2.0 * np.pi * integrate.trapezoid(r * f_sel, r)   # eV

    return E_rad, fluence

def get_max_E_per_ant_per_band(file, software):
    if not isinstance(software, str):
        raise TypeError("Software must be either 'coreas' or 'zhaires'.")
    software = software.lower()
    if software=='coreas':
        x_conv_factor = 1/100 # cm to m
        t_conv_factor = 1e9 # s to ns
        E_conv_factor = 29979.249999996064 * 1e6 # Conversion factor from statV/cm to uV/m
        software = 'CoREAS'
    elif software=='zhaires':
        x_conv_factor, t_conv_factor = 1, 1 # cm to m
        E_conv_factor = 1e6 # Conversion factor from V to uV/m
        software = 'ZHAireS'
    else:
        raise TypeError("Software must be either 'coreas' or 'zhaires'.")
    
    if not isinstance(file, str) or file[-3:]!='npy':
        raise TypeError("File must be path to npy array.")

    data = np.load(file, allow_pickle = True)
    x = data.item()['x']
    y = data.item()['y']
    t = data.item()['t']
    Ex = data.item()['Ex']
    Ey = data.item()['Ey']
    Ez = data.item()['Ez']

    # Keep all antennas (full 2D layout) for radiation energy calculation
    x_all = x * x_conv_factor
    y_all = y * x_conv_factor
    t_all = t * t_conv_factor
    Ex_all = Ex * E_conv_factor
    Ey_all = Ey * E_conv_factor
    Ez_all = Ez * E_conv_factor

    # Look only along one axis (y == 0) for E-field max plots
    x = x[y == 0] * x_conv_factor
    t = t[y == 0] * t_conv_factor
    Ex = Ex[y == 0] * E_conv_factor
    Ey = Ey[y == 0] * E_conv_factor
    Ez = Ez[y == 0] * E_conv_factor

    # Sort
    indices = np.argsort(x)
    x = x[indices]
    Ex = Ex[indices]
    Ey = Ey[indices]
    Ez = Ez[indices]
    E_tot = (Ex**2 + Ey**2 + Ez**2)**(1/2)

    # Filter time domain E field
    def freq_filter(t, Ex, Ey, Ez, band_l, band_u):
        dt = (t[1] - t[0])*1e-9
        Ex, Ey, Ez, t = np.array(Ex, dtype=float), np.array(Ey, dtype=float), np.array(Ez, dtype=float), np.array(t, dtype=float)
        f = np.fft.rfftfreq(len(t), dt)
        # Ex
        fft_Ex = np.fft.rfft(Ex)
        fft_Ex[np.logical_or(f<band_l, f>band_u)] = 0
        Ex = np.fft.irfft(fft_Ex)
        # Ey
        fft_Ey = np.fft.rfft(Ey)
        fft_Ey[np.logical_or(f<band_l, f>band_u)] = 0
        Ey = np.fft.irfft(fft_Ey)
        # Ez
        fft_Ez = np.fft.rfft(Ez)
        fft_Ez[np.logical_or(f<band_l, f>band_u)] = 0
        Ez = np.fft.irfft(fft_Ez)
        E_tot = (Ex**2 + Ey**2 + Ez**2)**(1/2)
        return E_tot
    full_band_maxes = np.max(E_tot, axis = 1)
    band_1_maxes = np.max(np.abs(np.array([freq_filter(t[i], Ex[i], Ey[i], Ez[i], 30e6, 100e6) for i in range(0, len(x))])), axis = 1)
    band_2_maxes = np.max(np.abs(np.array([freq_filter(t[i], Ex[i], Ey[i], Ez[i], 100e6, 300e6) for i in range(0, len(x))])), axis = 1)
    band_3_maxes = np.max(np.abs(np.array([freq_filter(t[i], Ex[i], Ey[i], Ez[i], 300e6, 500e6) for i in range(0, len(x))])), axis = 1)
    band_4_maxes = np.max(np.abs(np.array([freq_filter(t[i], Ex[i], Ey[i], Ez[i], 60e6, 500e6) for i in range(0, len(x))])), axis = 1)

    # --- 2D radiation energy using all antennas and shower coordinates ---
    # Full-band (unfiltered): direct time-domain sum of |E|², matching
    # energyfluence() in functions.py (e0 * c * bin_width * sum|Et|²).
    E_rad_full, _ = compute_radiation_energy_2D(
        t_all, x_all, y_all, Ex_all, Ey_all, Ez_all,
        zenith_deg=zenith, azimuth_deg=azimuth,
        Bx=Bx, By=By, Bz=Bz
    )
    
    print(E_rad_full)

    return x, full_band_maxes, band_1_maxes, band_2_maxes, band_3_maxes, band_4_maxes, software, E_rad_full

# zh_file = f'/work/icecube/users/pgalvezm/sims/PBRRadioSims/simulations/v{zh_version}/v{zh_version}_all_tasks_combined.npy'
# co_file = f'/work/icecube/users/pgalvezm/sims/corsika/corsika-78000/my_runs/CoREAS_vs_ZHAireS/simulations/v{co_version}/v{co_version}_all_antennas_combined.npy'
# x_zh, zh_full, zh_1, zh_2, zh_3, zh_4, zh, E_rad_zh = get_max_E_per_ant_per_band(zh_file, 'zhaires')
# x_co, co_full, co_1, co_2, co_3, co_4, co, E_rad_co = get_max_E_per_ant_per_band(co_file, 'coreas')


# #-------- Plot CoREAS and ZHAireS with ratio panel
# co_colors = ['#004D40', '#000000', '#FFB000', '#1E88E5', '#EA2A70']
# zh_colors = ['#004D40', '#000000', '#FFB000', '#1E88E5', '#EA2A70']

# # Create figure with two subplots
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5), 
#                                 gridspec_kw={'height_ratios': [1, 0.25]},
#                                 sharex=True)

# # Top panel: Electric field comparison
# # Plot Full Bandwidth
# ax1.plot(x_co, co_full, label='Full Bandwidth (CoREAS)', color=co_colors[0], linewidth=1, alpha=0.5)
# ax1.plot(x_zh, zh_full, label='Full Bandwidth (ZHAireS)', color=zh_colors[0], linestyle='--', linewidth=1, alpha=0.5)

# # Low Frequency Band
# ax1.plot(x_co, co_1, label='30MHz-100MHz (CoREAS)', color=co_colors[1], linewidth=1, alpha=0.5)
# ax1.plot(x_zh, zh_1, label='30MHz-100MHz (ZHAireS)', color=zh_colors[1], linestyle='--', linewidth=1, alpha=0.5)

# # Mid Frequency Band
# ax1.plot(x_co, co_2, label='100MHz-300MHz (CoREAS)', color=co_colors[2], linewidth=1, alpha=0.5)
# ax1.plot(x_zh, zh_2, label='100MHz-300MHz (ZHAireS)', color=zh_colors[2], linestyle='--', linewidth=1, alpha=0.5)

# # High Frequency Band
# ax1.plot(x_co, co_3, label='300MHz-500MHz (CoREAS)', color=co_colors[3], linewidth=1, alpha=0.5)
# ax1.plot(x_zh, zh_3, label='300MHz-500MHz (ZHAireS)', color=zh_colors[3], linestyle='--', linewidth=1, alpha=0.5)

# # PBR Band
# ax1.plot(x_co, co_4, label='PBR Bandwidth (60MHz-500MHz) (CoREAS)', color=co_colors[4], linewidth=1.5)
# ax1.plot(x_zh, zh_4, label='PBR Bandwidth (60MHz-500MHz) (ZHAireS)', color=zh_colors[4], linestyle='--', linewidth=1.5)

# ax1.set_ylabel(r'$E_{\mathrm{max}} (\mu V/m)$')
# ax1.grid(which='both', alpha=0.3)
# ax1.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False)

# # --- Cherenkov ring diameter: distance between the two peaks of the PBR E-field profile ---
# # For each software, find the x positions of the maximum E-field value along the arm.
# # The Cherenkov ring diameter is the distance between the positive and negative x peaks.
# co_pbr_pos = x_co[x_co >= 0]
# co_pbr_neg = x_co[x_co <  0]
# co_4_pos   = co_4[x_co >= 0]
# co_4_neg   = co_4[x_co <  0]
# co_peak_pos = co_pbr_pos[np.argmax(co_4_pos)] if len(co_pbr_pos) > 0 else np.nan
# co_peak_neg = co_pbr_neg[np.argmax(co_4_neg)] if len(co_pbr_neg) > 0 else np.nan
# co_cherenkov_diameter = co_peak_pos - co_peak_neg

# zh_pbr_pos = x_zh[x_zh >= 0]
# zh_pbr_neg = x_zh[x_zh <  0]
# zh_4_pos   = zh_4[x_zh >= 0]
# zh_4_neg   = zh_4[x_zh <  0]
# zh_peak_pos = zh_pbr_pos[np.argmax(zh_4_pos)] if len(zh_pbr_pos) > 0 else np.nan
# zh_peak_neg = zh_pbr_neg[np.argmax(zh_4_neg)] if len(zh_pbr_neg) > 0 else np.nan
# zh_cherenkov_diameter = zh_peak_pos - zh_peak_neg

# ax1.text(
#     1.12, 0.2,
#     f"Radiation Energy (full band) \n"
#     f"CoREAS:  {E_rad_co:.1e} eV\n"
#     f"ZHAireS: {E_rad_zh:.1e} eV\n\n"
#     f"Diameter Cherenkov ring CoREAS:  {co_cherenkov_diameter:.0f} m\n"
#     f"Diameter Cherenkov ring ZHAireS: {zh_cherenkov_diameter:.0f} m",
#     transform=ax1.transAxes,
#     verticalalignment='top',
#     fontsize=9
# )
# ax1.set_title(f'$10^{{{energy}}}$eV,  {zenith}°,  {primary}')

# # Bottom panel: Ratio plot (ZHAireS / CoREAS) for ALL bands
# # Store bands in lists for cleaner looping
# co_bands = [co_full, co_1, co_2, co_3, co_4]
# zh_bands = [zh_full, zh_1, zh_2, zh_3, zh_4]
# band_labels = [
#     'Full Bandwidth',
#     '30MHz-100MHz',
#     '100MHz-300MHz',
#     '300MHz-500MHz',
#     'PBR Bandwidth (60MHz-500MHz)'
# ]

# for i, (co_band, zh_band) in enumerate(zip(co_bands, zh_bands)):

#     # Interpolate ZHAireS onto CoREAS grid
#     zh_interp = interp1d(x_zh, zh_band, kind='linear', fill_value='extrapolate')
#     zh_on_co_x = zh_interp(x_co)

#     # Avoid division by zero
#     ratio = zh_on_co_x/co_band

#     # Use same linewidth logic as top panel
#     lw = 1.5 if i == 4 else 1
#     alpha = 1 if i == 4 else 0.5

#     ax2.plot(
#         x_co,
#         ratio,
#         color=co_colors[i],
#         linewidth=lw,
#         alpha=alpha,
#         #label=f'{band_labels[i]}'
#     )

# # Reference lines
# ax2.axhline(y=1.0, color='black', linestyle='--', linewidth=1, alpha=0.7, label='equal')
# #ax2.axhline(y=0.5, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
# ax2.axhline(y=2.0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)

# ax2.set_xlabel('X (m)')
# ax2.set_ylabel('Ratio (ZHAireS / CoREAS)')
# ax2.grid(which='both', alpha=0.3)
# ax2.legend(loc='upper right', bbox_to_anchor=(1.22, 1), frameon=False)
# ax2.set_yscale('log')
# ax2.set_ylim(bottom=0.1, top=10)

# # Adjust layout to prevent label overlap
# plt.tight_layout()
# #plt.savefig(f'{saveDir}E_per_band_comparison_CoREAS=v{co_version}_vs_ZHAireS=v{zh_version}_with_E2.pdf', bbox_inches='tight')

# # --- Average antenna spacing along the arm ---
# avg_spacing_co = np.mean(np.diff(np.sort(x_co)))
# avg_spacing_zh = np.mean(np.diff(np.sort(x_zh)))

# print(f"Average antenna spacing along arm:")
# print(f"  CoREAS:  {avg_spacing_co:.2f} m")
# print(f"  ZHAireS: {avg_spacing_zh:.2f} m")
# print(len(x_co))
