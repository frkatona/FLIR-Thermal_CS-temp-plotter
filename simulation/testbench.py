import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

#############################
###  UTILITY FUNCTIONS    ###
#############################
def compute_power_density_per_voxel(y_mid, dy, Q_dx, abs_coeff, height):
    y_start = y_mid - dy / 2
    y_end = y_mid + dy / 2
    power_integral, _ = quad(lambda y: Q_dx * np.exp(-abs_coeff * (height - y)), y_start, y_end)
    return power_integral / dy

def compute_power_distribution_2d(X, abs_coeff, beam_radius_m, dx, height, Q):
    beam_mask = X <= beam_radius_m
    q = np.zeros((Ny, Nx))
    Q_dx = Q / circular_crossSection / dx
    
    for i in range(Ny):
        y_mid = i * dy
        q_dxdy = compute_power_density_per_voxel(y_mid, dy, Q_dx, abs_coeff, height)
        for j in range(Nx):
            q[i, j] = beam_mask[i, j] * q_dxdy
    
    return q

def compute_laplacian_2d(T, dx, dy):
    d2Tdx2 = (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
    d2Tdy2 = (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2
    d2Tdx2[0, :] = d2Tdy2[0, :] = 0
    d2Tdx2[-1, :] = d2Tdy2[-1, :] = 0
    d2Tdx2[:, 0] = d2Tdy2[:, 0] = 0
    d2Tdx2[:, -1] = d2Tdy2[:, -1] = 0
    return d2Tdx2 + d2Tdy2

def advance_T_RK4(T, dt, dx, dy, thermal_diffusivity, q):
    def rate(T_current):
        return dt * (thermal_diffusivity * compute_laplacian_2d(T_current, dx, dy) + q)
    
    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

def apply_boundary_conditions(T_new, h_conv, T_air, dt, dy):
    T_new[0 , : ] = T_new[1 , : ]
    T_new[ : , 0] = T_new[ : , 1]
    T_new[ : , -1] = T_new[ : , -2]
    T_new[-1, : ] = T_new[-1, : ] - h_conv * (T_new[-1, : ] - T_air) * dt / dy 
    return T_new

#############################
###   PLOTTING FUNCTIONS  ###
#############################
def preview_power_distribution(q, Q, height):
    absorbed_power = np.sum(q)
    transmittance = absorbed_power / (Q / circular_crossSection) * 100 

    plt.figure(figsize=(8, 6))
    plt.imshow(q / 10000, extent=[0, height, 0, height], origin='lower', cmap='hot')
    plt.colorbar(label='Power Density (W/cm^2)')
    plt.title('Power Density Distribution')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()
    
    return transmittance, absorbed_power

def plot_temperature_slices(output_temperatures, output_times, height, Q, loading, beam_radius_m):
    ticks, Tmin, Tmax = 10, 20, 200
    num_times = len(output_times)
    fig, axes = plt.subplots(1, num_times, figsize=(4*num_times, 4))
    custom_cmap = plt.cm.get_cmap('hot', ticks-1)

    for i, T_out in enumerate(output_temperatures):
        ax = axes[i]
        c = ax.imshow(T_out, extent=[0, height, 0, height], origin='lower', vmin=Tmin, vmax=Tmax, cmap=custom_cmap)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=np.arange(Tmin, Tmax, (Tmax-Tmin)//(ticks)), fraction=.05)

    plt.tight_layout()
    fig.suptitle(f'temperature distribution for Q = {Q} W, loading ={loading} s\n beam radius = {beam_radius_m} m')
    plt.subplots_adjust(top=0.88)
    plt.show()

#############################
###    SIMULATION SETUP   ###
#############################
Nx, Ny = 100, 100
height = 0.03
dt = 0.01
output_times = [0.01, 0.5, 2, 4]
thermal_diffusivity = 1e-7
h_conv = 10
T_air = 20
T_initial = 20
beam_radius_m = 0.008
Q = 100
loading = 4
circular_crossSection = np.pi * beam_radius_m**2
X = np.linspace(0, height, Nx).reshape(1, -1)
Y = np.linspace(0, height, Ny).reshape(-1, 1)
abs_coeff = 100

dx, dy = height / Nx, height / Ny

q = compute_power_distribution_2d(X, abs_coeff, beam_radius_m, dx, height, Q)

T = np.full((Ny, Nx), T_initial)
output_temperatures = [T.copy()]

#############################
###      SIMULATION       ###
#############################
for i in range(int(loading / dt)):
    T_new = advance_T_RK4(T, dt, dx, dy, thermal_diffusivity, q)
    T_new = apply_boundary_conditions(T_new, h_conv, T_air, dt, dy)
    T = T_new.copy()

    if (i+1) * dt in output_times:
        output_temperatures.append(T.copy())

#############################
###      PLOTTING         ###
#############################
transmittance, absorbed_power = preview_power_distribution(q, Q, height)
plot_temperature_slices(output_temperatures, output_times, height, Q, loading, beam_radius_m)