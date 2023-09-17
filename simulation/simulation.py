import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

## materials system parameters ##
cube_length_m = 10  # Length of the cube, m
T_edge = 25.0  # Temperature at the edges, °C
T_air = 25.0  # Temperature of the surrounding air, °C
h_conv = 10.0  # Convective heat transfer coefficient, W/(m^2 K)
Q = 1e+4  # Total heat generation rate, W, artificially high for troubleshooting
loading = 1e-6  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + (loading * 0.1)  # W/mK
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

## heat source parameters ##
beam_radius_m = 0.01
abs_coeff = 0.1 + (loading * 1000)  # absorption coefficient of CB in PDMS, m^-1


## mesh-grid parameters ##
Nx, Ny = 101, 101  # Number of grid points in each dimension
dx = cube_length_m / (Nx - 1)  # Grid spacing in the x direction, m
dy = cube_length_m / (Ny - 1)  # Grid spacing in the z direction, m
dt = 0.05  # Time step, s
x = np.linspace(0, cube_length_m, Nx)
y = np.linspace(0, cube_length_m, Ny)
X, Y = np.meshgrid(x, y)

radius_squared = (X - cube_length_m / 2)**2 + (Y - cube_length_m / 2)**2
beam_mask = radius_squared <= beam_radius_m**2
q = np.zeros((Nx, Ny))
q[beam_mask] = Q / (np.pi * beam_radius_m**2)  # power density distributed across the entire surface

## troubleshooting power contained within system ## 
transmittance = (Q - np.sum(q * dx * dy)) / Q
print(f'Transmittance through material: {transmittance:.2%}')
print(f'max q (W): {np.max(q) * dx * dy:.2e}')

plt.figure(figsize=(8, 6))
plt.imshow(q, extent=[0, cube_length_m, 0, cube_length_m], origin='lower', cmap='hot')
plt.colorbar(label='Power Density (W/m^2)')
plt.title('Power Density Distribution')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()

def Runge_Kutta(T):
    """computes T at the the next time step with a 4th degree Runge-Kutta"""
    k1 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2
                + q)
    
    T_temp = T + k1 / 2
    k2 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                + q)
    
    T_temp = T + k2 / 2
    k3 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                + q)
    
    T_temp = T + k3
    k4 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                + q)
    
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

def Boundary_Conditions(T_new):
    '''applies boundary conditions and surface convection to T_new'''
    T_new[0, :] = T_edge
    T_new[-1, :] = T_edge
    T_new[:, 0] = T_edge
    T_new[:, -1] = T_new[:, -1] - h_conv * (T_new[:, -1] - T_air) * dt / dy # convection at top surface
    return T_new


def Compute_T(output_times):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny), T_edge)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    # progress bar
    total_iterations = max(output_indices) + 1
    one_percent_iterations = int(total_iterations * 0.01)
    pbar = tqdm(total=total_iterations, desc="Computing", postfix="0% done")

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T_new = Runge_Kutta(T)
        T = Boundary_Conditions(T_new)
        
        if n in output_indices:
            output_temperatures.append(T.copy())

        if np.any(T) > 26: # troubleshooting problem where T never changes but the plots show it does
            print(T)
           
    pbar.close()
    return output_temperatures

def Plot_T_Slices(X, Y, output_temperatures, output_times):
    '''plots slices of T at each output time'''
    num_times = len(output_times)
    fig, axes = plt.subplots(1, num_times, figsize=(4*num_times, 4))
    custom_cmap = plt.cm.get_cmap('hot', 20) # also, 'coolwarm'
    for i, T_out in enumerate(output_temperatures):
        ax = axes[i]
        c = ax.imshow(T_out, extent=[0, cube_length_m, 0, cube_length_m], origin='lower', 
                      vmin=20, vmax=115, cmap=custom_cmap)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=np.arange(20, 116, 10))
    plt.tight_layout()
    plt.show()

def save_output_temperatures(output_temperatures, filename):
    """
    Save the output temperatures to a file.
    
    Parameters:
    - output_temperatures: List of 2D numpy arrays containing temperature data.
    - filename: Name of the file to save the data to.
    """
    np.savez(filename, *output_temperatures)

output_times = [0, 1, 3, 5]
output_temperatures_RK = Compute_T(output_times)

## plot and/or save to file as npz ##
Plot_T_Slices(X, Y, output_temperatures_RK, output_times)
save_output_temperatures(output_temperatures_RK, "output_temperatures.npz")