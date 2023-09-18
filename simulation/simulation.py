import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from tqdm import tqdm


## materials system parameters ##
height = 10  # height of the simulation space, m
T_0 = 25.0  # Temperature at t = 0, °C
T_air = 25.0  # Temperature of the surrounding air, °C
h_conv = 10.0  # Convective heat transfer coefficient, W/(m^2 K)
Q = 1  # Total heat generation rate, W, artificially high for troubleshooting
loading = 1e-6  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + loading  # W/mK, should lerp 0.2 to 0.3 as loading increases from 0 to 0.1
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

## heat source parameters ##
beam_radius_m = 0.5
abs_coeff = 0.01 # 0.001 + (loading * 1000)  # absorption coefficient of CB in PDMS, m^-1, should lerp from transparent to opaque with loading

## mesh-grid parameters ##
Nx, Ny = 101, 101  # Number of grid points in each dimension
dx = height / (Nx - 1)  # Grid spacing in the x direction, m
dy = height / (Ny - 1)  # Grid spacing in the z direction, m
dt = 0.01  # Time step, s
x = np.linspace(0, height, Nx)
y = np.linspace(0, height, Ny)
X, Y = np.meshgrid(x, y)

beam_mask = X <= beam_radius_m
q = np.zeros((Nx, Ny))
Q_abs = Q / (np.pi * beam_radius_m**2 * height) # before implementing exponential decay
q[beam_mask] = Q / Q_abs / dx**3

def Preview_Decay(q, dx, dy, Q):
    '''troubleshooting power contained within system'''
    # transmittance equals Q times the integral of exp(-abs_coeff * y) from height to infinity
    transmittance = (Q - np.sum(q * dx**3)) / Q
    print(f'Transmittance through material: {transmittance:.2%}')
    print(f'max q (W): {np.max(q) * dx * dy:.2e}')

    plt.figure(figsize=(8, 6))
    plt.imshow(q, extent=[0, height, 0, height], origin='lower', cmap='hot')
    plt.colorbar(label='Power Density (W/m^2)')
    plt.title('Power Density Distribution')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

def laplacian_2D(T, dx, dy):
    """
    Compute the Laplacian of the temperature field T using central differences.
    """
    d2Tdx2 = (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
    d2Tdy2 = (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2
    return d2Tdx2 + d2Tdy2

def RK4(T, dt, dx, dy, thermal_diffusivity, q):
    """computes T at the next time step with a 4th degree Runge-Kutta"""
    
    def rate(T_current):
        return dt * (thermal_diffusivity * laplacian_2D(T_current, dx, dy) + q)
    
    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

def Boundary_Conditions(T_new):
    '''applies boundary conditions'''
    ## adiabatic edges ##
    T_new[Ny, :] = T_new[Ny-1, :]
    T_new[:, 0] = T_new[:, 1]
    T_new[:, Nx] = T_new[:, Nx-1]
    ## convective top ##
    T_new[0:, -1] = T_new[:, -1] - h_conv * (T_new[:, -1] - T_air) * dt / dy
    return T_new

def Compute_T(output_times):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny), T_0)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T_new = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q)
        T = Boundary_Conditions(T_new)
        
        if n in output_indices:
            output_temperatures.append(T.copy())

        if np.any(T) > 26: # troubleshooting problem where T never changes but the plots show it does
            print(T)
           
    return output_temperatures

def Plot_T_Slices(X, Y, output_temperatures, output_times):
    '''plots slices of T at each output time'''
    num_times = len(output_times)
    fig, axes = plt.subplots(1, num_times, figsize=(4*num_times, 4))
    custom_cmap = plt.cm.get_cmap('hot', 20) # also, 'coolwarm'
    for i, T_out in enumerate(output_temperatures):
        ax = axes[i]
        c = ax.imshow(T_out, extent=[0, height, 0, height], origin='lower', 
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

### MAIN CODE ###
# Preview_Decay(q, dx, dy, Q)

output_times = [0, 1, 3, 5]
output_temperatures_RK = Compute_T(output_times)

Plot_T_Slices(X, Y, output_temperatures_RK, output_times)
save_output_temperatures(output_temperatures_RK, "output_temperatures.npz")