import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

## material physical parameters ##
cube_length_m = 10  # Length of the cube, m
T_edge = 25.0  # Temperature at the edges, °C
T_air = 25.0  # Temperature of the surrounding air, °C
h_conv = 10.0  # Convective heat transfer coefficient, W/(m^2 K)
Q = 50  # Total heat generation rate, W
PDMS_thermal_conductivity_WpmK = 0.2e-18 # 0.2 W/mK, artificially low for troubleshooting
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

## mesh-grid parameters ##
Nx, Ny, Nz = 50, 50, 50  # Number of grid points in each dimension
dx = cube_length_m / (Nx - 1)  # Grid spacing in the x direction, m
dy = cube_length_m / (Ny - 1)  # Grid spacing in the y direction, m
dz = cube_length_m / (Nz - 1)  # Grid spacing in the z direction, m
dt = 0.01  # Time step, s
x = np.linspace(0, cube_length_m, Nx)
y = np.linspace(0, cube_length_m, Ny)
z = np.linspace(0, cube_length_m, Nz)
X, Y, Z = np.meshgrid(x, y, z)

## heat source parameters ##
# exponentially decaying power to emulate laser absorption #
beam_radius_m = 1.5e-2

loading = 1e-6 # mass fraction of CB in PDMS, g/g
abs_coeff = 1e+6 # absorption coefficient of CB in PDMS, m^-1
# P = abs_coeff * depth * loading # A = abc, a = abs_coeff, b = depth, c = loading

radius_squared = (X - cube_length_m / 2)**2 + (Y - cube_length_m / 2)**2
beam_mask = radius_squared <= beam_radius_m**2
q = np.zeros((Nx, Ny, Nz)) # initialize
q[beam_mask] = Q / (np.pi * beam_radius_m**2) # radial power density
q *= np.exp(-1 * abs_coeff * loading * (Z - cube_length_m)) # depth-dependent exponential decay

# Extracting the power density values along the central axis of the beam
q_center = q[Nx//2, Ny//2, :]

# show what fraction of power extends beyond the cube
print(f'Power transmitted through entire depth: {(Q - np.sum(q)) / Q:.2%}')

plt.figure(figsize=(8, 6))
plt.plot(z, q_center)
plt.xlabel('Depth (m)')
plt.ylabel('Power Density Fraction of Q')
plt.title('Power Density Distribution along Beam Central Axis')
plt.grid(True)
plt.show()



def Runge_Kutta(T):
    """computes the next time step with a 4th degree Runge-Kutta"""
    k1 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T, -1, axis=2) - 2*T + np.roll(T, 1, axis=2)) / dz**2
                + q)
    
    T_temp = T + k1 / 2
    k2 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=2) - 2*T_temp + np.roll(T_temp, 1, axis=2)) / dz**2
                + q)
    
    T_temp = T + k2 / 2
    k3 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=2) - 2*T_temp + np.roll(T_temp, 1, axis=2)) / dz**2
                + q)
    
    T_temp = T + k3
    k4 = dt * (PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=0) - 2*T_temp + np.roll(T_temp, 1, axis=0)) / dx**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=1) - 2*T_temp + np.roll(T_temp, 1, axis=1)) / dy**2
                + PDMS_thermal_diffusivity_m2ps * (np.roll(T_temp, -1, axis=2) - 2*T_temp + np.roll(T_temp, 1, axis=2)) / dz**2
                + q)
        
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

def Boundary_Conditions(T_new):
    '''applies boundary conditions and surface convection to T_new'''
    ## boundary conditions ##
    T_new[0, :, :] = T_edge
    T_new[-1, :, :] = T_edge
    T_new[:, 0, :] = T_edge
    T_new[:, -1, :] = T_edge
    T_new[:, :, 0] = T_edge
    T_new[:, :, -1] = T_new[:, :, -1] - h_conv * (T_new[:, :, -1] - T_air) * dt / (dz) # convection at top surface
    return T_new

def Compute_T(output_times):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution
    T = np.full((Nx, Ny, Nz), T_edge)

    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    for n in range(max(output_indices) + 1):
        T_new = Runge_Kutta(T)
        T = Boundary_Conditions(T_new)
        
        if n in output_indices:
            output_temperatures.append(T.copy())

        if np.any(T) > 26: # troubleshooting problem where T never changes but the plots show it does
            print(T)
           
    return output_temperatures

def Plot_T_Slices(X, Y, Z, output_temperatures, output_times):
    '''plots slices of T at each output time'''
    fig = plt.figure(figsize=(12, 8))
    custom_cmap = plt.cm.get_cmap('hot', 20) # also, 'coolwarm'
    for i, T_out in enumerate(output_temperatures):
        ax = fig.add_subplot(1, len(output_times), i+1, projection='3d')
        c = ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz//2], 
                        levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        ax.contourf(X[:, :, Nz//2], Y[:, :, Nz//2], T_out[:, :, Nz//2], zdir='z', offset=Z[0, 0, Nz-1], 
                        levels=np.linspace(20, 115, 20), cmap=custom_cmap, alpha=0.75)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('z (m)')
        ax.set_xlim([0, cube_length_m])
        ax.set_ylim([0, cube_length_m])
        ax.set_zlim([0, cube_length_m])
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=np.arange(20, 116, 10))
    plt.tight_layout()
    plt.show()

def save_output_temperatures(output_temperatures, filename="output_temperatures.npz"):
    """
    Save the output temperatures to a file.
    
    Parameters:
    - output_temperatures: List of 3D numpy arrays containing temperature data.
    - filename: Name of the file to save the data to.
    """
    np.savez(filename, *output_temperatures)

output_times = [0, 1, 3, 5]
output_temperatures_RK = Compute_T(output_times)

## plot and/or save to file as npz ##
Plot_T_Slices(X, Y, Z, output_temperatures_RK, output_times)
save_output_temperatures(output_temperatures_RK)