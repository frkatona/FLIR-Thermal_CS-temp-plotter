import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time
from numba import jit

def MakeLaserArrayGreatAgain(height, Nx, a, Q):
    depth_array = np.linspace(0, height, Nx)
    Abs_array = np.exp(-a * depth_array)
    P_array = []
    for i in range(0,len(Abs_array)-1):
        P_array.append(Abs_array[i] - Abs_array[i+1])
    P_array.append(P_array[-1])

    # graph Abs_array
    plt.plot(depth_array, Abs_array)
    plt.xlabel('depth (m)')
    plt.ylabel('absorption')
    plt.title('Absorption vs. Depth')

    # graph P array
    plt.plot(depth_array, P_array)
    plt.xlabel('depth (m)')
    plt.ylabel('power density (W/m^2)')
    plt.title('Power Density vs. Depth')
    plt.show()
    # incorporate Q
    # graph Q array

MakeLaserArrayGreatAgain(1, 50, 1, 100)

def MakeLaserArray(height, Nx, a, Q):
    '''create the 2D array of power source nodes with Beer's law decay''' 
    # what I want at each point is what is lost in absorption at each step, i.e., P(d) = A(d+1) - A(d)


    ## create a 1D linear space for the depth ##
    depth_array = np.linspace(0, height, Nx)

    ## apply Beer's law exp distribution to 1D array ##
    intensity = np.exp(-a * depth_array)
    print(f"max intensity: {intensity[0]}")

    ## tile into a 2D array for all beam nodes ##
    beam_array = np.tile(intensity, (Nx_beam, 1)).T

    ## calculate transmittance ##
    transmittance = intensity[-1]
    fraction_absorbed = 1 - transmittance
    print(f"fraction absorbed: {fraction_absorbed}")
    array_power_total = Q / (np.pi * r_beam**2)
    print(f"array power total: {array_power_total}")

    ## normalize the 2D array to beam power ##
    beam_array /= np.sum(beam_array)
    print(f"beam array sum: {beam_array.sum()}")
    print(f"beam array max: {beam_array.max()}")

    ## scale the normalized array to the extent not transmitted ##
    beam_array *= fraction_absorbed

    ## scale the normalized array to the power expected in this slice ##
    beam_array *= array_power_total
    print(f"beam array sum: {beam_array.sum()}")
    
    return beam_array, transmittance

def FillArray(Nx, beam_array):
    '''fill out the non-power-source nodes with zeros'''
    full_shape=(Nx, Nx)
    full_array = np.zeros(full_shape)
    full_array[:, :Nx_beam] = beam_array[:, :Nx_beam]
    
    return full_array

@jit(nopython=True)
def Laplacian_2d(T, dx, dy):
    """
    Compute the Laplacian of the temperature field T using central differences.
    """
    d2Tdx2 = (np.concatenate((T[1:], T[:1]), axis=0) - 2 * T + np.concatenate((T[-1:], T[:-1]), axis=0)) / dx**2
    d2Tdy2 = (np.concatenate((T[:, 1:], T[:, :1]), axis=1) - 2 * T + np.concatenate((T[:, -1:], T[:, :-1]), axis=1)) / dy**2

    # Reset boundary values to immediate inside neighbors (i.e., adiabatic edges)
    d2Tdx2[0, :] = d2Tdy2[0, :] = 0
    d2Tdx2[-1, :] = d2Tdy2[-1, :] = 0
    d2Tdx2[:, 0] = d2Tdy2[:, 0] = 0
    d2Tdx2[:, -1] = d2Tdy2[:, -1] = 0
    
    return d2Tdx2 + d2Tdy2


@jit(nopython=True)
def RK4(T, dt, dx, dy, thermal_diffusivity, q):
    """computes T at the next time step with a 4th degree Runge-Kutta"""
    
    def rate(T_current):
        return dt * (thermal_diffusivity * Laplacian_2d(T_current, dx, dy) + q)
    
    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

@jit(nopython=True)
def Boundary_Conditions(T_new, h_conv, T_air, dt, dy):
    '''applies boundary conditions'''
    ## adiabatic edges ##
    T_new[0 , : ] = T_new[1 , : ]
    T_new[ : , 0] = T_new[ : , 1]
    T_new[ : , -1] = T_new[ : , -2]
    # T_new[-1, : ] = T_new[-2, : ]

    ## convective top ##
    T_new[-1, : ] = (T_new[-1, : ] + (h_conv * dy / PDMS_thermal_conductivity_WpmK) * T_air) / (1 + h_conv * dy / PDMS_thermal_conductivity_WpmK)

    return T_new

def Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny), T_0)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T_new = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q)
        T = Boundary_Conditions(T_new, h_conv, T_air, dt, dy)
        
        if n in output_indices:
            output_temperatures.append(T.copy())
            
            print(f'Computed T at t = {n * dt:.2f} s ({n} / {max(output_indices)})')
           
    return output_temperatures

def Preview_Decay(q, Q, height):
    '''graph of the power distribution from the laser beam power source decay'''

    absorbed_power = np.sum(q)
    print(f'dt: {dt:.2e}s (M={M}), Q: {Q}W, absorbed power: {absorbed_power:.1f}W, transmittance: {transmittance:.1f}%')

    plt.figure(figsize=(16, 6))  # Adjusted figure size for side-by-side plots

    # Plot the distribution of q to visualize the power distribution within the material
    plt.subplot(1, 2, 1)  # 1x2 grid, first plot
    plt.imshow(q, cmap='hot', vmin=0)
    plt.colorbar(label='power density (W/m^2))')
    plt.title(f"Beer's Law Decay for k = {abs_coeff:0.2f} (transmittance {transmittance*100:.1f}%) for height = {height} m")
    plt.xlabel('simulation nodes')

    # plot a function of exponential decay from the Beer-Lambert law
    plt.subplot(1, 2, 2)  # 1x2 grid, second plot
    plt.plot(np.linspace(0, height, 100), Q * np.exp(-abs_coeff * np.linspace(0, height, 100)), label='power density')
    plt.xlabel('depth (m)')
    plt.ylabel('power density (W/m^2)')
    plt.title('Power Density Decay')
    plt.ylim(0, Q*1.1)
    plt.legend()

    plt.tight_layout()  # Adjusts spacing between plots for better layout
    plt.show()

def Plot_T_Slices(output_temperatures, output_times, height, Q, loading, r_beam, discretize=True):  
    '''plots slices of T at each output time'''
    
    num_times = len(output_times)
    fig, axes = plt.subplots(1, num_times, figsize=(4*num_times, 4))
    
    if discretize:
        max_temp = np.max(output_temperatures[-1])
        Tmax = int(np.ceil(max_temp / 20.0) * 20)
        Tmin = 20
        ticks = (Tmax - Tmin) // 20 + 1  # Number of discrete colors
        tick_values = np.linspace(Tmin, Tmax, ticks)  # Tick marks at the center of each discrete color
        custom_cmap = plt.cm.get_cmap('hot', ticks-1)
    else:
        custom_cmap = plt.cm.get_cmap('hot')
        tick_values = None  # Let Matplotlib select the tick marks
    
    for i, T_out in enumerate(output_temperatures):
        ax = axes[i]
        c = ax.imshow(T_out, extent=[0, height, 0, height], origin='lower', cmap=custom_cmap)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=tick_values, fraction = .05)
    
    plt.tight_layout()
    fig.suptitle(f'temperature distribution for Q = {Q} W, loading = {loading:.0e}, and beam radius = {r_beam} m')
    plt.show()

#############################
###       PARAMETERS      ###
#############################

## physical constants ##
height = 1 # height of the simulation space, m
T_0 = 20.0  # Temperature at t = 0, °C
T_air = 20.0  # Temperature of the surrounding air, °C
h_conv = 5 # Convective heat transfer coefficient, W/(m^2 K)
Q = 100  # Total heat generation rate, W, (i.e., laser power)
loading = 1e-4  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + loading # TC lerps between 0.2 and 0.3 over the loading range 0% to 10%
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

r_beam = 0.1

# circular_crossSection = np.pi * r_beam**2
# Q_2D_top = Q / circular_crossSection # power density at the top of the simulation space, W/m^2
abs_coeff = 0.01 + (loading * 2300) # abs lerps between 0.01 and ~230 over the loading range of 0% to 10%

## simulation parameters ##
Nx = Ny = 50
Nx_beam = int(Nx * (r_beam / height))
dx = dy = height / (Nx - 1)
M = 4
dt = (dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)) / (M/4)  # time step, s; CFL condition for conduction
dt_CFL_conv = (dx**2 / (2 * PDMS_thermal_diffusivity_m2ps * ((h_conv * dx / PDMS_thermal_conductivity_WpmK) + 1)))  # time step, s
if dt_CFL_conv < dt:
    dt = dt_CFL_conv
x = y = np.linspace(0, height, Nx)
X, Y = np.meshgrid(x, y)

output_times = [0, 1, 3, 5]


#############################
###          MAIN         ###
#############################
if __name__ == '__main__':
    q_beam, transmittance = MakeLaserArray(height, Nx, abs_coeff, Q)
    q = FillArray(Nx, q_beam)
    Preview_Decay(q, Q, height)
    output_temperatures_RK = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air)
    Plot_T_Slices(output_temperatures_RK, output_times, height, Q, loading, r_beam, discretize=False)