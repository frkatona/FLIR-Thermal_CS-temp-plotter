import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time
from numba import jit

def MakeLaserArrayGreatAgain(height, Nx, Nx_beam, atten, Q, r_beam):
    '''power source array generator v2'''

    dx = height / (Nx - 1)
    depth_array = np.linspace(0, height, Nx)
    Abs_array_norm = np.exp(-atten * depth_array)
    P_array_norm = np.array([])
    
    for i in range(0,len(Abs_array_norm)-1):
        # append the difference between the current and next value to
        P_array_norm = np.append(P_array_norm, Abs_array_norm[i] - Abs_array_norm[i+1])

    P_array_norm = np.append(P_array_norm, P_array_norm[-1])

    # distribute 1D P_array_norm into 2D with tiling and dividing by beam nodes
    P_array_norm = np.tile(P_array_norm, (Nx_beam, 1)).T
    P_array_norm /= Nx_beam
    
    # normalize P_array_norm to 1 as a sum of all elements
    P_array_norm /= P_array_norm.sum()

    # incorporate Q

     # find transmittance
    T = np.exp(-atten * height)

     # use T to find the fraction of Q absorbed by the cylinder (i.e., total Q not transmitted)
    Q_abs = (1-T) * Q

     # use the Q in that 3D space to find the volumetric P (i.e., P_vol = Q_abs / volume_cylinder)
    V_cylinder = np.pi * r_beam**2 * height
    P_cylinder = Q_abs / V_cylinder

     # use P_vol to get the power contained in the simulation space parallel to the beam path (i.e., P_slice = P_vol / circular_crossSection))
    V_slice = height * r_beam * dx
    P_slice = P_cylinder * (V_cylinder / V_slice)

    # distribute that power by scaling the normalized Beer's Law decay array
    P_array = P_array_norm * (P_slice)



    V_node = dx*dx*dx
    cm_convert = 1000000
    print(f"sum of P_array_norm: {P_array.sum()}")
    print(f"Q: {Q:.2f} W")
    print(f"Q_abs: {Q_abs:.2f} W")
    print(f"P_slice: {P_slice:.2f} W")
    print(f"slice % of V_cyl: {V_slice / V_cylinder * 100:.2e} %")
    print(f"node volume: {V_node / cm_convert:.2e} cm^3")
    print(f"max intensity: {P_array.max() / cm_convert:.2f} W/cm^3")
    print(f"sum of P_array: {P_array.sum() * V_node:.2f} W")
    
    # plot the 2D array
    plt.imshow(P_array/cm_convert, cmap='hot', vmin=0)
    plt.colorbar(label='power density (W/cm^3))')
    plt.show()

    return P_array, T

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
height = 0.05 # height of the simulation space, m
T_0 = 20.0  # Temperature at t = 0, °C
T_air = 20.0  # Temperature of the surrounding air, °C
h_conv = 5 # Convective heat transfer coefficient, W/(m^2 K)
Q = 100  # Total heat generation rate, W, (i.e., laser power)
loading = 1e-4  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + loading # TC lerps between 0.2 and 0.3 over the loading range 0% to 10%
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

r_beam = 0.01

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

output_times = [0, 1, 3]


#############################
###          MAIN         ###
#############################
if __name__ == '__main__':
    q_beam, transmittance = MakeLaserArrayGreatAgain(height, Nx, Nx_beam, abs_coeff, Q, r_beam)
    q = FillArray(Nx, q_beam)
    Preview_Decay(q, Q, height)
    output_temperatures_RK = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air)
    Plot_T_Slices(output_temperatures_RK, output_times, height, Q, loading, r_beam, discretize=False)