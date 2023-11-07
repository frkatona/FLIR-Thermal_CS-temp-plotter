import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import time

# def Integrate_Area_Power(y_mid, dy, Q_dx, abs_coeff, height):
#     """
#     Compute the power for a discrete unit area based on its vertical position
#     using the Beer-Lambert law and integrating across the depth.
#     """
#     y_start = y_mid - dy / 2
#     y_end = y_mid + dy / 2
#     power_integral, _ = quad(lambda y: np.exp(-abs_coeff * (height - y)), y_start, y_end)
#     P_area = power_integral * Q_dx
#     return P_area

def Compute_Power_Source_Distribution(X, k, r, dx, height, Q):
    """
    Returns the power distribution 2D array within material
    considering Beer's law decay and laser position.
    """
    Nbeam = height/r
    beamArray = np.linspace(0,r, Nbeam)

    q = (X <= r).astype(float) # 1 where beam is, 0 where it isn't

    # Create a linear space for the depth
    depth = np.linspace(0, height, Nx)
    
    # Apply Beer's law
    intensity = np.exp(-k * height)
    
    # Create a 2D array
    array_2d = np.tile(intensity, (Nx, 1)).T
    
    A = k * depth 
    transmittance = np.exp(-A)
    fraction_absorbed = 1 - transmittance

    # Normalize the entire 2D array
    array_2d /= np.sum(array_2d)
    array_2d *= fraction_absorbed

    return array_2d, transmittance

def MakeLaserArray(height, Nx, k, Q):
    '''create the 2D array of power source nodes with Beer's law decay''' 

    ## create a 1D linear space for the depth ##
    depth_array = np.linspace(0, height, Nx)

    ## apply Beer's law exp distribution to 1D array ##
    intensity = np.exp(-k * depth_array)

    ## tile into a 2D array for all beam nodes ##
    beam_array = np.tile(intensity, (Nx_beam, 1)).T

    ## calculate transmittance ##
    transmittance = intensity[-1]
    fraction_absorbed = 1 - transmittance
    array_power_total = Q / (np.pi * r_beam**2)

    ## normalize the 2D array to beam power ##
    beam_array /= np.sum(beam_array)
    beam_array *= fraction_absorbed
    beam_array *= array_power_total
    
    return beam_array, transmittance

def FillArray(height, Nx, beam_array):
    '''fill out the non-power-source nodes with zeros'''
    full_shape=(Nx, Nx)
    full_array = np.zeros(full_shape)
    full_array[:, :Nx_beam] = beam_array[:, :Nx_beam]
    
    return full_array

    ## power distributed across horizontal elements at surface ##
    circular_crossSection = np.pi * beam_radius_m**2
    n_beam_elements = beam_radius_m / dx
    Q_dx2_0 = (Q / circular_crossSection) / n_beam_elements 
    
    ## power decay ##
    
    for i in range(Nx):
        q_decayed = Q_dx2_0 * np.exp(-abs_coeff * i * dx)
        print(q_decayed)
        q[i, :] *= q_decayed # exp. decay

    # Nx_beam = Nx * (height / beam_radius_m)
    # Q_dx = (Q / circular_crossSection) / Nx_beam  * cheatNumber # * dx # power contained within a length dx of a radial slice of the beam
    
    # for every column, determine Beer's Law power contained in its row elements and apply it to the truthy beam mask
    # for i in range(Ny):
    #     # y_mid = i * dy
    #     # q_dxdy = Integrate_Area_Power(y_mid, dy, Q_dx, abs_coeff, height)
    #     for j in range(Nx):
    #         q[i, j] = (beam_mask[i, j] * q_dxdy)
    
    return  q


def Laplacian_2d(T, dx, dy):
    """
    Compute the Laplacian of the temperature field T using central differences.
    """
    d2Tdx2 = (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
    d2Tdy2 = (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2

    # Reset boundary values to immediate inside neighbors (i.e., adiabatic edges)
    d2Tdx2[0, :] = d2Tdy2[0, :] = 0
    d2Tdx2[-1, :] = d2Tdy2[-1, :] = 0
    d2Tdx2[:, 0] = d2Tdy2[:, 0] = 0
    d2Tdx2[:, -1] = d2Tdy2[:, -1] = 0
    
    return d2Tdx2 + d2Tdy2

def RK4(T, dt, dx, dy, thermal_diffusivity, q):
    """computes T at the next time step with a 4th degree Runge-Kutta"""
    
    def rate(T_current):
        return dt * (thermal_diffusivity * Laplacian_2d(T_current, dx, dy) + q)
    
    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + (k1 + 2*k2 + 2*k3 + k4) / 6

def Boundary_Conditions(T_new, h_conv, T_air, dt, dy):
    '''applies boundary conditions'''
    ## adiabatic edges ##
    T_new[0 , : ] = T_new[1 , : ]
    T_new[ : , 0] = T_new[ : , 1]
    T_new[ : , -1] = T_new[ : , -2]
    ## convective top ##
    T_new[-1, : ] = (T_new[-1, : ] + (h_conv * dy / PDMS_thermal_conductivity_WpmK) * T_air) / (1 + h_conv * dy / PDMS_thermal_conductivity_WpmK)

    return T_new

def Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny), T_0)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    ## track time for progress print ##
    start_time = time.time()
    first_output_time = None

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T_new = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q)
        T = Boundary_Conditions(T_new, h_conv, T_air, dt, dy)
        
        if n in output_indices:
            current_time = time.time()
            output_temperatures.append(T.copy())
            
            if first_output_time is None:
                first_output_time = current_time
            else:
                elapsed_time = current_time - first_output_time
                remaining_indices = len(output_indices) - output_indices.index(n)
                estimated_time_left = (elapsed_time / output_indices.index(n)) * remaining_indices
                print(f"Estimated time to completion: {estimated_time_left:.0f} seconds")
                
            print(f'Computed T at t = {n * dt:.2f} s ({n} / {max(output_indices)})')
           
    return output_temperatures

def Preview_Decay(q, Q, height):
    '''graph of the power distribution from the laser beam power source decay'''
    absorbed_power = np.sum(q)
    transmittance = (Q / circular_crossSection) * 100 # only 1 radial slice of the circular beam
    print(f'dt: {dt:.2e}s (M={M}), Q: {Q}W, absorbed power: {absorbed_power:.1f}W, transmittance: {transmittance:.1f}%')

    plt.figure(figsize=(16, 6))  # Adjusted figure size for side-by-side plots

    # Plot the distribution of q to visualize the power distribution within the material
    plt.subplot(1, 2, 1)  # 1x2 grid, first plot
    plt.imshow(q / 10000, extent=[0, height, 0, height], origin='lower', cmap='hot')
    plt.colorbar(label='Power Density (W/cm^2)')
    plt.title('Power Density Distribution')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')

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

def Plot_T_Slices(output_temperatures, output_times, height, Q, loading, beam_radius_m):  # added parameters Q, loading, and beam_radius_m
    '''plots slices of T at each output time'''
    # Find the maximum temperature value from the final output
    max_temp = np.max(output_temperatures[-1])
    # Compute the highest temperature value that is divisible by 20
    Tmax = int(np.ceil(max_temp / 20.0) * 20)
    Tmin = 20
    ticks = (Tmax - Tmin) // 20 + 1  # Number of discrete colors
    tick_values = np.linspace(Tmin, Tmax, ticks)  # Tick marks at the center of each discrete color

    num_times = len(output_times)
    fig, axes = plt.subplots(1, num_times, figsize=(4*num_times, 4))
    custom_cmap = plt.cm.get_cmap('hot', ticks-1)
    for i, T_out in enumerate(output_temperatures):
        ax = axes[i]
        c = ax.imshow(T_out, extent=[0, height, 0, height], origin='lower', 
                      vmin=Tmin, vmax=Tmax, cmap=custom_cmap)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=tick_values, fraction = .05)
    plt.tight_layout()
    fig.suptitle(f'temperature distribution for Q = {Q} W, loading = {loading:.0e}, and beam radius = {beam_radius_m} m')
    plt.show()


#############################
###       PARAMETERS      ###
#############################

## physical constants ##
height = 1 # height of the simulation space, m
T_0 = 20.0  # Temperature at t = 0, °C
T_air = 20.0  # Temperature of the surrounding air, °C
h_conv = 5 # Convective heat transfer coefficient, W/(m^2 K)
Q = 1000  # Total heat generation rate, W, (i.e., laser power)
loading = 1e-4  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + loading # TC lerps between 0.2 and 0.3 over the loading range 0% to 10%
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

beam_radius_m = r_beam = 0.1 #r_beam is borrowed from testbench implementation

circular_crossSection = np.pi * beam_radius_m**2
Q_2D_top = Q / circular_crossSection # power density at the top of the simulation space, W/m^2
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

output_times = [0, 5, 15, 30]

## test parameters ##
print(Nx, Nx_beam)

#############################
###          MAIN         ###
#############################
if __name__ == '__main__':
    # q = Compute_Power_Source_Distribution(X, abs_coeff, beam_radius_m, dy, height, Q)
    q_beam, transmittance = MakeLaserArray(height, Nx, abs_coeff, Q)
    q = FillArray(height, Nx, q_beam)

    ## plotting ##
    plt.imshow(q, cmap='hot', vmin=0)
    plt.colorbar(label='Normalized Intensity')
    plt.title(f"Beer's Law Decay for k = {abs_coeff} (transmittance {transmittance*100:.1f}%)")
    plt.show()

    Preview_Decay(q, Q, height)
    output_temperatures_RK = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air)
    Plot_T_Slices(output_temperatures_RK, output_times, height, Q, loading, beam_radius_m)