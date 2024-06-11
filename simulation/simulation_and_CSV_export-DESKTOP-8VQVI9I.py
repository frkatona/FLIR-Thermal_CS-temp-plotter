import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import halfnorm
import time
import pandas as pd
from pathlib import Path

'''
This script is for just running the simulation once, showing the associated graphs, and exporting the data to CSVs
'''

def GeneratePowerMask(height, Nx, Nx_beam, atten, Q, r_beam, power_offset):
    '''
    construct array of power density values for the irradiated component of the simulation space
    (the remaining portion of the simulation space with 0 power density is filled in later)
    '''

    ## create normalized array ##
    dx = height / (Nx - 1)
    depth_array = np.linspace(0, height, Nx)
    abs_array_norm = np.exp(-atten * depth_array)
    P_array_norm = np.array([])
    
    for i in range(len(abs_array_norm) - 1):
        # append the difference between the current and next value to the array
        P_array_norm = np.append(P_array_norm, abs_array_norm[i] - abs_array_norm[i+1])

    P_array_norm = np.append(P_array_norm, P_array_norm[-1])

    # distribute 1D P_array_norm into 2D with tiling and dividing by beam nodes
    P_array_norm = np.tile(P_array_norm, (Nx_beam, 1)).T

    ### HALF-NORMAL ###

    # Generate half-normal coefficients for columns
    x = np.linspace(0, 1, Nx_beam)
    half_normal_coeffs = halfnorm.pdf(x, scale=1/3)
    half_normal_coeffs = half_normal_coeffs / half_normal_coeffs[0]  # Normalize to keep the leftmost value as 1

    # Apply the coefficients to each column and replicate across rows
    half_normal_array = np.tile(half_normal_coeffs, (Nx, 1))

    # Multiply P_array_norm with the half-normal distribution array
    P_array_norm *= half_normal_array

    ### /HALF-NORMAL ###

    # normalize P_array_norm to 1 as a sum of all elements
    P_array_norm /= P_array_norm.sum()

    ## incorporate Q ##
     # find transmittance
    transmittance = np.exp(-atten * height)
     # use T to find the fraction of Q absorbed by the cylinder (i.e., total Q not transmitted)
    Q_abs = (1-transmittance) * Q
     # use the Q in that 3D space to find the volumetric P (i.e., P_vol = Q_abs / volume_cylinder)
    V_cylinder = np.pi * r_beam**2 * height
    P_density_cylinder = Q_abs / V_cylinder
     # use P_vol to get the power contained in the simulation space parallel to the beam path
    V_slice = height * r_beam * dx
    P_slice = P_density_cylinder * V_slice
     # distribute that power by scaling the normalized Beer's Law decay array
    P_array = P_array_norm * (P_slice) / (dx*dx*dx) * power_offset

    ## troubleshooting ##
    V_node = dx*dx*dx
    cm_convert = 1000000
    print(f"sum of P_array_norm: {P_array.sum()}")
    print(f"Q: {Q:.2f} W")
    print(f"Q_abs: {Q_abs:.2f} W")
    print(f"P_slice: {P_slice:.2f} W")
    print(f"slice % of V_cyl: {V_slice / V_cylinder * 100:.2e} %")
    print(f"node volume: {V_node / cm_convert:.2e} cm^3")
    print(f"max intensity: {P_array.max() / cm_convert:.4e} W/cm^3")
    print(f"sum of P_array: {P_array.sum() * V_node:.2f} W")
    
    return P_array, transmittance

def FillPowerMask(Nx, beam_array, Nx_beam):
    '''fill out the non-power-source nodes with zeros'''
    full_shape=(Nx, Nx)
    full_array = np.zeros(full_shape)
    full_array[:, :Nx_beam] = beam_array[:, :Nx_beam]
    
    # multiplier = height / Nx * 100

    # # Get current tick locations and labels
    # xticks, xlabels = plt.xticks()
    # yticks, ylabels = plt.yticks()

    # # Set new tick locations and labels
    # plt.xticks(xticks, xticks * multiplier)
    # plt.yticks(yticks, yticks * multiplier)

    plt.xlabel('d$_{radial}$ /cm')
    plt.ylabel('d$_{height}$ /cm')

    ## plot the 2D array ##
    cm_convert = 1000000
    plt.imshow(full_array/cm_convert, cmap='hot', vmin=0)
    plt.colorbar(label='power density (W/cm$^3$)')  # Modified label with superscript 3
    full_array = np.flip(full_array, axis=0)

    return full_array

def Laplacian_2d(T, dx, dy):
    """
    Compute the Laplacian of the temperature field T using central differences.
    """
    T_up = np.roll(T, shift=-1, axis=0)
    T_down = np.roll(T, shift=1, axis=0)
    T_left = np.roll(T, shift=-1, axis=1)
    T_right = np.roll(T, shift=1, axis=1)
    
    d2Tdx2 = (T_up - 2 * T + T_down) / dx**2
    d2Tdy2 = (T_left - 2 * T + T_right) / dy**2

    return d2Tdx2 + d2Tdy2

def Boundary_Conditions(T_new, h_conv, T_air, dt, dx, PDMS_thermal_diffusivity_m2ps, PDMS_thermal_conductivity_WpmK):
    '''
    applies boundary conditions of a convective top and adiabatic edges
    (increasing indices are down in y and to the right in x)
    '''
    ## convective top ##
    T_left = np.roll(T_new[-1, :], shift=-1)
    T_right = np.roll(T_new[-1, :], shift=1)
    T_below = T_new[-2, :]

    T_new[-1, 1:-1] = (
        (PDMS_thermal_diffusivity_m2ps * dt / (dx**2)) * 
        (
            (2 * (h_conv * dx / PDMS_thermal_conductivity_WpmK) * T_air) +
            (T_left[1:-1] + T_right[1:-1] + (2 * T_below[1:-1])) +
            (T_new[-1, 1:-1] * ((dx**2) / (PDMS_thermal_diffusivity_m2ps * dt) - 2 * (h_conv * dx / PDMS_thermal_conductivity_WpmK) - 4))
        )
    )

    ## adiabatic edges ##
    T_new[0 , : ] = T_new[1 , : ] # bottom
    T_new[ : , 0] = T_new[ : , 1] # left
    T_new[ : , -1] = T_new[ : , -2] # right
    # T_new[-1, : ] = T_new[-2, : ] # top
    
    return T_new

def RK4(T, dt, dx, dy, thermal_diffusivity, q, PDMS_heat_capacity_V_Jpm3K, h_conv, T_air, PDMS_thermal_diffusivity_m2ps, PDMS_thermal_conductivity_WpmK):
    """computes T at the next time step with a 4th degree Runge-Kutta"""
    # note to self: should I try applying heat outside of RK4 ?

    def rate(T_current):
        T_RK = dt * (thermal_diffusivity * Laplacian_2d(T_current, dx, dy) + (q / PDMS_heat_capacity_V_Jpm3K))
        T_RK = Boundary_Conditions(T_RK, h_conv, T_air, dt, dx, PDMS_thermal_diffusivity_m2ps, PDMS_thermal_conductivity_WpmK) # accidentally using last T in this function...consolidate these

        return T_RK

    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + ((k1 + 2*k2 + 2*k3 + k4) / 6)

def Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air, PDMS_thermal_conductivity_WpmK, PDMS_heat_capacity_V_Jpm3K):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny), T_0)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, PDMS_heat_capacity_V_Jpm3K, h_conv, T_air, PDMS_thermal_diffusivity_m2ps, PDMS_thermal_conductivity_WpmK)

        if n in output_indices:
            output_temperatures.append(T.copy())
            print(f'Computed T at t = {n * dt:.2f} s ({n} / {max(output_indices)})')
    return output_temperatures

def Preview_Decay(q, Q, height, transmittance, abs_coeff, dt, M):
    '''graph of the power distribution from the laser beam power source decay'''

    absorbed_power = np.sum(q)
    print(f'dt: {dt:.2e}s (M={M}), Q: {Q}W, absorbed power: {absorbed_power:.1f}W, transmittance: {transmittance:.1f}%')

    plt.figure(figsize=(16, 6))  # Adjusted figure size for side-by-side plots

    # plot a function of exponential decay from the Beer-Lambert law
    plt.plot(np.linspace(0, height, 100), Q * np.exp(-abs_coeff * np.linspace(0, height, 100)), label='power density')
    plt.xlabel('depth (m)')
    plt.ylabel('power density (W/m^2)')
    # plt.title('Power Density Decay')
    plt.ylim(0, Q*1.1)
    plt.legend()

    plt.tight_layout()  # Adjusts spacing between plots for better layout
    plt.show()

def Plot_T_Slices(output_temperatures, output_times, height, Q, loading, r_beam, dx, FLIR_measure_depth):  
    '''plots slices of T at each output time and the average temperature of the top 5 rows'''
    
    num_times = len(output_times)
    fig, axes = plt.subplots(2, num_times, figsize=(4*num_times, 8))  # Adjust layout to have two rows of plots

    total_width_cm = 5  # Total width in centimeters
    num_points = output_temperatures[0].shape[1]  # Number of points across the width
    point_spacing = total_width_cm / num_points  # Distance between each point

    custom_cmap = plt.cm.get_cmap('hot')
    tick_values = None
    
    measure_depth_index = int(FLIR_measure_depth / dx) * -1
    print(f'measure_depth_index: {measure_depth_index}')

    max_temp_positions = {}

    for i, T_out in enumerate(output_temperatures):
        # Temperature distribution plot
        # ax = axes[0, i]  # First row for the temperature distribution plot
        plt.figure(figsize=(4, 4))  # New figure for each temperature distribution plot
        ax = plt.gca()  # Get the current axes
        c = ax.imshow(T_out, extent=[0, height, 0, height], origin='lower', cmap=custom_cmap)

        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=tick_values, fraction = 0.05)

        # Add dotted horizontal line at the highest temperature
        max_temp = np.max(T_out)
        max_temp_index = np.unravel_index(np.argmax(T_out), T_out.shape)
        max_temp_positions[output_times[i]] = max_temp_index
        ax.axhline(max_temp_index[0] * dx, color='white', linestyle='dotted', label=f'max T = {max_temp:.2f} °C')

        # Average temperature across some depth from the top
        # ax2 = axes[1, i]
        plt.figure(figsize=(4, 4))  # New figure for each average temperature plot
        ax2 = plt.gca()  # Get the current axes
    
        avg_depth_temp = np.mean(T_out[measure_depth_index:, :], axis=0)
        ax2.plot(avg_depth_temp, label=f'avg T across {FLIR_measure_depth*1000} mm (top {measure_depth_index * -1} rows)')
        # ax2.set_title(f'T_ave across top {FLIR_measure_depth*100} cm at t = {output_times[i]} s')
        ax2.set_xlabel('depth index')
        # ax2.set_ylabel('Average Temperature (°C)')
        ax2.legend()

        distances = [i * point_spacing for i in range(num_points)]
        if i == 0:
            avg_temp_data = pd.DataFrame(avg_depth_temp, index=distances, columns=[f'Time_{output_times[i]}_s'])
        else:
            avg_temp_data[f'Time_{output_times[i]}_s'] = avg_depth_temp

    plt.tight_layout()
    # fig.suptitle(f'Temperature Distribution and Average of Top {FLIR_measure_depth * 100} cm for Q = {Q} W, Loading = {loading:.0e}, and Beam Radius = {r_beam} m')

    # Include the distances as the first column of the DataFrame
    avg_temp_data['Distance_cm'] = distances
    avg_temp_data.set_index('Distance_cm', inplace=True)

    # Save the DataFrame to a CSV file with distances included
    export_path = Path(f'exports/CSVs/simulated_toprow/{Q}W_{loading}_top-row.csv')
    avg_temp_data.to_csv(export_path, mode='w')

    # Scatterplot of max temperature positions over time
    max_temp_times = list(max_temp_positions.keys())
    max_temp_indices = [pos[0] for pos in max_temp_positions.values()]
    # depth_max_temp_cm = (height * 100) - (max_temp_indices / Nx * height * 100)
    export_path = Path(f'exports/CSVs/max_temp_positions/{Q}W_{loading}_max-temp-positions_{Nx}.csv')
    max_temp_data = pd.DataFrame({'Time_s': max_temp_times, 'Position': max_temp_indices})
    max_temp_data.to_csv(export_path, mode='w')

    plt.figure(figsize=(8, 6))
    plt.scatter(max_temp_times, max_temp_indices)
    plt.xlabel('t /s')
    plt.ylabel('Tₘₐₓ depth /cm')
    # plt.title('Max Temperature Positions over Time')
    plt.legend()

    plt.show()

#############################
###       PARAMETERS      ###
#############################

h_conv = 5 #5
conductivity_modifier_inner = 100 #20
conductivity_modifier_outer =  10 #10
abs_modifier_inner = 1e7 #1e7
abs_modifier_outer =  10 #10
power_offset = 1.3 #1.84

## physical constants ##
height = 0.05 # height of the simulation space, m
T_0 = T_air = 25 # Temperature at t = 0, °C
Q = 70  # Total heat generation rate, W, (i.e., laser power)
loading = 1e-5 # mass fraction of CB in PDMS, g/g
r_beam = 0.0125 #0.0125
FLIR_measure_depth = 0.005 # FLIR thermal camera image depth, m

PDMS_thermal_conductivity_WpmK = conductivity_modifier_outer * (0.2 + (loading * conductivity_modifier_inner)) # TC theoretically should lerp between 0.2 and 0.3 over the loading range 0% to 10%
PDMS_density_gpmL = 1.02
PDMS_density_gpm3 = PDMS_density_gpmL * 1e6
PDMS_heat_capacity_m_JpgK = 1.67
PDMS_heat_capacity_V_Jpm3K = PDMS_heat_capacity_m_JpgK * PDMS_density_gpm3
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_heat_capacity_V_Jpm3K)
abs_coeff = abs_modifier_outer * (0.01 + (loading * abs_modifier_inner)) # abs theoretically should lerp between 0.01 and ~500 over the loading range of 0% to 10%

## simulation parameters ##
Nx = Ny = 200
Nx_beam = int(Nx * (r_beam / height))
dx = dy = height / (Nx - 1)
M = 4e1
dt = (dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)) / (M/4)  # time step, s; CFL condition for conduction
dt_CFL_convection = (dx**2 / (2 * PDMS_thermal_diffusivity_m2ps * ((h_conv * dx / PDMS_thermal_conductivity_WpmK) + 1)))  # time step, s
if dt_CFL_convection < dt:
    dt = dt_CFL_convection

# output_times = [0, 1, 2, 5, 10, 20, 40, 70, 90, 120]
output_times = [0, 5, 15, 20, 30, 60]

#############################
###          MAIN         ###
#############################

t_0 = time.time()
q_beam, transmittance = GeneratePowerMask(height, Nx, Nx_beam, abs_coeff, Q, r_beam, power_offset)
q = FillPowerMask(Nx, q_beam, Nx_beam)
Preview_Decay(q, Q, height, transmittance, abs_coeff, dt, M)
output_temperatures_RK = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air, PDMS_thermal_conductivity_WpmK, PDMS_heat_capacity_V_Jpm3K)
t_elapsed = time.time() - t_0
print(f'elapsed time: {t_elapsed:.2f} s for {len(output_times)} time steps with {Nx}x{Ny} nodes ({t_elapsed / output_times[-1]:.2f} s_irl/s_sim)')
Plot_T_Slices(output_temperatures_RK, output_times, height, Q, loading, r_beam, dx, FLIR_measure_depth)