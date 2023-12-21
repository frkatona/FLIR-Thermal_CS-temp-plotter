import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import halfnorm
from lmfit import minimize, Parameters, create_params, fit_report
import pandas as pd
from pathlib import Path

'''
This script is for running the simulation repeatedly to minimize residuals with experimental data and reading the lmfit report.
'''

# [ ] interpolate simulation temps into experimental positions instead of vice versa
# [x] instead of holding constant, tell lmfit not to vary parameter

def MakeLaserArrayGreatAgain(height, Nx, Nx_beam, atten, Q, r_beam, power_offset):
    '''construct array of power density values for an irradiated cross-section'''

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

    return P_array, transmittance

def FillArray(Nx, beam_array, Nx_beam):
    '''fill out the non-power-source nodes with zeros'''
    full_shape=(Nx, Nx)
    full_array = np.zeros(full_shape)
    full_array[:, :Nx_beam] = beam_array[:, :Nx_beam]
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
    '''applies boundary conditions where increasing indices are down in y and to the right in x'''
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

def Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air, time_index_map, FLIR_measure_depth):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    T = np.full((Nx, Ny), T_0)

    output_temperatures = []
    side_temperatures = {}
    top_temperatures = {}

    for n in range(max(time_index_map.values()) + 1):
        T = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, PDMS_heat_capacity_V_Jpm3K, h_conv, T_air, PDMS_thermal_diffusivity_m2ps, PDMS_thermal_conductivity_WpmK)

        # Check if the current index matches one of the output times
        for time, index in time_index_map.items():
            if n == index:
                measure_depth_index = int(FLIR_measure_depth / dx) * -1
                top_temperatures[time] = np.average(T[0:measure_depth_index, :], axis=0)
                # top_temperatures[time] = T[0, :]
                side_temperatures[time] = T[:, 0]
                
                break 

    return output_temperatures, side_temperatures, top_temperatures

def main(h_conv, Nx, conductivity_modifier_inner, conductivity_modifier_outer, abs_modifier_inner, abs_modifier_outer, power_offset, r_beam, T_air, FLIR_measure_depth):

    #############################
    ###       PARAMETERS      ###
    #############################

    ## physical constants ##
    global height
    height = 0.05 # height of the simulation space, m

    ## physical variables ##
    T_0 = T_air  # Temperature at t = 0, °C
    Q = 70  # Total heat generation rate, W, (i.e., laser power)
    loading = 1e-6 # mass fraction of CB in PDMS, g/g

    ## system properties and conversions ##
    global cm_m_convert
    cm_m_convert = 1e6
    global mL_m3_convert
    mL_m3_convert = 1e6

    global PDMS_thermal_conductivity_WpmK
    PDMS_thermal_conductivity_WpmK = conductivity_modifier_outer * (0.2 + (loading * conductivity_modifier_inner)) # TC theoretically should lerp between 0.2 and 0.3 over the loading range 0% to 10%
    PDMS_density_gpmL = 1.02
    global PDMS_density_gpm3
    PDMS_density_gpm3 = PDMS_density_gpmL * mL_m3_convert
    PDMS_heat_capacity_m_JpgK = 1.67
    global PDMS_heat_capacity_V_Jpm3K
    PDMS_heat_capacity_V_Jpm3K = PDMS_heat_capacity_m_JpgK * PDMS_density_gpm3
    global PDMS_thermal_diffusivity_m2ps
    PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_heat_capacity_V_Jpm3K)
    global abs_coeff
    abs_coeff = abs_modifier_outer * (0.01 + (loading * abs_modifier_inner)) # abs theoretically should lerp between 0.01 and ~500 over the loading range of 0% to 10%

    ## simulation parameters ##
    global Nx_beam
    Nx_beam = int(Nx * (r_beam / height))
    dx = dy = height / (Nx - 1)
    M = 4e1
    dt = (dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)) / (M/4)  # time step, s; CFL condition for conduction
    dt_CFL_convection = (dx**2 / (2 * PDMS_thermal_diffusivity_m2ps * ((h_conv * dx / PDMS_thermal_conductivity_WpmK) + 1)))  # time step, s
    if dt_CFL_convection < dt:
        dt = dt_CFL_convection

    global output_times
    output_times = [0, 5, 15]


    #############################
    ###       SIM MAIN        ###
    #############################

    q_beam, _ = MakeLaserArrayGreatAgain(height, Nx, Nx_beam, abs_coeff, Q, r_beam, power_offset)
    q = FillArray(Nx, q_beam, Nx_beam)
    time_index_map = create_time_index_map(output_times, dt)
    _, side_temperatures, top_temperatures = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air, time_index_map, FLIR_measure_depth)

    return side_temperatures, top_temperatures, time_index_map

def interpolate_exp_data(exp_data, simulation_positions, time):
    '''interpolate to for compatibility with simulation array'''
    interpolated_values = np.interp(
        simulation_positions, 
        exp_data[exp_data['Time_s'] == time]['X_cm' or 'Y_cm'], 
        exp_data[exp_data['Time_s'] == time]['Temperature_C']
    )
    return interpolated_values

def objective(params):
    global loop_n
    loop_n += 1

    # Extract parameters
    h_conv = params['h_conv'].value
    conductivity_modifier_inner = params['conductivity_modifier_inner'].value # REAL1 = 50
    conductivity_modifier_outer = params['conductivity_modifier_outer'].value # REAL1 = 20
    abs_modifier_inner = params['abs_modifier_inner'].value #prev. 1e7
    abs_modifier_outer = params['abs_modifier_outer'].value #prev. 10
    power_offset = params['power_offset'].value # REAL1 = 1.77
    r_beam = params['r_beam'].value # REAL1 = 0.0125
    T_air = params['T_air'].value
    FLIR_measure_depth = params['FLIR_measure_depth'].value

    print(params)

    # Calculate residuals
    side_temperatures, top_temperatures, time_index_map  = main(h_conv, Nx, conductivity_modifier_inner, conductivity_modifier_outer, abs_modifier_inner, abs_modifier_outer, power_offset, r_beam, T_air, FLIR_measure_depth)
    
    residuals = {}
    for time in output_times:
        index = time_index_map.get(time)
        if index is not None and time in top_temperatures:
            residuals[time] = {
                'top': np.subtract(top_temperatures[time], interpolated_exp_data[time]['top']),
                # 'side': np.subtract(side_temperatures[time], interpolated_exp_data[time]['side'])
            }

    combined_residuals = []
    for time in output_times:
        if time in residuals:
            combined_residuals.extend(residuals[time]['top'])
            # combined_residuals.extend(residuals[time]['side'])

    if not combined_residuals:
        print("No residuals found for any of the output times.")
        return np.array([0])  # Return an array with a default value

    # Convert to a numpy array for further processing
    combined_residuals = np.array(combined_residuals)
    flattened_residuals = np.array(combined_residuals)
    if loop_n > 1:
        lastMean = np.mean(flattened_residuals)
    else:
        lastMean = 0
    mean_residual = np.mean(flattened_residuals)
    std_residual = np.std(flattened_residuals)

    print(f"{loop_n}) Residuals - Mean Change: {lastMean/mean_residual:.2e} | Mean: {mean_residual:.2e} | Std: {std_residual:.2e}")
    
    global global_residuals
    global_residuals.append(combined_residuals)

    # print shape of combined_residuals
    print(f"shape of combined_residuals: {combined_residuals.shape}")
    return combined_residuals

def create_time_index_map(output_times, dt):
    time_index_map = {}
    for time in output_times:
        index = round(time / dt)  # Round to the nearest index
        time_index_map[time] = index
    print(f"time index map: {time_index_map}, dt: {dt}")
    return time_index_map


#############################
###       FIT MAIN        ###
#############################

global_residuals = []

global Nx, Ny, output_times
Nx = Ny = 150
output_times = [0, 5, 15]

# import experimental data
exp_data_path = Path('exports/CSVs/lmfit_consolidated/1e-6_70W_temperature_profile.csv')
exp_data = pd.read_csv(exp_data_path)

simulation_x_positions = np.linspace(0, 5, Nx)  # Assuming your grid spans 5 cm
simulation_y_positions = np.linspace(0, 5, Ny)

interpolated_exp_data = {}
for time in output_times:
    interpolated_exp_data[time] = {
        'top': interpolate_exp_data(exp_data, simulation_x_positions, time)
    }

global loop_n
loop_n = 0

params = Parameters()
params.add('h_conv', value=9, min=5, max=20, vary=True)
params.add('conductivity_modifier_inner', value=990, min=100, max=2000, vary=True)
params.add('conductivity_modifier_outer', value=39.3, min=1, max=100, vary=True)
params.add('abs_modifier_inner', value=1.1e7, min=1e7, max=1e8, vary=True)
params.add('abs_modifier_outer', value=9.88, min=5, max=11, vary=True)
params.add('power_offset', value=1, min=0.75, max=1.5, vary=True)
params.add('r_beam', value=0.0100, min=0.009, max=0.011, vary=True)
params.add('T_air', value=24, min=15, max=25, vary=True)
params.add('FLIR_measure_depth', value=1e-2, min=1e-3, max=2e-2, vary=True)

result = minimize(objective, params)

print(fit_report(result))

plt.figure(figsize=(10, 6))
for i, res in enumerate(global_residuals):
    plt.plot(res, label=f'Iteration {i+1}')
plt.xlabel('Data Point Index')
plt.ylabel('Residual')
plt.title('Residuals at Each Iteration')
plt.legend()
plt.show()