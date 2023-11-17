import numpy as np
import matplotlib.pyplot as plt
import time
from lmfit import minimize, Parameters, create_params, fit_report
import pandas as pd
from scipy.stats import halfnorm

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

    ### ATTEMPTING HALF-NORMAL DISTRIBUTION ###

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
    # print(f"sum of P_array_norm: {P_array.sum()}")
    # print(f"Q: {Q:.2f} W")
    # print(f"Q_abs: {Q_abs:.2f} W")
    # print(f"P_slice: {P_slice:.2f} W")
    # print(f"slice % of V_cyl: {V_slice / V_cylinder * 100:.2e} %")
    # print(f"node volume: {V_node / cm_convert:.2e} cm^3")
    # print(f"max intensity: {P_array.max() / cm_convert:.4e} W/cm^3")
    # print(f"sum of P_array: {P_array.sum() * V_node:.2f} W")
    
    return P_array, transmittance

def FillArray(Nx, beam_array):
    '''fill out the non-power-source nodes with zeros'''

    full_shape=(Nx, Nx)
    full_array = np.zeros(full_shape)
    full_array[:, :Nx_beam] = beam_array[:, :Nx_beam]
    
    ## plot the 2D array ##
    # plt.imshow(full_array/cm_m_convert, cmap='hot', vmin=0)
    # plt.colorbar(label='power density (W/cm^3))')

    return full_array

def Laplacian_2d(T, dx, dy):
    """Compute the Laplacian of the temperature field T using central differences."""

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
    T_left = np.roll(T_new[-1, :], shift=-1)  # roll to the left
    T_right = np.roll(T_new[-1, :], shift=1)  # roll to the right
    T_below = T_new[-2, :]

    T_new[-1, 1:-1] = (
        (PDMS_thermal_diffusivity_m2ps * dt / (dx**2)) * (
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

def RK4(T, dt, dx, dy, thermal_diffusivity, q, h_conv):
    """computes T at the next time step with a 4th degree Runge-Kutta"""
    # note to self: apply heat outside of RK4 ?
    # accidentally using last T in this function?...consolidate these
    # should apply power separately?
    # are these time steps?  shouldn't they be?  should dt be an argument for rate?

    def rate(T_current):
        '''conduction then boundaries'''
        T_RK = dt * (thermal_diffusivity * Laplacian_2d(T_current, dx, dy) + (q / PDMS_heat_capacity_V_Jpm3K))
        T_RK = Boundary_Conditions(T_RK, h_conv, T_air, dt, dx, PDMS_thermal_diffusivity_m2ps, PDMS_thermal_conductivity_WpmK)
        return T_RK

    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + ((k1 + 2*k2 + 2*k3 + k4) / 6)

def Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air, time_index_map):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    T = np.full((Nx, Ny), T_0)

    output_temperatures = []
    side_temperatures = {}
    top_temperatures = {}

    for n in range(max(time_index_map.values()) + 1):
        T = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv)

        # Check if the current index matches one of the output times
        for time, index in time_index_map.items():
            if n == index:
                side_temperatures[time] = T[:, 0]
                top_temperatures[time] = T[0, :]
                break 

    return output_temperatures, side_temperatures, top_temperatures


def Preview_Decay(q, Q, height, dt, transmittance, M, abs_coeff):
    '''graph of the power distribution from the laser beam power source decay'''

    absorbed_power = np.sum(q)
    # print(f'dt: {dt:.2e}s (M={M}), Q: {Q}W, absorbed power: {absorbed_power:.1f}W, transmittance: {transmittance:.1f}%')

    plt.figure(figsize=(16, 6))  # Adjusted figure size for side-by-side plots

    # plot a function of exponential decay from the Beer-Lambert law
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

def main(h_conv, conductivity_modifier_inner, conductivity_modifier_outer, abs_modifier_inner, abs_modifier_outer, power_offset, r_beam):
    #############################
    ###       PARAMETERS      ###
    #############################

    ## physical constants ##
    global T_air
    T_air = 20.0  # Temperature of the surrounding air, °C
    global height
    height = 0.05 # height of the simulation space, m

    ## physical variables ##
    T_0 = 20.0  # Temperature at t = 0, °C
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
    Nx = Ny = 50
    global Nx_beam
    Nx_beam = int(Nx * (r_beam / height))
    dx = dy = height / (Nx - 1)
    M = 4e2
    dt = (dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)) / (M/4)  # time step, s; CFL condition for conduction
    dt_CFL_convection = (dx**2 / (2 * PDMS_thermal_diffusivity_m2ps * ((h_conv * dx / PDMS_thermal_conductivity_WpmK) + 1)))  # time step, s
    if dt_CFL_convection < dt:
        dt = dt_CFL_convection
    x = y = np.linspace(0, height, Nx)
    X, Y = np.meshgrid(x, y)

    global output_times
    output_times = [5, 15, 20, 30, 60]

    #############################
    ###       SIM MAIN        ###
    #############################

    # t_0 = time.time()
    q_beam, transmittance = MakeLaserArrayGreatAgain(height, Nx, Nx_beam, abs_coeff, Q, r_beam, power_offset)
    q = FillArray(Nx, q_beam)
    q = np.flip(q, axis=0) # flip around y to match the T array
    # Preview_Decay(q, Q, height, dt, transmittance, M, abs_coeff)
    # t_elapsed = time.time() - t_0
    # print(f'elapsed time: {t_elapsed:.2f} s for {len(output_times)} time steps with {Nx}x{Ny} nodes ({t_elapsed / output_times[-1]:.2f} s_irl/s_sim)')
    # Plot_T_Slices(output_temperatures_RK, output_times, height, Q, loading, r_beam, discretize=False)
    
    time_index_map = create_time_index_map(output_times, dt)
    output_temperatures, side_temperatures, top_temperatures = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air, time_index_map)

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
    conductivity_modifier_inner = params['conductivity_modifier_inner'].value
    conductivity_modifier_outer = params['conductivity_modifier_outer'].value
    abs_modifier_inner = params['abs_modifier_inner'].value
    abs_modifier_outer = params['abs_modifier_outer'].value
    power_offset = params['power_offset'].value
    r_beam = params['r_beam'].value

    print(params)

    # Calculate residuals
    side_temperatures, top_temperatures, time_index_map  = main(h_conv, conductivity_modifier_inner, conductivity_modifier_outer, abs_modifier_inner, abs_modifier_outer, power_offset, r_beam)
    
    residuals = []
    for time in output_times:
        index = time_index_map.get(time)
        top_resid = np.subtract(top_temperatures[time], interpolated_exp_data[time]['top'])
        # 'side': np.subtract(side_temperatures[time], interpolated_exp_data[time]['side'])

    residuals.extend(top_resid.ravel())


    print(f"{loop_n}) Residuals, mean: {np.mean(residuals):.8e}, std: {np.std(residuals):.2e}")
    print()

    return residuals

def create_time_index_map(output_times, dt):
    time_index_map = {}
    for time in output_times:
        index = round(time / dt)
        time_index_map[time] = index
    print(f"time index map: {time_index_map}, dt: {dt}")
    return time_index_map

#############################
###       FIT MAIN        ###
#############################

global Nx, Ny, output_times
Nx = Ny = 50
output_times = [5, 15, 20, 30, 60]

# import experimental data
exp_data = pd.read_csv(r'exports\CSVs\lmfit_consolidated\1e-6_70W_temperature_profile.csv')

simulation_x_positions = np.linspace(0, 5, Nx)
simulation_y_positions = np.linspace(0, 5, Ny)

interpolated_exp_data = {}
for time in output_times:
    interpolated_exp_data[time] = {
        'top': interpolate_exp_data(exp_data, simulation_x_positions, time),
        # 'side': interpolate_exp_data(exp_data, simulation_y_positions, time)
    }

global loop_n
loop_n = 0

params = Parameters()
params = create_params(h_conv = {'value':15, 'min':1, 'max':30}, 
                       conductivity_modifier_inner = {'value':20, 'min':1, 'max':100}, 
                       conductivity_modifier_outer = {'value':15, 'min':1, 'max':100}, 
                       abs_modifier_inner = {'value':1e7, 'min':1e5, 'max':1e9}, 
                       abs_modifier_outer = {'value':10, 'min':1, 'max':1e3}, 
                       power_offset = {'value':1, 'min':0.5, 'max':10},
                       r_beam = {'value':0.01, 'min':0.005, 'max':0.015})
result = minimize(objective, params)

print(fit_report(result))