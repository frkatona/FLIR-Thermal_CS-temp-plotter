import numpy as np
import matplotlib.pyplot as plt
import time

def Boundary_Conditions(T, T_new, h_conv, T_air, dt, dx):
    '''applies boundary conditions where increasing indices are down in y and to the right in x'''
    ## convective top ##
    T_new[-1, : ] = (T_new[-1, : ] + (h_conv * dx / PDMS_thermal_conductivity_WpmK) * T_air) / (1 + h_conv * dx / PDMS_thermal_conductivity_WpmK)
    
    ## adiabatic edges ##
    T_new[0 , : ] = T_new[1 , : ] # bottom
    T_new[ : , 0] = T_new[ : , 1] # left
    T_new[ : , -1] = T_new[ : , -2] # right
    
    return T_new

def RK4(T, dt, dx, dy, thermal_diffusivity, q):
    """computes T at the next time step with a 4th degree Runge-Kutta"""
    def Laplacian_2d(T, dx, dy):
        """
        Compute the Laplacian of the temperature field T using central differences.
        """
        d2Tdx2 = (np.roll(T, -1, axis=0) - 2 * T + np.roll(T, 1, axis=0)) / dx**2
        d2Tdy2 = (np.roll(T, -1, axis=1) - 2 * T + np.roll(T, 1, axis=1)) / dy**2
        return d2Tdx2 + d2Tdy2

    def rate(T_current):
        T_RK = dt * (thermal_diffusivity * Laplacian_2d(T_current, dx, dy) + (q / PDMS_heat_capacity_V_Jpm3K))
        T_RK = Boundary_Conditions(T_current, T_RK, h_conv, T_air, dt, dy)
        return T_RK

    k1 = rate(T)
    k2 = rate(T + 0.5 * k1)
    k3 = rate(T + 0.5 * k2)
    k4 = rate(T + k3)
    
    return T + ((k1 + 2*k2 + 2*k3 + k4) / 6)

def Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air):
    '''computes T at each grid point at each time step and returns data for each value in output_time'''
    # initialize temperature distribution and output structures
    T = np.full((Nx, Ny), T_0)
    output_indices = [int(t / dt) for t in output_times] # times enumerated as indices based on time step
    output_temperatures = []

    # loop across each necessary time step
    for n in range(max(output_indices) + 1):
        T = RK4(T, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q)

        if n in output_indices:
            output_temperatures.append(T.copy())
            print(f'Computed T at t = {n * dt:.2f} s ({n} / {max(output_indices)})')
           
    return output_temperatures

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
loading = 1e-4 # mass fraction of CB in PDMS, g/g
conductivity_modifier_inner = 1
conductivity_modifier_outer = 1
abs_modifier_inner = 5000
abs_modifier_outer = 100
r_beam = 0.01

cm_m_convert = 1e6
mL_m3_convert = 1e6

PDMS_thermal_conductivity_WpmK = conductivity_modifier_outer * (0.2 + (loading * conductivity_modifier_inner)) # TC theoretically should lerp between 0.2 and 0.3 over the loading range 0% to 10%
PDMS_density_gpmL = 1.02
PDMS_density_gpm3 = PDMS_density_gpmL * mL_m3_convert
PDMS_heat_capacity_m_JpgK = 1.67
PDMS_heat_capacity_V_Jpm3K = PDMS_heat_capacity_m_JpgK * PDMS_density_gpm3
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_heat_capacity_V_Jpm3K)
abs_coeff = abs_modifier_outer * (0.01 + (loading * abs_modifier_inner)) # abs theoretically should lerp between 0.01 and ~500 over the loading range of 0% to 10%

## simulation parameters ##
Nx = Ny = 100
Nx_beam = int(Nx * (r_beam / height))
dx = dy = height / (Nx - 1)
M = 4e2
dt = (dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)) / (M/4)  # time step, s; CFL condition for conduction
dt_CFL_convection = (dx**2 / (2 * PDMS_thermal_diffusivity_m2ps * ((h_conv * dx / PDMS_thermal_conductivity_WpmK) + 1)))  # time step, s
if dt_CFL_convection < dt:
    dt = dt_CFL_convection
x = y = np.linspace(0, height, Nx)
X, Y = np.meshgrid(x, y)

output_times = [0, 10, 20]

#############################
###          MAIN         ###
#############################

if __name__ == '__main__':
    t_0 = time.time()
    q = np.zeros((Nx, Ny))
    q[:, :Nx_beam] = 1
    q = np.flip(q, axis=0) # flip around y to match the T array
    output_temperatures_RK = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air)
    t_elapsed = time.time() - t_0
    print(f'elapsed time: {t_elapsed:.2f} s for {len(output_times)} time steps with {Nx}x{Ny} nodes ({t_elapsed / output_times[-1]:.2f} s_irl/s_sim)')
    Plot_T_Slices(output_temperatures_RK, output_times, height, Q, loading, r_beam, discretize=False)
