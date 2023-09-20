import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def Integrate_Power_dy(y_mid, dy, Q_dx, abs_coeff, height):
    """
    Compute the volumetric power for a voxel based on its vertical position using the Beer-Lambert law.
    Integrate the power distribution across the voxel height.
    """
    y_start = y_mid - dy / 2
    y_end = y_mid + dy / 2
    power_integral, _ = quad(lambda y: Q_dx * np.exp(-abs_coeff * (height - y)), y_start, y_end)
    return power_integral / dy

def Computer_Power_2d(X, abs_coeff, beam_radius_m, dx, height, Q):
    """
    Sets up the volumetric power distribution using the Beer-Lambert law.
    """
    beam_mask = X <= beam_radius_m
    q = np.zeros((Ny, Nx))
    Q_dx = (Q / circular_crossSection) / dx
    
    # for every y, distribute area power across its x if in beam
    for i in range(Ny):
        y_mid = i * dy
        q_dxdy = Integrate_Power_dy(y_mid, dy, Q_dx, abs_coeff, height)
        for j in range(Nx):
            q[i, j] = (beam_mask[i, j] * q_dxdy) / Nx
    
    return q

def Preview_Decay(q, Q, height):
    '''graph of the power distribution from the laser beam power source decay'''
    absorbed_power = np.sum(q)
    transmittance = absorbed_power / (Q / circular_crossSection) * 100 # only 1 radial slice of the circular beam

    # Plot the distribution of q to visualize the power distribution within the material
    plt.figure(figsize=(8, 6))
    plt.imshow(q / 10000, extent=[0, height, 0, height], origin='lower', cmap='hot')
    plt.colorbar(label='Power Density (W/cm^2)')
    plt.title('Power Density Distribution')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()
    
    return transmittance, absorbed_power

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
    T_new[-1, : ] = T_new[-1, : ] - h_conv * (T_new[-1, : ] - T_air) * dt / dy 
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
           
    return output_temperatures

def Plot_T_Slices(output_temperatures, output_times, height):
    '''plots slices of T at each output time'''
    ticks = 10
    Tmin = 20
    Tmax = 200
    num_times = len(output_times)
    fig, axes = plt.subplots(1, num_times, figsize=(4*num_times, 4))
    custom_cmap = plt.cm.get_cmap('hot', ticks-1) # also, 'coolwarm'
    for i, T_out in enumerate(output_temperatures):
        ax = axes[i]
        c = ax.imshow(T_out, extent=[0, height, 0, height], origin='lower', 
                      vmin=Tmin, vmax=Tmax, cmap=custom_cmap)
        ax.set_title(f't = {output_times[i]} s')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        fig.colorbar(c, ax=ax, label='temperature (°C)', ticks=np.arange(Tmin, Tmax, (Tmax-Tmin)//(ticks)))
    plt.tight_layout()
    # plt.title(f'temperature distribution for Q = {Q} W, loading = {loading*100}%, and beam radius = {beam_radius_m} m')
    plt.show()

#############################
###       PARAMETERS      ###
#############################

## physical constants ##
height = 0.2  # height of the simulation space, m
T_0 = 25.0  # Temperature at t = 0, °C
T_air = 25.0  # Temperature of the surrounding air, °C
h_conv = 10.0  # Convective heat transfer coefficient, W/(m^2 K)
Q = 30 # Total heat generation rate, W, (i.e., laser power)
loading = 1e-4  # mass fraction of CB in PDMS, g/g

PDMS_thermal_conductivity_WpmK = 0.2 + loading # TC lerps between 0.2 and 0.3 over the loading range 0% to 10%
PDMS_density_gpmL = 1.02
PDMS_heat_capacity_JpgK = 1.67
PDMS_thermal_diffusivity_m2ps = PDMS_thermal_conductivity_WpmK / (PDMS_density_gpmL * PDMS_heat_capacity_JpgK)

beam_radius_m = 0.05
circular_crossSection = np.pi * beam_radius_m**2
abs_coeff = 0.01 + (loading * 2300) # abs lerps between 0.01 and ~230 over the loading range of 0% to 10%

## simulation parameters ##
Nx = Ny = 150
dx = dy = height / (Nx - 1)
M = 4
dt = (dx**2 / (PDMS_thermal_diffusivity_m2ps * 4)) / (M/4)  # time step, s 
x = y = np.linspace(0, height, Nx)
X, Y = np.meshgrid(x, y)

output_times = [0, 1, 5, 10]

#############################
###          MAIN         ###
#############################
if __name__ == '__main__':
    q = Computer_Power_2d(X, abs_coeff, beam_radius_m, dy, height, Q)
    transmittance, absorbed_power = Preview_Decay(q, Q, height)
    print(f'dt: {dt:2e}s (M={M}), transmittance: {transmittance:3g}%, Q: {Q}W, absorbed power: {absorbed_power:3g}W')
    output_temperatures_RK = Compute_T(output_times, Nx, Ny, T_0, dt, dx, dy, PDMS_thermal_diffusivity_m2ps, q, h_conv, T_air)
    Plot_T_Slices(output_temperatures_RK, output_times, height)