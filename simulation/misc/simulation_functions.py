import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def Integrate_Power_dy(y_mid, dy, Q_0, abs_coeff, height):
    """
    Compute the volumetric power for a voxel based on its vertical position using the Beer-Lambert law.
    Integrate the power distribution across the voxel height.
    """
    y_start = y_mid - dy / 2
    y_end = y_mid + dy / 2
    power_integral, _ = quad(lambda y: Q_0 * np.exp(-abs_coeff * (height - y)), y_start, y_end)
    return power_integral / dy

def Volumetric_Power(beam_mask, abs_coeff, beam_radius_m, dy, height, Q):
    """
    Sets up the volumetric power distribution using the Beer-Lambert law.
    """
    Ny, Nx = beam_mask.shape
    q = np.zeros((Ny, Nx))
    
    Q_0 = Q / (np.pi * beam_radius_m**2)  # Average power density over the beam area
    
    for j in range(Nx):
        for i in range(Ny):
            if beam_mask[i, j]:
                y_mid = i * dy
                q[i, j] = Integrate_Power_dy(y_mid, dy, Q_0, abs_coeff, height)
                
    return q

def Preview_Decay(q, Q, height):
    '''troubleshooting power contained within system'''
    absorbed_power = np.sum(q)
    transmittance = (Q - absorbed_power) / Q * 100

    # Plot the distribution of q to visualize the power distribution within the material
    plt.figure(figsize=(8, 6))
    plt.imshow(q, extent=[0, height, 0, height], origin='lower', cmap='hot')
    plt.colorbar(label='Power Density (W/m^2)')
    plt.title('Power Density Distribution')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()
    
    return transmittance, absorbed_power

def laplacian_2D(T, dx, dy):
    """
    Compute the Laplacian of the temperature field T using central differences.
    """
    d2Tdx2 = (np.roll(T, -1, axis=0) - 2*T + np.roll(T, 1, axis=0)) / dx**2
    d2Tdy2 = (np.roll(T, -1, axis=1) - 2*T + np.roll(T, 1, axis=1)) / dy**2

    # Reset boundary values
    d2Tdx2[0, :] = d2Tdy2[0, :] = 0
    d2Tdx2[-1, :] = d2Tdy2[-1, :] = 0
    d2Tdx2[:, 0] = d2Tdy2[:, 0] = 0
    d2Tdx2[:, -1] = d2Tdy2[:, -1] = 0
    
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
    plt.show()
