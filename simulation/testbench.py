import numpy as np
import matplotlib.pyplot as plt

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


######################
## input parameters ##
######################

height = 0.4
r_beam = 0.04
Nx = 1000
k = 5
P_tot = 10

Nx_beam = int(Nx * (r_beam / height))


##########
## main ##
##########

## construct array of simulation nodes with local power source ##
beam_array, transmittance = MakeLaserArray(height, Nx, k, P_tot)
full_array = FillArray(height, Nx, beam_array)

## plotting ##
plt.imshow(full_array, cmap='hot', vmin=0)
plt.colorbar(label='Normalized Intensity')
plt.title(f"Beer's Law Decay for k = {k} (transmittance {transmittance*100:.1f}%)")
plt.show()

## print the absorbed power for verification ##
absorbed_power = np.sum(full_array)
print("absorbed_power:", absorbed_power)

"""
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
    """