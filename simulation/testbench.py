import numpy as np
import matplotlib.pyplot as plt

def MakeLaserArray(depth_array, Nx, beam_radius, k):    
    ## apply Beer's law exp distribution ##
    intensity = np.exp(-k * depth_array)

    ## tile into a 2D array ##
    array_2d = np.tile(intensity, (Nx, 1)).T

    ## calculate transmittance ##
    transmittance = intensity[-1]
    fraction_absorbed = 1 - transmittance

    ## normalize the 2D array ##
    array_2d /= np.sum(array_2d)
    array_2d *= fraction_absorbed
    
    return array_2d, transmittance

def FillArray(height, Nx, k, Nx_beam):
    # Create a 1D linear space for the depth
    depth_array = np.linspace(0, height, Nx)

    power_array, transmittance = MakeLaserArray(depth_array, Nx, beam_radius, k)
    
    # Create the full-sized array initialized to zeros
    full_array = np.zeros(full_shape)
    
    # Copy the values from the small array to the leftmost columns of the large array
    full_array[:, :Nx_beam] = power_array[:, :Nx_beam]
    
    return full_array, transmittance

## input parameters ##
height = 0.4
beam_radius = 0.04
Nx = 1000
Nx_beam = int(Nx * (beam_radius / height))
full_shape=(Nx, Nx)
k = 5

# Apply the function to the larger array with only leftmost 100 columns influenced by Beer's law
resultant_partial_large_array, transmittance = FillArray(height, Nx, k, Nx_beam)

# Plot the resultant 2D array of the larger shape
plt.imshow(resultant_partial_large_array, cmap='hot', vmin=0)
plt.colorbar(label='Normalized Intensity')
plt.title(f"Beer's Law Decay for k = {k} (transmittance {transmittance*100:.1f}%)")
plt.show()

# Print the absorbed power for verification
absorbed_power_partial_large = np.sum(resultant_partial_large_array)
print("absorbed_power:", absorbed_power_partial_large)