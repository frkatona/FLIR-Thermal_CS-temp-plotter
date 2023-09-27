import numpy as np
import matplotlib.pyplot as plt

def apply_beers_law_corrected(y_dim, num_elements, k):
    # Create a linear space for the depth
    depth = np.linspace(0, y_dim, num_elements)
    
    # Apply Beer's law
    intensity = np.exp(-k * depth)
    
    # Create a 2D array
    array_2d = np.tile(intensity, (num_elements, 1)).T
    
    A = k * y_dim
    transmittance = np.exp(-A)
    fraction_absorbed = 1 - transmittance

    print("fraction absorbed", fraction_absorbed)

    # Normalize the entire 2D array
    array_2d /= np.sum(array_2d)
    array_2d *= fraction_absorbed
    
    return array_2d, transmittance

# Given values
y_dim = 0.4  # example physical depth
num_elements = 100  # example number of elements
k = 200  # example absorption coefficient

# Apply the corrected function
resultant_array, transmittance = apply_beers_law_corrected(y_dim, num_elements, k)

absorbed_power = np.sum(resultant_array)
print("absorbed_power", absorbed_power)

# Plot the resultant corrected 2D array
plt.imshow(resultant_array, cmap='hot', vmin = 0)
plt.colorbar(label='normalized intensity')
plt.title(f"BL Decay for k = {k} (transmittance {transmittance*100:.1f}%)")
plt.show()