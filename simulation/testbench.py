import numpy as np
import matplotlib.pyplot as plt

def apply_decay(y_dim, num_elements, k):
    # Create a linear space for the depth
    depth = np.linspace(0, y_dim, num_elements)
    
    # Apply Beer's law
    intensity = np.exp(-k * depth)
    
    # Create a 2D array
    array_2d = np.tile(intensity, (num_elements, 1)).T
    
    A = k * y_dim
    transmittance = np.exp(-A)
    fraction_absorbed = 1 - transmittance

    # Normalize the entire 2D array
    array_2d /= np.sum(array_2d)
    array_2d *= fraction_absorbed
    
    return array_2d, transmittance

def partial_array_fill(y_dim, num_elements, k, full_shape=(500, 500), apply_to_columns=100):
    # Generate the 2D array with Beer's law applied to the specified number of rows
    small_array, transmittance = apply_decay(y_dim, num_elements, k)
    
    # Create the full-sized array initialized to zeros
    large_array = np.zeros(full_shape)
    
    # Copy the values from the small array to the leftmost columns of the large array
    large_array[:, :apply_to_columns] = small_array[:, :apply_to_columns]
    
    return large_array, transmittance

# Given values
y_dim = 0.4  # example physical depth
num_elements = 500  # updated number of elements for the larger array
k = 200  # example absorption coefficient

# Apply the function to the larger array with only leftmost 100 columns influenced by Beer's law
resultant_partial_large_array, transmittance = partial_array_fill(y_dim, num_elements, k)

# Plot the resultant 2D array of the larger shape
plt.imshow(resultant_partial_large_array, cmap='hot', vmin=0)
plt.colorbar(label='Normalized Intensity')
plt.title(f"Beer's Law Decay for k = {k} (transmittance {transmittance*100:.1f}%)")
plt.show()

# Print the absorbed power for verification
absorbed_power_partial_large = np.sum(resultant_partial_large_array)
print("absorbed_power:", absorbed_power_partial_large)
