import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def extract_and_compute_stress(file_path, area):
    data = pd.read_csv(file_path, delimiter='\t', skiprows=4)
    
    load = pd.to_numeric(data['Load '], errors='coerce')
    strain = pd.to_numeric(data['Strain 1 '], errors='coerce')
    
    stress = (load / area) * 0.1  # Convert from N/cm^2 to N/mm^2
    return stress, strain

def smooth_data(data, window_size=5):
    """Apply simple moving average smoothing."""
    return data.rolling(window=window_size).mean()

def compute_youngs_modulus(mask, strain, stress):
    strain_range = strain[mask]
    stress_range = stress[mask]
    slope, intercept = np.polyfit(strain_range, stress_range, 1)
    return slope, intercept

def compute_toughness(strain, stress):
    return np.trapz(stress, strain)

# Given data
thickness = 0.2725
width_inch = 1
width_cm = width_inch * 2.54
area = thickness * width_cm

# Folder path where the files are located
folder_path = r'C:\Users\antho\Desktop\sandra_tensile'

# Get all .txt files from the specified folder
all_files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f)) and f.endswith('.txt')]

plt.figure(figsize=(12, 8))
colors = plt.cm.Dark2(np.linspace(0, 1, len(all_files)))

# Loop through each file and process
for idx, file in enumerate(all_files):
    file_path = os.path.join(folder_path, file)
    sample_name = os.path.splitext(file)[0]
    
    # Extract and compute stress and strain
    stress, strain = extract_and_compute_stress(file_path, area)
    
    # Smooth data
    stress_smoothed = smooth_data(stress)
    strain_smoothed = smooth_data(strain)
    
    # Compute Young's modulus
    mask = (strain_smoothed >= 0.001) & (strain_smoothed <= 0.003)
    slope, _ = compute_youngs_modulus(mask, strain_smoothed, stress_smoothed)
    
    # Compute toughness
    toughness = compute_toughness(strain_smoothed, stress_smoothed)
    
    # Adjust sample_name to include YM and toughness
    sample_name += f" (YM: {slope:.2f}, Toughness: {toughness:.2f})"
    
    # Plotting
    plt.plot(strain_smoothed, stress_smoothed, label=sample_name, color=colors[idx], linewidth=3)
    plt.plot(strain_smoothed[mask], slope * strain_smoothed[mask], '--', color='black', linewidth=3)

plt.xlabel('Strain (mm/mm)')
plt.ylabel('Stress (N/mm^2)')  # Adjusted ylabel to show N/mm^2
plt.title('Stress-Strain Curves')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
