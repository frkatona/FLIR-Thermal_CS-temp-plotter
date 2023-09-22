import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from matplotlib.colors import ListedColormap

# Load data from a folder and sort by extracted time
def load_data_from_folder(folder_path):
    files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    files.sort(key=lambda x: int(x.split('.')[0]))  # Sort by time extracted from filename
    data_list = [pd.read_csv(os.path.join(folder_path, file), sep='\t', header=None).values for file in files]
    times = [int(file.split('.')[0]) for file in files]
    return data_list, times

# Extract temperatures along a horizontal line centered at the maximum position
def extract_temperatures_along_line(data, max_pos, span=15):
    start = max(0, max_pos[1] - span // 2)
    end = min(data.shape[1], max_pos[1] + span // 2 + 1)
    return data[max_pos[0], start:end]

# Generate a smooth curve using spline interpolation
def smooth_curve(x, y, points=100):
    x_smooth = np.linspace(min(x), max(x), points)
    spline = make_interp_spline(x, y, k=3)
    y_smooth = spline(x_smooth)
    return x_smooth, y_smooth

# Plot the data
def plot_data(data_list, times):
    # Define color maps and temperature range
    base_cmap = plt.cm.inferno
    new_colors_10 = [base_cmap(i) for i in np.linspace(0, 1, 10)]
    cmap_10_stops = ListedColormap(new_colors_10)
    vmin_fixed = 20
    vmax_rounded = (np.max(data_list[-1]) // 20 + 1) * 20
    
    # Plot thermal images side by side
    fig, axs = plt.subplots(1, len(data_list), figsize=(15, 5))
    for ax, data, time in zip(axs, data_list, times):
        im = ax.imshow(data, cmap=cmap_10_stops, vmin=vmin_fixed, vmax=vmax_rounded)
        ax.set_title(f"{time} s")
        ax.axis('off')
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Temperature (°C)")
    
    # Plot smoothed lines with original data points
    plt.figure(figsize=(10, 6))
    color_palette = plt.cm.Reds(np.linspace(0, 1, len(data_list)))
    for data, label, color in zip(data_list, times, color_palette):
        max_pos = np.unravel_index(data.argmax(), data.shape)
        temps = extract_temperatures_along_line(data, max_pos)
        x_smooth, y_smooth = smooth_curve(range(len(temps)), temps)
        plt.plot(x_smooth, y_smooth, label=f"{label} s", color=color, linewidth=2)
        plt.scatter(range(len(temps)), temps, color=color, marker='o')
    plt.title("Temperatures along the horizontal lines (Smoothed)")
    plt.xlabel("Position along the line")
    plt.ylabel("Temperature (°C)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Run the functions
folder_path = r'C:\Users\antho\Desktop\Desktop_tmp\FLIR_tmp\1e-4\text'  # Replace with your folder path
data_list, times = load_data_from_folder(folder_path)
plot_data(data_list, times)