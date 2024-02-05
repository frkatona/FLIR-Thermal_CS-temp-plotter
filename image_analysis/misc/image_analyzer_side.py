import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from matplotlib.colors import ListedColormap
from skimage import feature

# Load data from a folder and sort by extracted time
def load_data_from_folder(folder_path):
    files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    files.sort(key=lambda x: int(x.split('.')[0]))  # Sort by time extracted from filename
    data_list = [pd.read_csv(os.path.join(folder_path, file), sep='\t', header=None).values for file in files]
    times = [int(file.split('.')[0]) for file in files]
    return data_list, times

def extract_temperatures_along_line(data, start, end):
    averaged_data = []
    for x in range(start, end):
        values = []
        for dy in [-1, 0, 1]:
            y = start[0] + dy
            if 0 <= y < data.shape[0]:
                values.append(data[y, x])
        averaged_data.append(np.mean(values))
    return np.array(averaged_data)

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
    
    # Initialize variables to store the line parameters for the 0 s and 60 s files
    line_params_0s = None
    line_params_60s = None
    
    # Find the line parameters from the 60 s file
    for data, time in zip(data_list, times):
        if time == 60:
            max_pos = np.unravel_index(data.argmax(), data.shape)
            half_width = data.shape[1] // 2
            start_x = max_pos[1]
            end_x = start_x + half_width
            if end_x > data.shape[1]:
                start_x -= (end_x - data.shape[1])
                end_x = data.shape[1]
            line_params_60s = ((start_x, max_pos[0]), (end_x, max_pos[0]))
            break
    
    # Plot thermal images side by side
    fig, axs = plt.subplots(1, len(data_list), figsize=(15, 5))
    for ax, data, time in zip(axs, data_list, times):
        im = ax.imshow(data, cmap=cmap_10_stops, vmin=vmin_fixed, vmax=vmax_rounded)
        ax.set_title(f"{time} s")
        ax.axis('off')

        # Draw line where data is taken
        if time == 0 and line_params_60s is not None:
            line_params_0s = line_params_60s
            start, end = line_params_0s
        else:
            max_pos = np.unravel_index(data.argmax(), data.shape)
            half_width = data.shape[1] // 2
            start_x = max_pos[1]
            end_x = start_x + half_width
            if end_x > data.shape[1]:
                start_x -= (end_x - data.shape[1])
                end_x = data.shape[1]
            start = (start_x, max_pos[0])
            end = (end_x, max_pos[0])
        
        ax.plot([start[1], end[1]], [start[0], end[0]], color='white', linewidth=2)

    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Temperature (°C)")
    
    # Plot smoothed lines with original data points
    plt.figure(figsize=(10, 6))
    color_palette = plt.cm.Reds(np.linspace(0, 1, len(data_list)))
    for data, label, color in zip(data_list, times, color_palette):
        if label == 0 and line_params_0s is not None:
            start, end = line_params_0s
            temps = extract_temperatures_along_line(data, start[1], end[1])
        else:
            max_pos = np.unravel_index(data.argmax(), data.shape)
            half_width = data.shape[1] // 2
            start_x = max_pos[1]
            end_x = start_x + half_width
            if end_x > data.shape[1]:
                start_x -= (end_x - data.shape[1])
                end_x = data.shape[1]
            temps = extract_temperatures_along_line(data, start_x, end_x)
            
        x_smooth, y_smooth = smooth_curve(range(len(temps)), temps)
        plt.plot(x_smooth, y_smooth, label=f"{label} s", color=color, linewidth=2)
        plt.scatter(range(len(temps)), temps, color=color, marker='o')
        
    plt.title(f"side-view temperature profile")
    plt.xlabel("Position along the line")
    plt.ylabel("Temperature (°C)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Run the functions
folder_path = r'input_txts\lmfit_txts\txts_side\1e-6_70W'  # Replace with your folder path
data_list, times = load_data_from_folder(folder_path)
plot_data(data_list, times)