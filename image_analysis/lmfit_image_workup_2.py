import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def extract_temperature_profile(thermal_image, start_point, end_point, num_points=100):
    if start_point[0] == end_point[0]:  # Handling vertically aligned points
        x = np.full(num_points, start_point[0])
        y = np.linspace(start_point[1], end_point[1], num_points)
    else:
        x = np.linspace(start_point[0], end_point[0], num_points)
        y = np.linspace(start_point[1], end_point[1], num_points)
        f = interp1d([start_point[0], end_point[0]], [start_point[1], end_point[1]], kind='linear')
        y = f(x)
    temperature_profile = [thermal_image[int(round(y_i)), int(round(x_i))] for x_i, y_i in zip(x, y)]
    distances = np.linspace(0, 5, num_points)
    return distances, temperature_profile

def plot_temperature_profiles(data_directory):
    txt_files = [f for f in os.listdir(data_directory) if f.endswith('.txt') and f != "line.txt"]
    line_file_path = os.path.join(data_directory, "line.txt")

    temperature_data = {}
    for file in txt_files:
        file_path = os.path.join(data_directory, file)
        time = int(file.split('.')[0])
        temperature_data[time] = np.loadtxt(file_path)

    line_data = pd.read_csv(line_file_path)

    plt.figure(figsize=(10, 5))
    colors = plt.cm.viridis(np.linspace(0, 1, len(temperature_data)))
    for i, (time, thermal_image) in enumerate(sorted(temperature_data.items())):
        line_row = line_data[line_data['time'] == time]
        if not line_row.empty:
            start_point = line_row[['x1', 'y1']].values.flatten()
            end_point = line_row[['x2', 'y2']].values.flatten()
            distances, temperature_profile = extract_temperature_profile(thermal_image, start_point, end_point, num_points=100)
            smoothed_temperature_profile = savgol_filter(temperature_profile, window_length=11, polyorder=2)
            plt.plot(distances, smoothed_temperature_profile, label=f'{time}s', color=colors[i])

    plt.title('Smoothed Temperature Profiles along the Line')
    plt.xlabel('Distance (cm)')
    plt.ylabel('Temperature')
    plt.legend()
    plt.grid(True)
    plt.show()

# Specify the directory containing the data files
data_directory = r"txt-inputs\lmfit_txts\txts_side\0cb_70W"
plot_temperature_profiles(data_directory)
