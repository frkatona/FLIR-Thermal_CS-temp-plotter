import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from colorsys import rgb_to_hls, hls_to_rgb

def reduce_saturation(rgb, fraction):
    h, l, s = rgb_to_hls(*rgb[:3])
    s *= fraction
    return hls_to_rgb(h, l, s)

def generate_colormap(base_color, num_colors):
    color_sub = num_colors - 1
    if color_sub == 0:
        return [base_color]
    return [reduce_saturation(base_color, 0.5 + 0.5 * (i / (color_sub))) for i in range(num_colors)]

def plot_concentric_circles(csv_folder_path):
    csv_files = [file for file in os.listdir(csv_folder_path) if file.endswith('.csv')]
    csv_files.sort()

    # Identify the shared times across all CSVs
    shared_times = None
    for csv_file in csv_files:
        file_path = os.path.join(csv_folder_path, csv_file)
        data = pd.read_csv(file_path)
        if shared_times is None:
            shared_times = set(data['Time (s)'])
        else:
            shared_times &= set(data['Time (s)'])

    total_csvs = len(csv_files)
    jet_colors = plt.get_cmap('jet', total_csvs)

    # For each shared time, plot the circles in separate plots
    for time in shared_times:
        plt.figure(figsize=(10, 10))
        ax = plt.gca()

        max_threshold = 0
        for idx, csv_file in enumerate(csv_files):
            file_path = os.path.join(csv_folder_path, csv_file)
            data = pd.read_csv(file_path)
            distances_at_time = data[data['Time (s)'] == time]["68.2% Distance (cm)"]

            base_color = jet_colors(idx / (total_csvs - 1))
            colors = generate_colormap(base_color, len(distances_at_time))

            for dist_idx, dist in enumerate(distances_at_time):
                circle = plt.Circle((0, 0), dist, color=colors[dist_idx], fill=False, linewidth=3)
                ax.add_patch(circle)
                max_threshold = max(max_threshold, dist)

            ax.plot([], [], color=base_color, linewidth=3, label=f"{os.path.splitext(csv_file)[0]}")

        # Set up the plot layout
        tick_interval = 0.5  # change to 0.2 if you prefer
        ticks = np.arange(0, max_threshold + tick_interval, tick_interval)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlim(0, max_threshold)
        ax.set_ylim(0, max_threshold)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.set_aspect('equal', 'box')
        ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1))
        plt.title(f"Time: {time}s")
        plt.tight_layout()
        plt.show()

folder_path = 'exports\csv_outputs'
plot_concentric_circles(folder_path)
