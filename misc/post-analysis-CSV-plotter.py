import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the path to your local folder containing the CSV files
folder_path = 'YOUR_LOCAL_PATH_TO_FOLDER'

# Logarithmic function for fitting
def log_func(x, a, b, c):
    return a + b * np.log(x + c)

# List CSV files
csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]
csv_files.sort()

# Use the 'jet' colormap to assign colors to each dataset
jet_colors = plt.cm.jet(np.linspace(0, 1, len(csv_files)))

plt.figure(figsize=(12, 8))

for idx, csv_file in enumerate(csv_files):
    file_path = os.path.join(folder_path, csv_file)
    data = pd.read_csv(file_path)
    
    # Scatter plot of the data
    plt.scatter(data["Time (s)"], data["Max Temperature (C)"], label=os.path.basename(csv_file).replace('_scatterplot_data.csv', ''),
                color=jet_colors[idx], marker='o', s=4**2*4)
    
    # Logarithmic fit and plot
    popt, _ = curve_fit(log_func, data["Time (s)"], data["Max Temperature (C)"], maxfev=5000)
    x_vals = np.linspace(min(data["Time (s)"]), max(data["Time (s)"]), 1000)
    y_vals = log_func(x_vals, *popt)
    plt.plot(x_vals, y_vals, color=jet_colors[idx], linewidth=1)

plt.xlabel("Time (s)")
plt.ylabel("Max Temperature (C)")
plt.title("Max Temperature vs Time (Log Fit)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
