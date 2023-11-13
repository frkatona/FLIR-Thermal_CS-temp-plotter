import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Replace with the actual paths to your CSV files
experimental_path = r'exports\CSVs\lmfit_consolidated\1e-6_70W_temperature_profile.csv'
simulation_path = r'exports\CSVs\simulated_toprow\70W_1e-06_top-row.csv'

# Read the data
df_experimental = pd.read_csv(experimental_path)
df_simulation = pd.read_csv(simulation_path)

# Process exp. data
df_experimental_top = df_experimental[df_experimental['View'] == 'Top'].pivot(index='Time_s', columns='X_cm', values='Temperature_C')
num_points_file1 = len(df_experimental_top.columns)
x_vales_experimental = np.linspace(0, 5, num_points_file1)
# df_experimental_top['cm'] = x_vales_experimental
# df_experimental_top.to_csv('experimental_top_1e-6_70W.csv')

# Process sim. data
top_data2 = df_simulation.set_index('Distance_cm').T
num_points_file2 = len(top_data2.columns)
x_values_file2 = np.linspace(0, 5, num_points_file2)

# Set 'rainbow' color map
cmap = cm.Reds
colors = cmap(np.linspace(0, 1, max(len(df_experimental_top.index), len(top_data2.index))))

# Create the plot
plt.figure(figsize=(14, 8))

# Scatter plot for data from file1
for i, (time, row) in enumerate(df_experimental_top.iterrows()):
    color = colors[i % len(colors)]
    plt.scatter(x_vales_experimental, row.values, color=color, alpha=0.9, label=f'experimental {time}s', marker='o')

# Line plot for data from file2 with thicker lines
for i, (time, row) in enumerate(top_data2.iterrows()):
    color = colors[i % len(colors)]
    plt.plot(x_values_file2, row.values, color=color, alpha=0.7, linewidth=3, label=f'simulation {time}')

# Customize the plot
plt.xlabel('Position (cm)')
plt.ylabel('Temperature (°C)')
plt.title('Top Row Temperature Profiles from Two Files (Rainbow Color Map)')
plt.legend()
plt.show()