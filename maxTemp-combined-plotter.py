import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

# Load the datasets
df1 = pd.read_csv('exports\CSVs\max_temp_positions/70W_1e-07_max-temp-positions_200.csv')
df2 = pd.read_csv('exports\CSVs\max_temp_positions/70W_1e-05_max-temp-positions_200.csv')
df3 = pd.read_csv('exports\CSVs\max_temp_positions/70W_1e-06_max-temp-positions_200.csv')  # Updated file path

# Subtract 5 from the y-axis values
df1['Position_cm'] = 5 - df1['Position'] / 40
df2['Position_cm'] = 5 - df2['Position'] / 40
df3['Position_cm'] = 5 - df3['Position'] / 40

def exp_decay(x, a, b, c):
    return -a * np.exp(-b * x) + c

# Initial guesses: a = max of Position_cm, b = 0.01 (assuming a slow decay), c = min of Position_cm (as the baseline)

# # Adjusting initial guesses based on the dataset
initial_guess1 = [0.15, 0.15, 0.2]
# initial_guess2 = [df2['Position_cm'].max(), 0.01, df2['Position_cm'].min()]
# initial_guess3 = [df3['Position_cm'].max(), 0.01, df3['Position_cm'].min()]

# # Fit exponential decay to each dataset with initial guesses
# params1, _ = curve_fit(exp_decay, df1['Time_s'], df1['Position_cm'], p0=initial_guess1)
# params2, _ = curve_fit(exp_decay, df2['Time_s'], df2['Position_cm'], p0=initial_guess2)
# params3, _ = curve_fit(exp_decay, df3['Time_s'], df3['Position_cm'], p0=initial_guess3)


# Fit exponential decay to each dataset
params1, _ = curve_fit(exp_decay, df1['Time_s'], df1['Position_cm'])
params2, _ = curve_fit(exp_decay, df2['Time_s'], df2['Position_cm'])
params3, _ = curve_fit(exp_decay, df3['Time_s'], df3['Position_cm'])

# Generate x values for plotting the fitted curves
x_values = np.linspace(0, max(df1['Time_s'].max(), df2['Time_s'].max(), df3['Time_s'].max()), 100)

# Plotting
plt.figure(figsize=(16, 10))

color_low = '#fb5607'
color_mid = '#0096c7'
color_high = '#90be6d'
size = 600

plt.scatter(df1['Time_s'], df1['Position_cm'], color=color_low, label='70W 1e-07', marker='o', s=size)
plt.scatter(df3['Time_s'], df3['Position_cm'], color=color_mid, label='70W 1e-06', marker='o', s=size)
plt.scatter(df2['Time_s'], df2['Position_cm'], color=color_high, label='70W 1e-05', marker='o', s=size)

# Plot the fitted curves
linewidth = 5
plt.plot(x_values, exp_decay(x_values, *params1), color=color_low, linestyle='-', linewidth = linewidth)
plt.plot(x_values, exp_decay(x_values, initial_guess1[0], initial_guess1[1], initial_guess1[2]), color=color_high, linestyle='-', linewidth = linewidth)
# plt.plot(x_values, exp_decay(x_values, *params2), color=color_high, linestyle='--')
plt.plot(x_values, exp_decay(x_values, *params3), color=color_mid, linestyle='-', linewidth = linewidth)

plt.xlabel('time /s', fontsize=40)
plt.ylabel('depth of max T (cm)', fontsize=40)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(direction='out', length=6, width=2)

plt.legend(fontsize=20)
plt.grid(False)

print("low: ", params1)
print("mid: ", params3)
print("high: ", params2)
# print("high: ", initial_guess1[0], initial_guess1[1], initial_guess1[2])

plt.show()