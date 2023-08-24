import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap
import csv

## folder paths ##
folder_path = r'C:\Users\antho\OneDrive\1 - Thesis Data - Master\1_FLIR Images\001_ThesisRevisions\top-view-attempts\1e-6_high-P\text'
folder_name = os.path.basename(os.path.normpath(folder_path))
output_csv_path = f'{folder_name}_scatterplot_data.csv'

## colormap and data range customization ##
cmap = "inferno"
vmin = 20
vmax = 140

colors = ["#251645", "#0773B1"]
cm = LinearSegmentedColormap.from_list("custom", colors, N=len(os.listdir(folder_path)))

## convert values from text files to thermal images ##
def display_individual_image(local_file_path, cmap, vmin, vmax):
    data = np.loadtxt(local_file_path)
    max_position = np.unravel_index(data.argmax(), data.shape)
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax)
    title_text = os.path.basename(local_file_path).split('.')[0] + "s"
    ax.text(0.05, 0.95, title_text, transform=ax.transAxes, color='white', verticalalignment='top', fontsize=12)
    fig.colorbar(im, ax=ax, ticks=[vmin, vmax])
    ax.axhline(max_position[0], color='white', xmin=1/3-1/60, xmax=2/3+1/60)
    ax.axis('off')
    plt.tight_layout()

txt_files = [file for file in os.listdir(folder_path) if file.endswith('.txt')]
txt_files.sort(key=lambda x: int(x.split('.')[0]))

for txt_file in txt_files:
    file_path = os.path.join(folder_path, txt_file)
    display_individual_image(file_path, cmap, vmin, vmax)

def interpolate_data(x, y, num_points=300):
    '''smooth line plot'''
    f = interp1d(x, y, kind='cubic')
    xnew = np.linspace(x.min(), x.max(), num_points)
    ynew = f(xnew)
    return xnew, ynew

def log_func(x, a, b, c):
    '''log fit'''
    return a + b * np.log(x + c)

## create text file list ##
txt_files = [file for file in os.listdir(folder_path) if file.endswith('.txt')]
txt_files.sort(key=lambda x: int(x.split('.')[0]))

## plot temperature profiles ##
fig, ax = plt.subplots(figsize=(10, 6))
for idx, txt_file in enumerate(txt_files):
    file_path = os.path.join(folder_path, txt_file)
    data = np.loadtxt(file_path)
    y_values = data[np.unravel_index(data.argmax(), data.shape)[0], data.shape[1]//3-1:2*data.shape[1]//3+1]
    x_values_cm = np.linspace((y_values.shape[0]//3 - 1)*5/30, (2*y_values.shape[0]//3 + 1)*5/30, len(y_values))
    x_smooth, y_smooth = interpolate_data(x_values_cm, y_values)
    ax.plot(x_smooth, y_smooth, color=cm(idx), label=os.path.basename(txt_file).split('.')[0] + "s")

ax.set_xlabel("Position (cm)")
ax.set_ylabel("Temperature (C)")
ax.legend()

## find max temperature, plot vs time, export as csv ##
max_values = [np.max(np.loadtxt(os.path.join(folder_path, txt_file))) for txt_file in txt_files]
times = [int(txt_file.split('.')[0]) for txt_file in txt_files]
popt, _ = curve_fit(log_func, times, max_values)
y_pred = log_func(np.array(times), *popt)

with open(output_csv_path, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    csvwriter.writerow(["Time (s)", "Max Temperature (C)"]) # header
    
    for t, max_val in zip(times, max_values): # data
        csvwriter.writerow([t, max_val])

print(f"Data exported to {output_csv_path}")

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(times, max_values, color='blue', marker='o', label='Data')
ax.plot(times, y_pred, color='red', label=f'Log Fit: y = {popt[0]:.2f} + {popt[1]:.2f} * log(x + {popt[2]:.2f})')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Max Temperature (C)")
ax.legend()

## print fit equation and error ##
fit_equation = f"y = {popt[0]:.2f} + {popt[1]:.2f} * log(x + {popt[2]:.2f})"
fit_error = np.sqrt(np.mean((max_values - y_pred)**2))
print(f"Fit Equation: {fit_equation}")
print(f"Fit Error: ±{fit_error:.2f} C")

plt.show()