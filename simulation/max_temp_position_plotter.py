import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from lmfit import Model
import numpy as np

# Define the exponential model
def exponential(x, amplitude, decay, C):
    return amplitude * np.exp(-decay * x) + C

csv_path = Path('exports/CSVs/max_temp_positions/70W_1e-06_max-temp-positions_50.csv')
df = pd.read_csv(csv_path)
indices = int(csv_path.stem.split('_')[-1])
df['Position'] = df['Position'] / indices * 5 # convert from m to cm
sample = csv_path.stem.split("_")[0:2]
print(f"number of indices in file = {indices}")

# Create the model
model = Model(exponential)

# Set initial parameter values
params = model.make_params(amplitude=18, decay=0.46, C=181)

# Fit the model to the data
result = model.fit(df['Position'], x=df['Time_s'], params=params)

# Generate new x-values for interpolation
x_interp = np.linspace(df['Time_s'].min(), df['Time_s'].max(), 100)

# Calculate interpolated y-values using the fitted model
y_interp = result.eval(x=x_interp)

# Print the equation and fit report
print(result.model)
print(result.fit_report())

# Plot the data, fitted curve, and interpolated curve
plt.figure(figsize=(8, 6))
plt.scatter(df['Time_s'], df['Position'], label='Data')
plt.plot(x_interp, y_interp, 'r--', label='Interpolated')
plt.xlabel('t (s)')
plt.ylabel('height (cm)')
plt.title(f'Max Temperature Position vs. Irradiation Time for P = {sample[0]}, loading = {sample[1]}')

# Add the exponential equation and decay error to the graph
equation = f"Equation: {result.params['amplitude'].value:.2f} * exp(-{result.params['decay'].value:.2f} * x) + {result.params['C'].value:.2f}"
decay_error = f"Decay Error: {result.params['decay'].stderr / result.params['decay'] * 100:.2f}%"

plt.text(0.1, 0.9, equation, transform=plt.gca().transAxes)
plt.text(0.1, 0.85, decay_error, transform=plt.gca().transAxes)

plt.legend()
plt.show()