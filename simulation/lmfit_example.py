from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt


# Define the Gaussian function
def gaussian(x, amplitude, center, sigma):
    return amplitude * np.exp(-(x - center) ** 2 / (2 * sigma ** 2))

# Create a Model based on the Gaussian function
gaussian_model = Model(gaussian)

# Generate synthetic data
np.random.seed(0)  # for repeatability
x = np.linspace(-10, 10, 200)
true_amplitude = 1
true_center = 0
true_sigma = 2
noise = np.random.normal(scale=0.1, size=x.size)
y = gaussian(x, true_amplitude, true_center, true_sigma) + noise

# Initialize parameters: give initial guesses for the fit.
params = gaussian_model.make_params(amplitude=1, center=0, sigma=1)

# Perform the fit
result = gaussian_model.fit(y, params, x=x)

# Plot the data and the fit
plt.figure(figsize=(10, 5))
plt.plot(x, y, 'b', label='data with noise')
plt.plot(x, result.best_fit, 'r-', label='best fit')
plt.legend(loc='best')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Gaussian Model Fit')
plt.show()

# Print out the fit statistics and optimized parameters
fit_report = result.fit_report()