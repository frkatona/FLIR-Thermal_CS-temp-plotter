import numpy as np
from lmfit import Model, Parameters
import matplotlib.pyplot as plt

# Define the function to be fitted
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Generate synthetic data
x = np.linspace(0, 15, 301)
y = func(x, 2.5, 1.3, 0.3)
np.random.seed(0)
y_noise = y + 0.2 * np.random.randn(len(x))

# Create a model from the function
gmodel = Model(func)

# Create a parameter set with initial guesses
params = Parameters()
params.add('a', value=1)
params.add('b', value=1)
params.add('c', value=1)

# Fit the model to the data
result = gmodel.fit(y_noise, params, x=x)

# Print the fit report
print(result.fit_report())

# Plot the data and the fit
plt.figure()
plt.plot(x, y_noise, 'bo')
plt.plot(x, result.best_fit, 'r-')
plt.show()