# create a 5x5 array of random numbers
import numpy as np

# create a 5x5 array of random numbers
array = np.random.rand(5, 5)
print(array)

# print a 1x5 array of the average of the top 3 rows
print(np.average(array[0:3], axis=0))