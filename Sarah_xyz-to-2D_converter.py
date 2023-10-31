import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_xyz_data(file_path):
    xyz_data = []
    with open(file_path, 'r') as file:
        xyz_data_started = False
        for line in file:
            if line.startswith("#"):
                xyz_data_started = True
                continue
            if xyz_data_started:
                parts = line.split()
                if len(parts) >= 3:
                    x, y, z = parts[:3]
                    try:
                        x, y, z = float(x), float(y), float(z)
                        xyz_data.append((x, y, z))
                    except ValueError:
                        print("Could not convert line to floats:", line)
                        pass
    return xyz_data

def convert_to_rectangular_array(xyz_data):
    if not xyz_data:
        print("No data to convert")
        return None
    
    xyz_array = np.array(xyz_data)
    x_min, x_max = int(np.min(xyz_array[:, 0])), int(np.max(xyz_array[:, 0]))
    y_min, y_max = int(np.min(xyz_array[:, 1])), int(np.max(xyz_array[:, 1]))
    rows, cols = y_max - y_min + 1, x_max - x_min + 1
    rectangular_array = np.full((rows, cols), np.nan)
    for x, y, z in xyz_data:
        x_idx, y_idx = int(x - x_min), int(y - y_min)
        rectangular_array[y_idx, x_idx] = z
    return rectangular_array

def visualize_3d_data(rectangular_array):
    if rectangular_array is None:
        print("No data to visualize")
        return
    
    x = np.arange(rectangular_array.shape[1])
    y = np.arange(rectangular_array.shape[0])
    X, Y = np.meshgrid(x, y)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, rectangular_array, cmap='viridis', edgecolor='none')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Surface Plot of XYZ Data')
    plt.colorbar(surf)
    plt.show()

# Path to the .xyz file
file_path = r"C:\Users\antho\Desktop\Desktop_tmp\laser_thin.csv"

# Read the XYZ data
xyz_data = read_xyz_data(file_path)

if not xyz_data:
    print("No XYZ data read from file")

# Convert the XYZ data to a rectangular array
rectangular_array = convert_to_rectangular_array(xyz_data)

# Save the rectangular array as a .csv file
if rectangular_array is not None:
    output_csv_path = "sarahtest.csv"
    np.savetxt(output_csv_path, rectangular_array, delimiter=",")
    print("Data saved to", output_csv_path)

# Visualize the data in 3D
visualize_3d_data(rectangular_array)