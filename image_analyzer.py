import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os

def moving_average(data, window_size=5):
    """Compute moving average."""
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def extract_temp_values_from_filename(filename):
    values = filename.split("_")
    time_heated = int(values[0])
    max_temp = int(values[1]) / 10
    min_temp = int(values[2].split(".")[0]) / 10
    return time_heated, max_temp, min_temp

def get_distance_reversed(color):
    differences = np.sqrt(np.sum((average_colors - color) ** 2, axis=1))
    x_position = np.argmin(differences)
    return 1 - (x_position / (cropped_img.width - 1))

def get_mapped_distance_reversed(color, left_val, right_val):
    relative_distance = get_distance_reversed(color)
    return left_val + relative_distance * (right_val - left_val)

def get_averaged_temperatures_along_y(img_path, y_value, max_temp, min_temp):
    """
    Get the temperatures along a specific y value in the image and average with its neighbors.
    """
    img = Image.open(img_path)
    
    # Extract the colors along the specified y value and its neighboring pixels
    main_colors = [img.getpixel((x, y_value)) for x in range(img.width)]
    upper_colors = [img.getpixel((x, y_value - 1)) for x in range(img.width)]
    lower_colors = [img.getpixel((x, y_value + 1)) for x in range(img.width)]
    
    # Calculate the average colors
    avg_colors = [(np.array(main) + np.array(upper) + np.array(lower)) / 3 for main, upper, lower in zip(main_colors, upper_colors, lower_colors)]
    
    # Convert the average colors to temperatures using the reversed scale bar data structure
    temperatures = [get_mapped_distance_reversed(color, min_temp, max_temp) for color in avg_colors]
    return temperatures
    

# Assuming the scale bar is from the image named "300_905_212.jpg"
scale_bar_image = Image.open("0cb_s184/300_905_212.jpg")
left, top, right, bottom = 65, 453, 574, 458
cropped_img = scale_bar_image.crop((left, top, right, bottom))
img_array = np.array(cropped_img)
average_colors = np.mean(img_array, axis=0)

# Read all images in the folder
folder_path = "0cb_s184"
image_files = [f for f in os.listdir(folder_path) if f.endswith(".jpg")]
image_files.sort()  # To ensure consistent ordering

# Extract temperatures along y=245 for each image
temperatures_data_reversed = []
for img_file in image_files:
    img_path = os.path.join(folder_path, img_file)
    time_heated, max_temp, min_temp = extract_temp_values_from_filename(img_file)
    temperatures = get_averaged_temperatures_along_y(img_path, 245, max_temp, min_temp)  # <-- This line is changed
    temperatures_data_reversed.append((time_heated, temperatures))

# Sort data by heating time
temperatures_data_reversed.sort(key=lambda x: x[0])

# Plot with smoothed data
cmap = plt.get_cmap("viridis", len(temperatures_data_reversed))
plt.figure(figsize=(12, 6))

for index, (time_heated, temps) in enumerate(temperatures_data_reversed):
    smoothed_temps = moving_average(temps, window_size=5)
    plt.plot(smoothed_temps, label=f"{time_heated}s", color=cmap(index))

plt.title(folder_path)
plt.ylim(20, 160)
plt.xlabel("Pixel Position")
plt.ylabel("temperature (°C)")
plt.legend(title="Irradiation Time")
plt.grid(True)
plt.tight_layout()
plt.show()