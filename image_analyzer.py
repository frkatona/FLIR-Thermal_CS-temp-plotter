import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os

def moving_average(data, window_size=5):
    """Compute moving average."""
    window = np.ones(int(window_size))/float(window_size) # boxcar kernel
    # window = np.exp(-np.linspace(-1, 1, window_size)**2) # Gaussian kernel (broken)
    return np.convolve(data, window, mode='valid')

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
scale_bar_image = Image.open(r"image-inputs\0cb_s182\000_310_197.jpg")
left, top, right, bottom = 65, 453, 574, 458
cropped_img = scale_bar_image.crop((left, top, right, bottom))
img_array = np.array(cropped_img)
average_colors = np.mean(img_array, axis=0)

# Read all images in the folder
folder_path = r"image-inputs\0cb_s184"
image_files = [f for f in os.listdir(folder_path) if f.endswith(".jpg")]
image_files.sort()  # To ensure consistent ordering

# Extract temperatures along y=245 for each image
temperatures_data_reversed = []
temperatures_at_262 = []

for img_file in image_files:
    img_path = os.path.join(folder_path, img_file)
    time_heated, max_temp, min_temp = extract_temp_values_from_filename(img_file)
    temperatures = get_averaged_temperatures_along_y(img_path, 245, max_temp, min_temp)
    temperatures_data_reversed.append((time_heated, temperatures))
    
    # Extracting temperature at x=262
    temperature_262 = temperatures[262]
    temperatures_at_262.append((time_heated, temperature_262))

# Sort data by heating time
temperatures_data_reversed.sort(key=lambda x: x[0])

# Plot smoothed temperature profiles along y=245
cmap = plt.get_cmap("viridis", len(temperatures_data_reversed))
plt.figure(figsize=(12, 6))

for index, (time_heated, temps) in enumerate(temperatures_data_reversed):
    smoothed_temps = moving_average(temps, window_size=5)
    plt.plot(smoothed_temps, label=f"{time_heated}s", color=cmap(index))

plt.title(folder_path)
plt.xlim(100, 700)
plt.ylim(20, 120)
plt.xlabel("pixel position")
plt.ylabel("temperature (°C)")
plt.legend(title="irradiation time")
plt.grid(True)
plt.tight_layout()
plt.show()

# Scatterplot of temperatures at x=262 against irradiation time
times = [item[0] for item in temperatures_at_262]
temps_262 = [item[1] for item in temperatures_at_262]

# Compute the linear trendline
slope, intercept = np.polyfit(times[0:7], temps_262[0:7], 1)
trendline = [slope * x + intercept for x in times]

plt.figure(figsize=(8, 6))

# Plot the points used for the trendline in a darker color
plt.scatter(times[0:7], temps_262[0:7], color='red', label="linear fit points")

# Plot the other points in a lighter color
plt.scatter(times[7:], temps_262[7:], color='blue', label="excluded points")

# Plotting the trendline
plt.plot(times, trendline, color='red', linestyle='--')  

# Display the equation on the plot
equation = f"y = {slope:.3f}x + {intercept:.3f}"
plt.text(min(times), max(temps_262), equation, color='red', verticalalignment='top')

plt.ylim(20, 90)
plt.xlim(0, 800)
plt.title(f"temperature at pixel position 262 for {folder_path}")
plt.xlabel("irradiation time (s)")
plt.ylabel("temperature (°C)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Displaying each image with the overlaid cross-sectional line
fig, ax = plt.subplots(figsize=(14, 10))
y_value = 245  # This is the y-value of the cross-section

for img_file in image_files:
    img_path = os.path.join(folder_path, img_file)
    img = Image.open(img_path)
    ax.imshow(img, alpha=0.7)
    
    # Draw the cross-sectional line
    ax.plot([0, img.width], [y_value, y_value], color='red', linestyle='--')

# Setting title and other properties
ax.set_title(f"Images with Cross-sectional Line from {folder_path}")
ax.axis('off')  # Hide the axis values and ticks
plt.tight_layout()
plt.show()