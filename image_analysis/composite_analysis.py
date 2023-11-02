import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from scipy.ndimage import gaussian_filter

########################
######## UTILS #########
########################

def read_temperature_data(file_path):
    '''read temperature data from a text file'''
    with open(file_path, 'r') as file:
        data = [list(map(float, line.split())) for line in file]
    return np.array(data)

def display_thermal_images_with_lines(files, folder_path, direction, length, cm_per_pixel, window_title):
    '''display temperature data as thermal images with superimposed lines showing where the data will be extracted for plotting'''
    # Determine the number of images
    num_images = len(files)
    
    # Order the files by the time extracted from the filename (assumes format is "<time>.txt")
    files_sorted = sorted(files, key=lambda x: int(x.split('.')[0]))
    
    # Get the position of the line from the '5.txt' image
    file_5_path = os.path.join(folder_path, '5.txt')
    temp_data_5 = read_temperature_data(file_5_path)
    max_index_5 = np.unravel_index(np.argmax(temp_data_5, axis=None), temp_data_5.shape)

    # Setup the subplot grid
    ncols = 4  # You can adjust this as per your display requirements
    nrows = num_images // ncols + (num_images % ncols > 0)
    
    # Create a figure to hold all subplots
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 5, nrows * 5))
    fig.suptitle(window_title, fontsize=16)  # Set the overall title for the window
    
    # Flatten the axes array for easy iteration
    axes_flat = axes.flatten()
    
    # Loop through each file and plot it on the respective subplot
    for i, filename in enumerate(files_sorted):
        file_path = os.path.join(folder_path, filename)
        temp_data = read_temperature_data(file_path)
        
        # Apply a Gaussian filter for smoothing (optional)
        smoothed_data = gaussian_filter(temp_data, sigma=1)
        
        # Use the max index from '5.txt' for '0.txt' or find the new max index
        if filename == '0.txt':
            max_index = max_index_5
        else:
            max_index = np.unravel_index(np.argmax(temp_data, axis=None), temp_data.shape)
        
        # Calculate the end index based on the direction
        if direction == 'right':
            end_index = min(max_index[1] + length, temp_data.shape[1])
            profile_line = ([max_index[1], end_index], [max_index[0], max_index[0]])
        elif direction == 'down':
            end_index = min(max_index[0] + length, temp_data.shape[0])
            profile_line = ([max_index[1], max_index[1]], [max_index[0], end_index])
        
        # Plot the thermal image and the profile line
        ax = axes_flat[i]
        ax.imshow(smoothed_data, cmap='hot', interpolation='nearest')
        ax.plot(profile_line[0], profile_line[1], 'b-', linewidth=2)
        ax.set_title(f'{filename}')
    
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis('off')

    # Adjust layout to prevent overlapping
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)  # Adjust the top padding to make space for the overall title

def process_files(files, extract_path, direction, length, cm_per_pixel):
    '''process temperature FLIR image files extracted from ImageJ and extract temperature profiles'''
    profiles = []
    for filename in files:
        file_path = os.path.join(extract_path, filename)
        temp_data = read_temperature_data(file_path)
        temp_data_array = np.array(temp_data)
        
        max_index = np.unravel_index(np.argmax(temp_data_array, axis=None), temp_data_array.shape)
        if direction == 'right':
            end_index = min(max_index[1] + length, temp_data_array.shape[1])
            profile = temp_data_array[max_index[0], max_index[1]:end_index]
        elif direction == 'down':
            end_index = min(max_index[0] + length, temp_data_array.shape[0])
            profile = temp_data_array[max_index[0]:end_index, max_index[1]]
        
        cm_indices = np.arange(len(profile)) * cm_per_pixel
        profiles.append((cm_indices, profile))
    return profiles

def plot_temperature_profiles(files, profiles, title):
    '''plot the temperature profiles from the top and side views from each time as scatterplots--one for each view'''
    fig, ax = plt.subplots(figsize=(12, 6))
    cmap = LinearSegmentedColormap.from_list('custom blue', [(0, 'lightblue'), (1, 'darkblue')], N=256)
    cmap_adjusted = cmap(np.linspace(0.2, 1, 256))
    
    # Combine files and profiles into a list of tuples and sort by time extracted from filenames
    combined = sorted(zip(files, profiles), key=lambda x: int(x[0].split('.')[0]))
    
    # Use the sorted data for plotting
    for i, (filename, (cm_indices, profile)) in enumerate(combined):
        time = int(filename.split('.')[0])
        color = cmap_adjusted[int(i / len(combined) * 256)]
        ax.scatter(cm_indices, profile, s=50, color=color, label=f'{time}s', marker='o')
    
    ax.set_title(title)
    ax.set_xlabel('Distance (cm)')
    ax.set_ylabel('Temperature (°C)')
    ax.legend(title='Time (s)')
    ax.grid(True)

def combine_temperature_profiles(side_files, side_profiles, top_files, top_profiles, folder_name, cm_per_pixel_side, cm_per_pixel_top):
    '''combine the temperature profiles from the top and side views into a single CSV file for lmfit analysis alongside simulation data'''

    # Initialize lists to store data
    combined_data = []
    
    # Append data from side profiles
    for i, (cm_indices, profile) in enumerate(side_profiles):
        time = int(side_files[i].split('.')[0])  # Extract time from filename as integer
        for cm, temp in zip(cm_indices, profile):
            combined_data.append({
                'Sample_Type': folder_name,
                'View': 'Side',
                'X_cm': 0,   # X is constant for side view
                'Y_cm': cm * cm_per_pixel_side,  # Y is variable for side view
                'Time_s': time,
                'Temperature_C': temp
            })

    # Append data from top profiles
    for i, (cm_indices, profile) in enumerate(top_profiles):
        time = int(top_files[i].split('.')[0])  # Extract time from filename as integer
        for cm, temp in zip(cm_indices, profile):
            combined_data.append({
                'Sample_Type': folder_name,
                'View': 'Top',
                'X_cm': cm * cm_per_pixel_top,  # X is variable for top view
                'Y_cm': 0,   # Y is constant for top view
                'Time_s': time,
                'Temperature_C': temp
            })
    
    df = pd.DataFrame(combined_data)
    df_sorted = df.sort_values('Time_s')
    filename = f"{folder_name}_temperature_profile.csv"
    df_sorted.to_csv(filename, index=False)
    print(f"Combined CSV exported as {filename}")

########################
######## MAIN ##########
########################

# Paths to the folders containing the text files
side_path = r"txt-inputs\lmfit_txts\txts_side\1e-4_5W" 
top_path = r"txt-inputs\lmfit_txts\txts_top\1e-4_5W"   

# List all files in each folder
side_files = sorted(os.listdir(side_path))
top_files = sorted(os.listdir(top_path))

# Process the files and extract profiles
side_profiles = process_files(side_files, side_path, 'right', 50, 0.1)
top_profiles = process_files(top_files, top_path, 'down', 30, 5 / 30)

# Plot the temperature profiles
plot_temperature_profiles(side_files, side_profiles, "Side View Temperature Profiles")
plot_temperature_profiles(top_files, top_profiles, "Top View Temperature Profiles")

# Display the thermal images with lines
display_thermal_images_with_lines(side_files, side_path, 'right', 50, 0.1, "Side View")
display_thermal_images_with_lines(top_files, top_path, 'down', 30, 5 / 30, "Top View")

plt.show()

# Folder names (assuming they are the same as the sample types)
side_folder_name = os.path.basename(os.path.normpath(side_path))
top_folder_name = os.path.basename(os.path.normpath(top_path))

# Verify that folder names are the same for both views
if side_folder_name != top_folder_name:
    raise ValueError("Folder names (sample types) do not match.")
folder_name = side_folder_name  # This will be used as the sample type

# Combine the temperature profiles from both views and export to a single CSV
combine_temperature_profiles(
    side_files, side_profiles,
    top_files, top_profiles,
    folder_name,
    cm_per_pixel_side=0.1,
    cm_per_pixel_top=5 / 30
)

print("Scatterplots generated and CSV files exported.")