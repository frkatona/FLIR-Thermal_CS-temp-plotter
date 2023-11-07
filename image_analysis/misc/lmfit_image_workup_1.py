import os
import numpy as np
import matplotlib.pyplot as plt

"""
This script visualizes the temperature data in a folder of text files in a grid 
for the purposes of finding the start and endpoints of the lines to use for lmfitting simulation data.
Those start/end points are saved in a text file in the same folder as the temperature data
"""

def read_temperature_data(file_path):
    with open(file_path, 'r') as file:
        data = np.array([list(map(float, line.split('\t')[:-1])) for line in file.readlines()])
    return data

def visualize_temperature_images(folder_path):
    # List all .txt files in the folder
    txt_files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
    
    # Read and plot each image
    num_images = len(txt_files)
    cols = 3  # Number of columns in the grid of plots
    rows = (num_images + cols - 1) // cols  # Calculate the required number of rows
    
    fig, axs = plt.subplots(rows, cols, figsize=(15, 5 * rows))
    if not isinstance(axs, np.ndarray):
        axs = np.array([axs])
    axs = axs.flatten()
    
    for i, txt_file in enumerate(txt_files):
        data = read_temperature_data(os.path.join(folder_path, txt_file))
        axs[i].imshow(data, cmap='hot', aspect='auto')
        axs[i].set_title(txt_file)
        axs[i].axis('off')
    
    # Hide any remaining empty subplots
    for i in range(num_images, len(axs)):
        axs[i].axis('off')
    
    plt.tight_layout()
    plt.show()

# Specify the path to the folder containing the text files
folder_path = r"txt-inputs\lmfit_txts\txts_side\0cb_70W"

# Run the visualization
visualize_temperature_images(folder_path)