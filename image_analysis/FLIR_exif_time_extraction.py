import os
from PIL import Image, ExifTags
import csv

def extract_exif_from_image(image_path):
    """Extract all Exif data from an image."""
    with Image.open(image_path) as img:
        exif_data = img._getexif()
        if exif_data:
            return {ExifTags.TAGS.get(tag, tag): value for tag, value in exif_data.items()}
    return {}

def generate_exif_csv_from_images_in_folder(folder_path):
    """Generate a CSV of filenames and all their Exif data."""
    # List all JPEG files in the folder
    image_files = [f for f in os.listdir(folder_path) if f.lower().endswith('.jpg')]
    
    # Extract Exif data for all images
    all_exif_data = {}
    for image_file in image_files:
        image_path = os.path.join(folder_path, image_file)
        all_exif_data[image_file] = extract_exif_from_image(image_path)
    
    # Find all unique Exif tags across the images
    all_tags = set()
    for exif_data in all_exif_data.values():
        all_tags.update(exif_data.keys())
    all_tags = sorted(list(all_tags))
    
    # Write to CSV
    csv_filename = os.path.basename(folder_path) + ".csv"
    with open(csv_filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Filename"] + all_tags)
        for image_file, exif_data in all_exif_data.items():
            writer.writerow([image_file] + [exif_data.get(tag, "") for tag in all_tags])
    print(f"CSV saved as {csv_filename}")

# Run the function
folder_path = r"C:\Users\antho\Downloads\FLIR_tmp-20230922T215314Z-001\FLIR_tmp\1e-4"  # Replace with your folder path
generate_exif_csv_from_images_in_folder(folder_path)