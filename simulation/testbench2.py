from PIL import Image, ImageDraw, ImageFont
import os

# Path to your folder containing the images and the desired output GIF
image_folder_path = r'C:\Users\antho\Desktop\Desktop_tmp\FLIR_tmp\1e-4'
output_gif_path = 'output_gif_file.gif'

# Get the list of images from the folder
image_files = sorted([os.path.join(image_folder_path, f) for f in os.listdir(image_folder_path) if f.endswith('.jpg')])

# Load and annotate images
annotated_images = []
for image_file in image_files:
    image = Image.open(image_file)
    
    # Annotate the image with its filename (without the extension)
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()
    title = os.path.splitext(os.path.basename(image_file))[0]
    text_width, text_height = draw.textsize(title, font=font)
    draw.text(((image.width - text_width) / 2, 10), title, font=font, fill="white")
    
    annotated_images.append(image)

# Create a GIF using the annotated images
annotated_images[0].save(output_gif_path, save_all=True, append_images=annotated_images[1:], loop=0, duration=200)

print(f'GIF created at: {output_gif_path}')
