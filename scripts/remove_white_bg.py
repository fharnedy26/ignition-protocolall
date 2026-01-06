"""
Remove white background from PNG image, making it transparent.
"""
from pathlib import Path
from PIL import Image
import numpy as np

def remove_white_background(image: Image.Image, threshold: int = 240) -> Image.Image:
    """
    Remove white/light background by making white and near-white pixels transparent.
    
    Args:
        image: PIL Image object
        threshold: RGB threshold (0-255) - pixels above this are considered white/light
    
    Returns:
        Image with RGBA mode and transparent white background
    """
    # Convert to RGBA if not already
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    # Convert to numpy array for processing
    data = np.array(image)
    
    # Strategy 1: Remove pure white and near-white pixels
    white_mask = (data[:, :, 0] > threshold) & \
                 (data[:, :, 1] > threshold) & \
                 (data[:, :, 2] > threshold)
    
    # Strategy 2: Also remove pixels where the average RGB is very high
    avg_brightness = (data[:, :, 0].astype(np.float32) + 
                     data[:, :, 1].astype(np.float32) + 
                     data[:, :, 2].astype(np.float32)) / 3.0
    bright_mask = avg_brightness > (threshold + 20)
    
    # Strategy 3: Remove pixels that are very close to white (low saturation)
    rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
    low_sat_mask = (saturation < 30) & (rgb_max > threshold)
    
    # Combine all masks
    combined_mask = white_mask | bright_mask | low_sat_mask
    
    # Set alpha channel to 0 for white/light pixels
    data[:, :, 3] = np.where(combined_mask, 0, data[:, :, 3])
    
    # Convert back to PIL Image
    return Image.fromarray(data)

if __name__ == "__main__":
    import sys
    
    # Input file
    input_path = Path("exports/Ignition Protocol - Book Cover_page_1.png")
    
    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    
    if not input_path.exists():
        print(f"Error: File not found: {input_path}")
        sys.exit(1)
    
    # Output path (overwrite original)
    output_path = input_path
    
    print(f"Processing {input_path.name}...")
    print(f"  Input size: {Image.open(input_path).size}")
    
    # Load image
    image = Image.open(input_path)
    
    # Remove white background
    print("  Removing white background...")
    image = remove_white_background(image, threshold=240)
    
    # Save result (overwrite original)
    image.save(output_path, 'PNG')
    
    print(f"  Output size: {image.size}")
    print(f"  Saved: {output_path}")
    print(f"\nProcessing complete! White background removed.")



