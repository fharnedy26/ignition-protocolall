"""
Process book cover PNG: remove white background and resize to A4 size.
A4 dimensions: 210mm x 297mm
At 300 DPI: 2480 x 3508 pixels
"""
from pathlib import Path
from PIL import Image
import numpy as np

def remove_white_background(image: Image.Image, threshold: int = 200) -> Image.Image:
    """
    Remove white/light background by making white and near-white pixels transparent.
    Uses multiple strategies to ensure complete white removal.
    
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
    # Check if all RGB channels are above threshold (white/light colors)
    white_mask = (data[:, :, 0] > threshold) & \
                 (data[:, :, 1] > threshold) & \
                 (data[:, :, 2] > threshold)
    
    # Strategy 2: Also remove pixels where the average RGB is very high
    # This catches slightly off-white colors
    avg_brightness = (data[:, :, 0].astype(np.float32) + 
                     data[:, :, 1].astype(np.float32) + 
                     data[:, :, 2].astype(np.float32)) / 3.0
    bright_mask = avg_brightness > (threshold + 20)
    
    # Strategy 3: Remove pixels that are very close to white (low saturation)
    # Calculate saturation: max - min of RGB channels
    rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
    # Low saturation + high brightness = white/light gray
    low_sat_mask = (saturation < 30) & (rgb_max > threshold)
    
    # Combine all masks - remove any pixel that matches any condition
    combined_mask = white_mask | bright_mask | low_sat_mask
    
    # Set alpha channel to 0 for white/light pixels
    data[:, :, 3] = np.where(combined_mask, 0, data[:, :, 3])
    
    # Convert back to PIL Image
    return Image.fromarray(data)

def resize_to_a4(image: Image.Image, dpi: int = 300) -> Image.Image:
    """
    Resize image to A4 dimensions at specified DPI.
    
    Args:
        image: PIL Image object
        dpi: Resolution in DPI (default 300)
    
    Returns:
        Resized image to A4 dimensions
    """
    # A4 dimensions in mm
    a4_width_mm = 210
    a4_height_mm = 297
    
    # Convert to pixels
    mm_to_inch = 1 / 25.4
    a4_width_px = int(a4_width_mm * mm_to_inch * dpi)
    a4_height_px = int(a4_height_mm * mm_to_inch * dpi)
    
    # Create new A4-sized image with transparent background
    a4_image = Image.new('RGBA', (a4_width_px, a4_height_px), (0, 0, 0, 0))
    
    # Calculate scaling to fit image within A4 while maintaining aspect ratio
    img_width, img_height = image.size
    scale_w = a4_width_px / img_width
    scale_h = a4_height_px / img_height
    scale = min(scale_w, scale_h)  # Use smaller scale to fit within A4
    
    # Calculate new dimensions
    new_width = int(img_width * scale)
    new_height = int(img_height * scale)
    
    # Resize image
    resized = image.resize((new_width, new_height), Image.Resampling.LANCZOS)
    
    # Center the image on A4 canvas
    x_offset = (a4_width_px - new_width) // 2
    y_offset = (a4_height_px - new_height) // 2
    
    # Paste resized image onto A4 canvas
    a4_image.paste(resized, (x_offset, y_offset), resized)
    
    return a4_image

def process_book_cover(input_path: str, output_path: str = None, dpi: int = 300, white_threshold: int = 240):
    """
    Process book cover: remove white background and resize to A4.
    
    Args:
        input_path: Path to input PNG file
        output_path: Path to output PNG file (default: same name with '_processed' suffix)
        dpi: Resolution in DPI (default 300)
        white_threshold: RGB threshold for white detection (default 240)
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Image not found: {input_path}")
    
    if output_path is None:
        output_path = input_path.parent / f"{input_path.stem}_processed.png"
    else:
        output_path = Path(output_path)
    
    print(f"Processing {input_path.name}...")
    print(f"  Input size: {Image.open(input_path).size}")
    
    # Load image
    image = Image.open(input_path)
    
    # Remove white background
    print("  Removing white background...")
    image = remove_white_background(image, threshold=white_threshold)
    
    # Resize to A4
    print(f"  Resizing to A4 size at {dpi} DPI...")
    a4_image = resize_to_a4(image, dpi=dpi)
    
    # Save result
    a4_image.save(output_path, 'PNG', dpi=(dpi, dpi))
    
    print(f"  Output size: {a4_image.size}")
    print(f"  Saved: {output_path}")
    print(f"\nProcessing complete!")

if __name__ == "__main__":
    import sys
    
    # Default input path
    input_path = "exports/book-cover.png"
    
    # Allow override via command line
    if len(sys.argv) > 1:
        input_path = sys.argv[1]
    
    # Output path
    output_path = None
    if len(sys.argv) > 2:
        output_path = sys.argv[2]
    
    # DPI (default 300)
    dpi = 300
    if len(sys.argv) > 3:
        dpi = int(sys.argv[3])
    
    # White threshold (default 240)
    white_threshold = 240
    if len(sys.argv) > 4:
        white_threshold = int(sys.argv[4])
    
    process_book_cover(input_path, output_path, dpi, white_threshold)

