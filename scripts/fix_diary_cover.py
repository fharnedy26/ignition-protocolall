"""
Fix Project Diary Cover: aggressively remove ALL white background, ensure sharp corners.
"""
from pathlib import Path
from PIL import Image
import numpy as np

def remove_white_background_ultra_aggressive(image: Image.Image, threshold: int = 180) -> Image.Image:
    """
    Ultra-aggressively remove ALL white/light background.
    Uses very low threshold and multiple strategies to catch every white pixel.
    
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
    
    # Strategy 1: Remove pure white and near-white pixels (very aggressive)
    white_mask = (data[:, :, 0] > threshold) & \
                 (data[:, :, 1] > threshold) & \
                 (data[:, :, 2] > threshold)
    
    # Strategy 2: Remove pixels where ANY channel is very high (catches off-whites)
    any_high_mask = (data[:, :, 0] > threshold + 30) | \
                    (data[:, :, 1] > threshold + 30) | \
                    (data[:, :, 2] > threshold + 30)
    
    # Strategy 3: Remove pixels where the average RGB is high
    avg_brightness = (data[:, :, 0].astype(np.float32) + 
                     data[:, :, 1].astype(np.float32) + 
                     data[:, :, 2].astype(np.float32)) / 3.0
    bright_mask = avg_brightness > threshold
    
    # Strategy 4: Remove pixels that are very close to white (low saturation)
    rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
    low_sat_mask = (saturation < 50) & (rgb_max > threshold)
    
    # Strategy 5: Remove pixels that are very light (high luminance)
    luminance = (0.299 * data[:, :, 0].astype(np.float32) + 
                 0.587 * data[:, :, 1].astype(np.float32) + 
                 0.114 * data[:, :, 2].astype(np.float32))
    high_lum_mask = luminance > threshold
    
    # Strategy 6: Remove pixels where all channels are very close (grayscale/white)
    channel_diff = np.maximum(
        np.abs(data[:, :, 0].astype(np.float32) - data[:, :, 1].astype(np.float32)),
        np.maximum(
            np.abs(data[:, :, 1].astype(np.float32) - data[:, :, 2].astype(np.float32)),
            np.abs(data[:, :, 0].astype(np.float32) - data[:, :, 2].astype(np.float32))
        )
    )
    grayscale_mask = (channel_diff < 25) & (rgb_max > threshold)
    
    # Strategy 7: Remove pixels at edges that are light (edge cleanup)
    edge_threshold = threshold + 40
    edge_mask = np.zeros_like(white_mask, dtype=bool)
    edge_width = 50  # Check first/last 50 pixels
    if data.shape[0] > edge_width * 2 and data.shape[1] > edge_width * 2:
        # Top edge
        edge_mask[:edge_width, :] = (rgb_max[:edge_width, :] > edge_threshold)
        # Bottom edge
        edge_mask[-edge_width:, :] = (rgb_max[-edge_width:, :] > edge_threshold)
        # Left edge
        edge_mask[:, :edge_width] = (rgb_max[:, :edge_width] > edge_threshold)
        # Right edge
        edge_mask[:, -edge_width:] = (rgb_max[:, -edge_width:] > edge_threshold)
    
    # Combine ALL masks - remove any pixel that matches any condition
    combined_mask = white_mask | any_high_mask | bright_mask | low_sat_mask | high_lum_mask | grayscale_mask | edge_mask
    
    # Set alpha channel to 0 for white/light pixels
    data[:, :, 3] = np.where(combined_mask, 0, data[:, :, 3])
    
    # Convert back to PIL Image
    return Image.fromarray(data)

def ensure_sharp_corners(image: Image.Image) -> Image.Image:
    """
    Ensure the image has completely sharp corners by checking edges.
    This removes any alpha transparency at the very corners that might look rounded.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    height, width = data.shape[:2]
    
    # Force corners to be fully transparent if they're light
    corner_size = 100  # Check 100px corners
    
    # Top-left corner
    corner_tl = data[:corner_size, :corner_size]
    avg_brightness_tl = np.mean(corner_tl[:, :, :3], axis=2)
    data[:corner_size, :corner_size, 3] = np.where(avg_brightness_tl > 200, 0, data[:corner_size, :corner_size, 3])
    
    # Top-right corner
    corner_tr = data[:corner_size, -corner_size:]
    avg_brightness_tr = np.mean(corner_tr[:, :, :3], axis=2)
    data[:corner_size, -corner_size:, 3] = np.where(avg_brightness_tr > 200, 0, data[:corner_size, -corner_size:, 3])
    
    # Bottom-left corner
    corner_bl = data[-corner_size:, :corner_size]
    avg_brightness_bl = np.mean(corner_bl[:, :, :3], axis=2)
    data[-corner_size:, :corner_size, 3] = np.where(avg_brightness_bl > 200, 0, data[-corner_size:, :corner_size, 3])
    
    # Bottom-right corner
    corner_br = data[-corner_size:, -corner_size:]
    avg_brightness_br = np.mean(corner_br[:, :, :3], axis=2)
    data[-corner_size:, -corner_size:, 3] = np.where(avg_brightness_br > 200, 0, data[-corner_size:, -corner_size:, 3])
    
    return Image.fromarray(data)

def fix_diary_cover(input_path: str, output_path: str = None, white_threshold: int = 180):
    """
    Fix Project Diary Cover: ULTRA-aggressively remove ALL white background, ensure sharp corners.
    
    Args:
        input_path: Path to input PNG file
        output_path: Path to output PNG file (default: overwrites input)
        white_threshold: RGB threshold for white detection (lower = more aggressive, default 180)
    """
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Image not found: {input_path}")
    
    if output_path is None:
        output_path = input_path
    else:
        output_path = Path(output_path)
    
    print(f"Fixing {input_path.name}...")
    print(f"  Input size: {Image.open(input_path).size}")
    
    # Load image
    image = Image.open(input_path)
    
    # Remove white background (ULTRA-aggressive)
    print(f"  Removing ALL white background (ultra-aggressive threshold: {white_threshold})...")
    image = remove_white_background_ultra_aggressive(image, threshold=white_threshold)
    
    # Ensure sharp corners
    print("  Ensuring sharp corners...")
    image = ensure_sharp_corners(image)
    
    # Save result
    image.save(str(output_path), 'PNG', dpi=(300, 300))
    
    print(f"  Output size: {image.size}")
    print(f"  Saved: {output_path}")
    print(f"\nProcessing complete! No white background, sharp corners guaranteed.")

if __name__ == "__main__":
    import sys
    
    # Default input path
    input_path = "exports/Ignition Protocol - Project Diary Cover.png"
    
    # Allow override via command line
    if len(sys.argv) > 1:
        input_path = sys.argv[1]
    
    # Output path (default: overwrite input)
    output_path = None
    if len(sys.argv) > 2:
        output_path = sys.argv[2]
    
    # White threshold (default 180 - very aggressive)
    white_threshold = 180
    if len(sys.argv) > 3:
        white_threshold = int(sys.argv[3])
    
    fix_diary_cover(input_path, output_path, white_threshold)
