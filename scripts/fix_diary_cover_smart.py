"""
Smart fix for Project Diary Cover: remove white background while preserving ALL text.
"""
from pathlib import Path
from PIL import Image
import numpy as np
from scipy import ndimage

def remove_white_background_smart(image: Image.Image, threshold: int = 240) -> Image.Image:
    """
    Smart white background removal that preserves text and dark content.
    Only removes pixels that are clearly white background, not text.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    original_alpha = data[:, :, 3].copy()
    
    # Calculate brightness metrics
    rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    avg_brightness = (data[:, :, 0].astype(np.float32) + 
                     data[:, :, 1].astype(np.float32) + 
                     data[:, :, 2].astype(np.float32)) / 3.0
    saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
    
    # Strategy 1: Remove pure white pixels (very high threshold to preserve text)
    # Only remove if ALL channels are very high (pure white background)
    pure_white = (data[:, :, 0] > threshold) & \
                 (data[:, :, 1] > threshold) & \
                 (data[:, :, 2] > threshold)
    
    # Strategy 2: Remove very light gray/white with low saturation (background, not text)
    # Text usually has some color, so low saturation + high brightness = background
    low_sat_white = (saturation < 20) & (avg_brightness > threshold)
    
    # Strategy 3: Remove pixels at edges that are very light (edge cleanup)
    edge_threshold = 250
    edge_mask = np.zeros_like(pure_white, dtype=bool)
    edge_width = 50
    if data.shape[0] > edge_width * 2 and data.shape[1] > edge_width * 2:
        # Only remove if it's very bright at edges
        edge_mask[:edge_width, :] = (rgb_max[:edge_width, :] > edge_threshold) & (saturation[:edge_width, :] < 15)
        edge_mask[-edge_width:, :] = (rgb_max[-edge_width:, :] > edge_threshold) & (saturation[-edge_width:, :] < 15)
        edge_mask[:, :edge_width] = (rgb_max[:, :edge_width] > edge_threshold) & (saturation[:, :edge_width] < 15)
        edge_mask[:, -edge_width:] = (rgb_max[:, -edge_width:] > edge_threshold) & (saturation[:, -edge_width:] < 15)
    
    # Strategy 4: Use morphological operations to find large white regions (background)
    # This helps distinguish between small text pixels and large background areas
    white_region = (rgb_max > threshold) & (saturation < 25)
    
    # Dilate to find large connected white regions (background)
    # Text will be small, background will be large
    from scipy.ndimage import binary_dilation, binary_erosion
    kernel_size = 20
    kernel = np.ones((kernel_size, kernel_size), dtype=bool)
    large_white_regions = binary_dilation(white_region, structure=kernel, iterations=3)
    large_white_regions = binary_erosion(large_white_regions, structure=kernel, iterations=3)
    
    # Only remove if it's in a large white region (background) AND very bright
    background_mask = large_white_regions & (rgb_max > threshold) & (saturation < 30)
    
    # Combine: remove pure white, low-sat white, edge white, and large background regions
    # BUT preserve anything that's not in these categories (text)
    combined_mask = pure_white | (low_sat_white & (avg_brightness > 245)) | edge_mask | background_mask
    
    # Preserve pixels that are clearly content (dark or colored)
    # Don't remove if it's dark (likely text/content)
    dark_content = avg_brightness < 200
    colored_content = saturation > 30
    
    # Don't remove dark or colored pixels
    preserve_mask = dark_content | colored_content
    combined_mask = combined_mask & (~preserve_mask)
    
    # Apply removal
    data[:, :, 3] = np.where(combined_mask, 0, original_alpha)
    
    return Image.fromarray(data)

def ensure_sharp_corners_only(image: Image.Image) -> Image.Image:
    """
    Only clean corners if they're clearly white background, preserve everything else.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    height, width = data.shape[:2]
    
    # Only clean corners if they're very bright and low saturation (white background)
    corner_size = 100
    rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
    
    corners = [
        (0, 0, corner_size, corner_size),  # Top-left
        (0, width - corner_size, corner_size, width),  # Top-right
        (height - corner_size, 0, height, corner_size),  # Bottom-left
        (height - corner_size, width - corner_size, height, width),  # Bottom-right
    ]
    
    for y1, x1, y2, x2 in corners:
        corner_rgb_max = rgb_max[y1:y2, x1:x2]
        corner_sat = saturation[y1:y2, x1:x2]
        # Only make transparent if very bright AND low saturation (white background)
        mask = (corner_rgb_max > 245) & (corner_sat < 20)
        data[y1:y2, x1:x2, 3] = np.where(mask, 0, data[y1:y2, x1:x2, 3])
    
    return Image.fromarray(data)

if __name__ == "__main__":
    try:
        from scipy import ndimage
    except ImportError:
        print("Installing scipy...")
        import subprocess
        subprocess.check_call(["pip", "install", "scipy"])
        from scipy import ndimage
    
    input_path = Path("exports/Ignition Protocol - Project Diary Cover.png")
    
    if not input_path.exists():
        print(f"Error: {input_path} not found")
        exit(1)
    
    print(f"Smart processing of {input_path.name}...")
    image = Image.open(input_path)
    print(f"  Input size: {image.size}")
    
    print("  Removing white background (preserving text)...")
    image = remove_white_background_smart(image, threshold=240)
    
    print("  Cleaning corners (preserving content)...")
    image = ensure_sharp_corners_only(image)
    
    output_path = input_path
    image.save(str(output_path), 'PNG', dpi=(300, 300))
    
    print(f"  Saved: {output_path}")
    print("  Complete! White background removed, text preserved, sharp corners.")


