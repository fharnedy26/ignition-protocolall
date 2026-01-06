"""
EXTREME fix for Project Diary Cover: remove EVERY white pixel, ensure perfect sharp corners.
"""
from pathlib import Path
from PIL import Image
import numpy as np

def remove_all_white_extreme(image: Image.Image) -> Image.Image:
    """
    EXTREMELY aggressively remove ALL white/light pixels.
    Uses multiple passes and very low thresholds.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    
    # Multiple passes with different strategies
    for pass_num in range(3):
        # Pass 1: Very aggressive threshold (120)
        # Pass 2: Medium aggressive (150) 
        # Pass 3: Light aggressive (180)
        thresholds = [120, 150, 180]
        threshold = thresholds[min(pass_num, 2)]
        
        # Calculate brightness metrics
        rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
        rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
        avg_brightness = (data[:, :, 0].astype(np.float32) + 
                         data[:, :, 1].astype(np.float32) + 
                         data[:, :, 2].astype(np.float32)) / 3.0
        saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
        luminance = (0.299 * data[:, :, 0].astype(np.float32) + 
                     0.587 * data[:, :, 1].astype(np.float32) + 
                     0.114 * data[:, :, 2].astype(np.float32))
        
        # Multiple removal strategies
        mask1 = (data[:, :, 0] > threshold) & (data[:, :, 1] > threshold) & (data[:, :, 2] > threshold)
        mask2 = avg_brightness > threshold
        mask3 = luminance > threshold
        mask4 = (saturation < 60) & (rgb_max > threshold)
        mask5 = rgb_max > threshold + 20
        
        # Edge-specific removal (more aggressive at edges)
        edge_threshold = threshold + 30
        edge_mask = np.zeros_like(mask1, dtype=bool)
        edge_width = 100
        if data.shape[0] > edge_width * 2 and data.shape[1] > edge_width * 2:
            edge_mask[:edge_width, :] = rgb_max[:edge_width, :] > edge_threshold
            edge_mask[-edge_width:, :] = rgb_max[-edge_width:, :] > edge_threshold
            edge_mask[:, :edge_width] = rgb_max[:, :edge_width] > edge_threshold
            edge_mask[:, -edge_width:] = rgb_max[:, -edge_width:] > edge_threshold
        
        # Corner-specific removal (very aggressive)
        corner_size = 150
        corner_mask = np.zeros_like(mask1, dtype=bool)
        corner_threshold = threshold + 50
        if data.shape[0] > corner_size * 2 and data.shape[1] > corner_size * 2:
            # Top-left
            corner_mask[:corner_size, :corner_size] = rgb_max[:corner_size, :corner_size] > corner_threshold
            # Top-right
            corner_mask[:corner_size, -corner_size:] = rgb_max[:corner_size, -corner_size:] > corner_threshold
            # Bottom-left
            corner_mask[-corner_size:, :corner_size] = rgb_max[-corner_size:, :corner_size] > corner_threshold
            # Bottom-right
            corner_mask[-corner_size:, -corner_size:] = rgb_max[-corner_size:, -corner_size:] > corner_threshold
        
        # Combine all masks
        combined = mask1 | mask2 | mask3 | mask4 | mask5 | edge_mask | corner_mask
        
        # Remove pixels
        data[:, :, 3] = np.where(combined, 0, data[:, :, 3])
    
    return Image.fromarray(data)

def force_sharp_corners(image: Image.Image) -> Image.Image:
    """
    Force completely sharp corners by making corner regions transparent if they're light.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    height, width = data.shape[:2]
    
    # Very aggressive corner cleaning
    corner_size = 200
    
    # Check each corner and make transparent if average brightness is high
    corners = [
        (0, 0, corner_size, corner_size),  # Top-left
        (0, width - corner_size, corner_size, width),  # Top-right
        (height - corner_size, 0, height, corner_size),  # Bottom-left
        (height - corner_size, width - corner_size, height, width),  # Bottom-right
    ]
    
    for y1, x1, y2, x2 in corners:
        corner = data[y1:y2, x1:x2]
        if corner.size > 0:
            avg_brightness = np.mean(corner[:, :, :3], axis=2)
            # Make transparent if brightness > 150
            mask = avg_brightness > 150
            data[y1:y2, x1:x2, 3] = np.where(mask, 0, data[y1:y2, x1:x2, 3])
    
    return Image.fromarray(data)

if __name__ == "__main__":
    input_path = Path("exports/Ignition Protocol - Project Diary Cover.png")
    
    if not input_path.exists():
        print(f"Error: {input_path} not found")
        exit(1)
    
    print(f"EXTREME processing of {input_path.name}...")
    image = Image.open(input_path)
    print(f"  Input size: {image.size}")
    
    print("  Removing ALL white pixels (extreme mode)...")
    image = remove_all_white_extreme(image)
    
    print("  Forcing sharp corners...")
    image = force_sharp_corners(image)
    
    output_path = input_path
    image.save(str(output_path), 'PNG', dpi=(300, 300))
    
    print(f"  Saved: {output_path}")
    print("  Complete! All white removed, corners are sharp.")


