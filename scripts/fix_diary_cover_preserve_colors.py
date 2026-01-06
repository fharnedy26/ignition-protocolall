"""
Fix Project Diary Cover: remove ONLY pure white background, preserve ALL original colors including gray text.
"""
from pathlib import Path
from PIL import Image
import numpy as np

def remove_only_pure_white(image: Image.Image) -> Image.Image:
    """
    Remove ONLY pure white background pixels. Preserve everything else including gray text.
    Very conservative approach - only removes pixels that are clearly white background.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    original_alpha = data[:, :, 3].copy()
    
    # Calculate metrics
    rgb_max = np.maximum(np.maximum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    rgb_min = np.minimum(np.minimum(data[:, :, 0], data[:, :, 1]), data[:, :, 2])
    avg_brightness = (data[:, :, 0].astype(np.float32) + 
                     data[:, :, 1].astype(np.float32) + 
                     data[:, :, 2].astype(np.float32)) / 3.0
    saturation = rgb_max.astype(np.float32) - rgb_min.astype(np.float32)
    
    # VERY CONSERVATIVE: Only remove pixels that are:
    # 1. Very bright (all channels > 250)
    # 2. Very low saturation (almost no color, just white)
    # 3. This ensures we don't touch gray text or any colored content
    
    # Pure white background: all channels very high AND very low saturation
    pure_white = (data[:, :, 0] > 250) & \
                 (data[:, :, 1] > 250) & \
                 (data[:, :, 2] > 250) & \
                 (saturation < 10)
    
    # Also remove very light gray-white at edges (background artifacts)
    edge_threshold = 252
    edge_mask = np.zeros_like(pure_white, dtype=bool)
    edge_width = 30
    if data.shape[0] > edge_width * 2 and data.shape[1] > edge_width * 2:
        # Only at very edges, and only if it's pure white
        edge_mask[:edge_width, :] = (rgb_max[:edge_width, :] > edge_threshold) & \
                                    (saturation[:edge_width, :] < 5)
        edge_mask[-edge_width:, :] = (rgb_max[-edge_width:, :] > edge_threshold) & \
                                      (saturation[-edge_width:, :] < 5)
        edge_mask[:, :edge_width] = (rgb_max[:, :edge_width] > edge_threshold) & \
                                     (saturation[:, :edge_width] < 5)
        edge_mask[:, -edge_width:] = (rgb_max[:, -edge_width:] > edge_threshold) & \
                                      (saturation[:, -edge_width:] < 5)
    
    # Combine: only remove pure white background
    combined_mask = pure_white | edge_mask
    
    # Apply removal - preserve everything else (including gray text)
    data[:, :, 3] = np.where(combined_mask, 0, original_alpha)
    
    return Image.fromarray(data)

def clean_corners_conservative(image: Image.Image) -> Image.Image:
    """
    Only clean corners if they're clearly pure white background.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    height, width = data.shape[:2]
    
    # Only remove corners if they're pure white (all channels > 250, low saturation)
    corner_size = 50
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
        corner_r = data[y1:y2, x1:x2, 0]
        corner_g = data[y1:y2, x1:x2, 1]
        corner_b = data[y1:y2, x1:x2, 2]
        # Only remove if pure white (all channels > 250 AND low saturation)
        mask = (corner_r > 250) & (corner_g > 250) & (corner_b > 250) & (corner_sat < 5)
        data[y1:y2, x1:x2, 3] = np.where(mask, 0, data[y1:y2, x1:x2, 3])
    
    return Image.fromarray(data)

if __name__ == "__main__":
    import sys
    import io
    sys.path.insert(0, str(Path(__file__).parent.parent))
    
    # First, reprocess from PDF to get fresh image with original colors
    print("Reprocessing from PDF to get original colors...")
    try:
        import fitz  # PyMuPDF
        doc = fitz.open("Ignition Protocol - Project Diary Cover.pdf")
        zoom = 300 / 72.0
        mat = fitz.Matrix(zoom, zoom)
        page = doc[0]
        pix = page.get_pixmap(matrix=mat)
        img_data = pix.tobytes("png")
        fresh_image = Image.open(io.BytesIO(img_data))
        fresh_image.save("exports/Ignition Protocol - Project Diary Cover.png", 'PNG', dpi=(300, 300))
        doc.close()
        print("  Fresh image loaded from PDF")
    except Exception as e:
        print(f"  Warning: Could not reload from PDF: {e}")
        print("  Using existing PNG file")
    
    print()
    
    input_path = Path("exports/Ignition Protocol - Project Diary Cover.png")
    
    if not input_path.exists():
        print(f"Error: {input_path} not found")
        exit(1)
    
    print(f"Conservative processing of {input_path.name}...")
    image = Image.open(input_path)
    print(f"  Input size: {image.size}")
    
    print("  Removing ONLY pure white background (preserving ALL colors)...")
    image = remove_only_pure_white(image)
    
    print("  Cleaning corners (preserving content)...")
    image = clean_corners_conservative(image)
    
    output_path = input_path
    image.save(str(output_path), 'PNG', dpi=(300, 300))
    
    print(f"  Saved: {output_path}")
    print("  Complete! Only pure white removed, ALL original colors preserved.")

