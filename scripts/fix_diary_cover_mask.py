"""
Fix Project Diary Cover: detect the cover boundaries and remove everything outside.
"""
from pathlib import Path
from PIL import Image
import numpy as np

def detect_cover_boundaries(image: Image.Image, border_padding: int = 20) -> tuple:
    """
    Detect the boundaries of the actual cover content.
    Returns (left, top, right, bottom) coordinates with padding.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    data = np.array(image)
    height, width = data.shape[:2]
    
    # Calculate brightness for each pixel
    brightness = (data[:, :, 0].astype(np.float32) + 
                  data[:, :, 1].astype(np.float32) + 
                  data[:, :, 2].astype(np.float32)) / 3.0
    
    # The cover is dark, background is white
    # Find where content starts (non-white pixels)
    # Look for pixels that are NOT white (brightness < 240)
    
    # Find top boundary: scan from top, find first row with dark content
    top = 0
    for y in range(height):
        row_brightness = brightness[y, :]
        # If more than 1% of row is dark (not white), we found the top
        dark_pixels = np.sum(row_brightness < 240)
        if dark_pixels > width * 0.01:
            top = max(0, y - border_padding)
            break
    
    # Find bottom boundary: scan from bottom
    bottom = height
    for y in range(height - 1, -1, -1):
        row_brightness = brightness[y, :]
        dark_pixels = np.sum(row_brightness < 240)
        if dark_pixels > width * 0.01:
            bottom = min(height, y + 1 + border_padding)
            break
    
    # Find left boundary: scan from left
    left = 0
    for x in range(width):
        col_brightness = brightness[:, x]
        dark_pixels = np.sum(col_brightness < 240)
        if dark_pixels > height * 0.01:
            left = max(0, x - border_padding)
            break
    
    # Find right boundary: scan from right
    right = width
    for x in range(width - 1, -1, -1):
        col_brightness = brightness[:, x]
        dark_pixels = np.sum(col_brightness < 240)
        if dark_pixels > height * 0.01:
            right = min(width, x + 1 + border_padding)
            break
    
    return (left, top, right, bottom)

def remove_outside_cover(image: Image.Image, border_padding: int = 20, target_size: tuple = None) -> Image.Image:
    """
    Detect cover boundaries and remove everything outside the cover + border.
    If target_size is specified, crop to exactly that size centered on the detected cover.
    """
    if image.mode != 'RGBA':
        image = image.convert('RGBA')
    
    # Detect cover boundaries
    left, top, right, bottom = detect_cover_boundaries(image, 0)  # Detect without padding first
    
    print(f"  Detected cover boundaries: left={left}, top={top}, right={right}, bottom={bottom}")
    print(f"  Cover size: {right-left} x {bottom-top} pixels")
    
    data = np.array(image)
    height, width = data.shape[:2]
    
    if target_size:
        # Crop to exact target size, with more cropping from bottom
        target_width, target_height = target_size
        cover_center_x = (left + right) // 2
        
        # For vertical positioning: crop more from bottom, less from top
        # Shift crop down to remove more white from bottom
        crop_left = max(0, cover_center_x - target_width // 2)
        # Start crop slightly above detected top, but shift down more
        # This removes more from bottom while keeping same height
        crop_top = max(0, top - 20)  # Small margin from top
        crop_right = min(width, crop_left + target_width)
        crop_bottom = min(height, crop_top + target_height)
        
        # If there's room, shift down to crop more from bottom
        if crop_bottom < height - 50:
            # Shift down by 50px to remove more bottom white
            crop_top = min(crop_top + 50, height - target_height)
            crop_bottom = crop_top + target_height
        
        # Adjust if we hit edges
        if crop_right - crop_left < target_width:
            crop_left = max(0, crop_right - target_width)
        if crop_bottom - crop_top < target_height:
            # If we need more height, take it from top
            crop_top = max(0, crop_bottom - target_height)
        
        print(f"  Cropping to exact size: {target_width} x {target_height}")
        print(f"  Crop box: left={crop_left}, top={crop_top}, right={crop_right}, bottom={crop_bottom}")
        
        # Crop the image
        cropped_data = data[crop_top:crop_bottom, crop_left:crop_right]
        
        # Create new image with exact target size
        result = Image.new('RGBA', (target_width, target_height), (0, 0, 0, 0))
        result_data = np.array(result)
        
        # Copy cropped data to result
        copy_height = min(cropped_data.shape[0], target_height)
        copy_width = min(cropped_data.shape[1], target_width)
        result_data[:copy_height, :copy_width] = cropped_data[:copy_height, :copy_width]
        
        return Image.fromarray(result_data)
    else:
        # Original behaviour: mask based on detected boundaries with padding
        left, top, right, bottom = detect_cover_boundaries(image, border_padding)
        
        # Create mask (1 = keep, 0 = remove)
        mask = np.zeros((height, width), dtype=bool)
        mask[top:bottom, left:right] = True
        
        # Apply mask to alpha channel
        data[:, :, 3] = np.where(mask, data[:, :, 3], 0)
        
        return Image.fromarray(data)

if __name__ == "__main__":
    import sys
    import io
    
    # First, reload from PDF to get original
    print("Loading fresh image from PDF...")
    try:
        import fitz  # PyMuPDF
        doc = fitz.open("Ignition Protocol - Project Diary Cover.pdf")
        zoom = 300 / 72.0
        mat = fitz.Matrix(zoom, zoom)
        page = doc[0]
        pix = page.get_pixmap(matrix=mat)
        img_data = pix.tobytes("png")
        image = Image.open(io.BytesIO(img_data))
        doc.close()
        print("  Loaded from PDF")
    except Exception as e:
        print(f"  Error loading PDF: {e}")
        input_path = Path("exports/Ignition Protocol - Project Diary Cover.png")
        if not input_path.exists():
            print(f"Error: {input_path} not found")
            exit(1)
        image = Image.open(input_path)
        print("  Loaded from existing PNG")
    
    print(f"  Image size: {image.size}")
    print()
    
    # Target size (width x height) or border padding
    target_size = None
    border_padding = 20
    
    if len(sys.argv) > 1:
        # Check if it's a size specification (e.g., "3500x4900")
        arg = sys.argv[1]
        if 'x' in arg.lower():
            parts = arg.lower().split('x')
            if len(parts) == 2:
                target_size = (int(parts[0]), int(parts[1]))
                print(f"Cropping to exact size: {target_size[0]} x {target_size[1]} pixels")
        else:
            border_padding = int(arg)
    
    if target_size:
        print(f"Cropping to exact size: {target_size[0]} x {target_size[1]} pixels...")
    else:
        print(f"Removing everything outside cover (border padding: {border_padding}px)...")
    
    image = remove_outside_cover(image, border_padding, target_size)
    
    output_path = Path("exports/Ignition Protocol - Project Diary Cover.png")
    image.save(str(output_path), 'PNG', dpi=(300, 300))
    
    print(f"  Saved: {output_path}")
    print("  Complete! Everything outside cover removed, all original colors preserved.")

