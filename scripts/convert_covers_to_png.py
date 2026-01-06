"""
Convert PDF covers to PNG with white background removed.
Preserves dark cover and text while removing white background.
"""
from pathlib import Path
import sys
import io

# Add parent directory to path to import existing modules
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import fitz  # PyMuPDF
except ImportError:
    print("Installing PyMuPDF...")
    import subprocess
    subprocess.check_call(["pip", "install", "PyMuPDF"])
    import fitz

from PIL import Image
import numpy as np

def remove_white_background(image: Image.Image, threshold: int = 240) -> Image.Image:
    """
    Remove white/light background by making white and near-white pixels transparent.
    Preserves dark cover and text.
    
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

def convert_pdf_to_png_no_white(pdf_path: str, output_path: str = None, dpi: int = 300):
    """
    Convert PDF to PNG and remove white background.
    
    Args:
        pdf_path: Path to PDF file
        output_path: Path to output PNG file (default: same name with .png extension)
        dpi: Resolution in DPI (300 for print quality)
    """
    pdf_path = Path(pdf_path)
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF not found: {pdf_path}")
    
    if output_path is None:
        output_path = pdf_path.parent / f"{pdf_path.stem}.png"
    else:
        output_path = Path(output_path)
    
    print(f"Converting {pdf_path.name} to PNG...")
    print(f"DPI: {dpi}")
    
    try:
        # Open PDF
        doc = fitz.open(str(pdf_path))
        
        if len(doc) == 0:
            raise ValueError("PDF has no pages")
        
        # Calculate zoom factor for desired DPI
        zoom = dpi / 72.0
        mat = fitz.Matrix(zoom, zoom)
        
        # Convert first page (assuming covers are single page)
        page = doc[0]
        
        # Render page to pixmap
        pix = page.get_pixmap(matrix=mat)
        
        # Convert to PIL Image
        img_data = pix.tobytes("png")
        image = Image.open(io.BytesIO(img_data))
        
        print(f"  Original size: {image.size}")
        
        # Remove white background
        print("  Removing white background...")
        image = remove_white_background(image, threshold=240)
        
        # Save as PNG
        image.save(str(output_path), 'PNG', dpi=(dpi, dpi))
        
        print(f"  Output size: {image.size}")
        print(f"  Saved: {output_path}")
        
        doc.close()
        print(f"\nConversion complete!")
        
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    # PDF files to convert
    pdf_files = [
        "Ignition Protocol - Book Cover.pdf",
        "Ignition Protocol - Project Diary Cover.pdf"
    ]
    
    # Output directory
    output_dir = Path("exports")
    output_dir.mkdir(exist_ok=True)
    
    # DPI (300 for print quality)
    dpi = 300
    
    # Process each PDF
    for pdf_file in pdf_files:
        pdf_path = Path(pdf_file)
        if not pdf_path.exists():
            print(f"Warning: {pdf_file} not found, skipping...")
            continue
        
        output_path = output_dir / f"{pdf_path.stem}.png"
        
        try:
            convert_pdf_to_png_no_white(str(pdf_path), str(output_path), dpi)
            print()
        except Exception as e:
            print(f"Error processing {pdf_file}: {e}")
            print()

