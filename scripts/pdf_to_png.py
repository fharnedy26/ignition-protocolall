"""
Convert PDF book cover to high-resolution PNG using PyMuPDF.
"""
from pathlib import Path
try:
    import fitz  # PyMuPDF
except ImportError:
    print("Installing PyMuPDF...")
    import subprocess
    subprocess.check_call(["pip", "install", "PyMuPDF"])
    import fitz

def convert_pdf_to_png(pdf_path: str, output_dir: str = "exports", dpi: int = 300):
    """
    Convert PDF to PNG at specified DPI using PyMuPDF.
    
    Args:
        pdf_path: Path to PDF file
        output_dir: Output directory for PNG
        dpi: Resolution in DPI (300 for print quality)
    """
    pdf_path = Path(pdf_path)
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF not found: {pdf_path}")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print(f"Converting {pdf_path.name} to PNG...")
    print(f"DPI: {dpi}")
    
    try:
        # Open PDF
        doc = fitz.open(str(pdf_path))
        
        # Calculate zoom factor for desired DPI
        # Default PDF is 72 DPI, so zoom = desired_dpi / 72
        zoom = dpi / 72.0
        mat = fitz.Matrix(zoom, zoom)
        
        # Convert each page
        for page_num in range(len(doc)):
            page = doc[page_num]
            
            # Render page to pixmap
            pix = page.get_pixmap(matrix=mat)
            
            # Save as PNG
            output_path = output_dir / f"{pdf_path.stem}_page_{page_num+1}.png"
            pix.save(str(output_path))
            
            print(f"Saved: {output_path}")
            print(f"   Size: {pix.width} x {pix.height} pixels")
        
        # Get page count before closing
        page_count = len(doc)
        
        # If single page, also save as book-cover.png
        if page_count == 1:
            page = doc[0]
            pix = page.get_pixmap(matrix=mat)
            cover_path = output_dir / "book-cover.png"
            pix.save(str(cover_path))
            print(f"Also saved as: {cover_path}")
        
        doc.close()
        print(f"\nConversion complete! {page_count} page(s) converted.")
        
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    import sys
    
    # Default PDF path
    pdf_path = "Ignition Protocol - Book Cover.pdf"
    
    # Allow override via command line
    if len(sys.argv) > 1:
        pdf_path = sys.argv[1]
    
    # Output directory
    output_dir = "exports"
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    
    # DPI (default 300 for print quality)
    dpi = 300
    if len(sys.argv) > 3:
        dpi = int(sys.argv[3])
    
    convert_pdf_to_png(pdf_path, output_dir, dpi)
