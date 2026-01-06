"""
Convert HTML book cover to high-resolution PNG using Playwright.
"""
from pathlib import Path
import sys

def convert_html_to_png(html_path: str, output_path: str = None, dpi: int = 300):
    """
    Convert HTML to PNG at specified DPI using Playwright.
    
    Args:
        html_path: Path to HTML file
        output_path: Output path for PNG (default: same name as HTML with .png)
        dpi: Resolution in DPI (300 for print quality)
    """
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        print("Installing playwright...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "playwright"])
        subprocess.check_call([sys.executable, "-m", "playwright", "install", "chromium"])
        from playwright.sync_api import sync_playwright
    
    html_path = Path(html_path)
    if not html_path.exists():
        raise FileNotFoundError(f"HTML not found: {html_path}")
    
    if output_path is None:
        output_path = html_path.parent / f"{html_path.stem}.png"
    else:
        output_path = Path(output_path)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Converting {html_path.name} to PNG...")
    print(f"DPI: {dpi}")
    
    # Calculate pixel dimensions for A4 at specified DPI
    # A4: 210mm x 297mm = 8.27in x 11.69in
    width_px = int(8.27 * dpi)
    height_px = int(11.69 * dpi)
    
    try:
        with sync_playwright() as p:
            # Launch browser
            browser = p.chromium.launch()
            page = browser.new_page()
            
            # Set viewport to match A4 dimensions
            page.set_viewport_size({"width": width_px, "height": height_px})
            
            # Load HTML file
            file_url = f"file://{html_path.absolute()}"
            page.goto(file_url)
            
            # Wait for fonts and images to load
            page.wait_for_load_state("networkidle")
            
            # Inject CSS to ensure transparent background
            page.add_style_tag(content="body { background: transparent !important; }")
            
            # Wait a bit for styles to apply
            page.wait_for_timeout(100)
            
            # Take screenshot of just the book-cover element to preserve transparency
            book_cover = page.locator(".book-cover")
            book_cover.screenshot(
                path=str(output_path),
                type="png"
            )
            
            browser.close()
        
        print(f"Saved: {output_path}")
        print(f"   Size: {width_px} x {height_px} pixels")
        print(f"\nConversion complete!")
        
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    # Default HTML path
    html_path = "website/public/book-cover.html"
    
    # Allow override via command line
    if len(sys.argv) > 1:
        html_path = sys.argv[1]
    
    # Output path
    output_path = None
    if len(sys.argv) > 2:
        output_path = sys.argv[2]
    else:
        # Default to exports/book-cover-done.png
        output_path = "exports/book-cover-done.png"
    
    # DPI (default 300 for print quality)
    dpi = 300
    if len(sys.argv) > 3:
        dpi = int(sys.argv[3])
    
    convert_html_to_png(html_path, output_path, dpi)

