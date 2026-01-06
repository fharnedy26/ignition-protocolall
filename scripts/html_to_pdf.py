"""
Convert HTML book cover to PDF using Playwright.
"""
from pathlib import Path
import sys

def convert_html_to_pdf(html_path: str, output_path: str = None):
    """
    Convert HTML to PDF using Playwright.
    
    Args:
        html_path: Path to HTML file
        output_path: Output path for PDF (default: same name as HTML with .pdf)
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
        output_path = html_path.parent / f"{html_path.stem}.pdf"
    else:
        output_path = Path(output_path)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Converting {html_path.name} to PDF...")
    
    try:
        with sync_playwright() as p:
            # Launch browser
            browser = p.chromium.launch()
            page = browser.new_page()
            
            # Load HTML file
            file_url = f"file://{html_path.absolute()}"
            page.goto(file_url)
            
            # Wait for fonts and images to load
            page.wait_for_load_state("networkidle")
            
            # Inject CSS to ensure transparent background is handled
            page.add_style_tag(content="body { background: transparent !important; }")
            
            # Wait a bit for styles to apply
            page.wait_for_timeout(100)
            
            # Generate PDF with A4 size
            page.pdf(
                path=str(output_path),
                format="A4",
                print_background=True,
                margin={
                    "top": "0",
                    "right": "0",
                    "bottom": "0",
                    "left": "0"
                }
            )
            
            browser.close()
        
        print(f"Saved: {output_path}")
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
        # Default to exports/book-cover-done.pdf
        output_path = "exports/book-cover-done.pdf"
    
    convert_html_to_pdf(html_path, output_path)










