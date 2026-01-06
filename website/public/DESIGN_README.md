# Ignition Protocol - Poster & Book Cover Designs

## Quick Start

### Viewing the Designs

1. **Poster**: Open `poster.html` in your web browser
2. **Book Cover**: Open `book-cover.html` in your web browser

### Exporting for Print

#### Method 1: Browser Print to PDF
1. Open the HTML file in Chrome/Edge
2. Press `Ctrl+P` (or `Cmd+P` on Mac)
3. Select "Save as PDF"
4. Set paper size:
   - Poster: Custom 24" × 36"
   - Book Cover: Letter (8.5" × 11")
5. Set margins to "None"
6. Save the PDF

#### Method 2: Screenshot (High Quality)
1. Open HTML file in browser
2. Use browser zoom to set scale (100% = actual size)
3. Use full-page screenshot extension:
   - Chrome: "Full Page Screen Capture"
   - Firefox: "FireShot"
4. Save as PNG at 300 DPI

#### Method 3: Professional Print Service
1. Export as PDF using Method 1
2. Upload to print service (e.g., FedEx, Staples, local print shop)
3. Specify:
   - Poster: 24" × 36", matte finish
   - Book Cover: 8.5" × 11", cover stock, matte finish

## Design Flow

The visual flow is designed to guide viewers:

**Poster** (Exhibition Stand)
↓
**Website** (Interactive Exploration)
↓
**Book Cover** (Detailed Documentation)

All three share the same:
- Color palette (Flame & Darkness)
- Typography style
- Visual language
- Brand identity

## Customization

### Changing Colors
Edit the CSS color values in the HTML files:
- `#ff6b35` - Primary orange
- `#ff4500` - Red-orange
- `#ffd700` - Yellow
- `#0a0a0a` - Deep black

### Changing Text
Simply edit the HTML content:
- Title: `<h1 class="main-title">`
- Subtitle: `<p class="subtitle">`
- Author/Exhibition info: Bottom sections

### Adding QR Code
Replace the QR placeholder with an actual QR code:
1. Generate QR code linking to your website
2. Save as PNG
3. Replace `<div class="qr-placeholder">` with `<img src="qr-code.png">`

## File Structure

```
website/public/
├── poster.html          # Exhibition poster design
├── book-cover.html     # Project book cover design
├── design-guide.md    # Detailed design documentation
└── DESIGN_README.md   # This file
```

## Tips for Best Results

1. **Print Quality**: Always use 300 DPI for print
2. **Color Matching**: Use CMYK color profile for professional printing
3. **Paper Choice**: Matte finish works best for dark designs
4. **Testing**: Print a small test first to check colors
5. **Backup**: Keep digital copies of both HTML and exported PDFs

## Support

For questions or modifications, refer to `design-guide.md` for detailed specifications.














