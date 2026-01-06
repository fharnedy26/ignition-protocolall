"""
Generate A0 poster background as high-resolution PNG or PDF.

A0 size (landscape): 1189mm × 841mm
At 300 DPI: 14043 × 9933 pixels
At 150 DPI: 7022 × 4967 pixels (faster, still good quality)
"""

from PIL import Image, ImageDraw, ImageFilter
import math

def create_a0_background(output_path="poster_background_a0.png", dpi=300, format="PNG"):
    """
    Create an A0-sized poster background.
    
    Args:
        output_path: Output file path
        dpi: Resolution (300 for print, 150 for screen)
        format: "PNG" or "PDF"
    """
    # A0 dimensions in mm (landscape orientation)
    width_mm = 1189
    height_mm = 841
    
    # Convert to pixels
    mm_to_inch = 1 / 25.4
    width_px = int(width_mm * mm_to_inch * dpi)
    height_px = int(height_mm * mm_to_inch * dpi)
    
    print(f"Creating A0 poster background: {width_px} × {height_px} pixels ({dpi} DPI)")
    
    # Create image with Flame & Darkness theme background
    img = Image.new('RGB', (width_px, height_px), color='#0a0a0a')
    draw = ImageDraw.Draw(img)
    
    # Draw gradient background: Deep Red (#1a0000) to Black (#0a0a0a) - horizontal for landscape
    num_steps = 100
    for i in range(num_steps):
        ratio = i / num_steps
        # Interpolate from #1a0000 (26, 0, 0) to #0a0a0a (10, 10, 10)
        r = int(26 - ratio * 16)  # 26 -> 10
        g = int(ratio * 10)        # 0 -> 10
        b = int(ratio * 10)        # 0 -> 10
        
        x_start = int(width_px * i / num_steps)
        x_end = int(width_px * (i + 1) / num_steps)
        draw.rectangle([(x_start, 0), (x_end, height_px)], fill=(r, g, b))
    
    # Draw flame glow effects (radial gradients approximated)
    from PIL import ImageFilter
    
    # Create flame glow layers (adjusted for landscape)
    glow1 = Image.new('RGBA', (width_px, height_px), (0, 0, 0, 0))
    glow_draw1 = ImageDraw.Draw(glow1)
    # Left side glow (orange)
    glow_draw1.ellipse([-width_px*0.2, height_px*0.1, width_px*0.5, height_px*0.9], 
                      fill=(255, 107, 53, 38))  # rgba(255, 107, 53, 0.15)
    glow1 = glow1.filter(ImageFilter.GaussianBlur(radius=width_px * 0.08))
    
    glow2 = Image.new('RGBA', (width_px, height_px), (0, 0, 0, 0))
    glow_draw2 = ImageDraw.Draw(glow2)
    # Right side glow (red-orange)
    glow_draw2.ellipse([width_px*0.7, height_px*0.1, width_px*1.4, height_px*0.9], 
                      fill=(255, 69, 0, 31))  # rgba(255, 69, 0, 0.12)
    glow2 = glow2.filter(ImageFilter.GaussianBlur(radius=width_px * 0.08))
    
    glow3 = Image.new('RGBA', (width_px, height_px), (0, 0, 0, 0))
    glow_draw3 = ImageDraw.Draw(glow3)
    # Center glow (yellow)
    glow_draw3.ellipse([width_px*0.25, height_px*0.1, width_px*0.75, height_px*0.9], 
                      fill=(255, 215, 0, 20))  # rgba(255, 215, 0, 0.08)
    glow3 = glow3.filter(ImageFilter.GaussianBlur(radius=width_px * 0.12))
    
    # Composite glow effects
    img = Image.alpha_composite(img.convert('RGBA'), glow1).convert('RGB')
    img = Image.alpha_composite(img.convert('RGBA'), glow2).convert('RGB')
    img = Image.alpha_composite(img.convert('RGBA'), glow3).convert('RGB')
    draw = ImageDraw.Draw(img)
    
    # Draw subtle grid pattern
    grid_spacing = int(50 * mm_to_inch * dpi)  # 50mm grid
    grid_color = (255, 255, 255, 8)  # Very subtle white (0.03 opacity)
    
    for x in range(0, width_px, grid_spacing):
        draw.line([(x, 0), (x, height_px)], fill=grid_color, width=1)
    for y in range(0, height_px, grid_spacing):
        draw.line([(0, y), (width_px, y)], fill=grid_color, width=1)
    
    # Draw decorative molecular structure elements with flame colors (scaled)
    node_radius = int(8 * mm_to_inch * dpi)
    line_width = int(1.5 * mm_to_inch * dpi)
    
    # Flame orange color for molecules
    molecule_color = (255, 107, 53, 64)  # rgba(255, 107, 53, 0.25)
    bond_color = (255, 107, 53, 102)    # rgba(255, 107, 53, 0.4)
    
    # Top left cluster
    nodes = [
        (int(150 * mm_to_inch * dpi), int(200 * mm_to_inch * dpi)),
        (int(250 * mm_to_inch * dpi), int(180 * mm_to_inch * dpi)),
        (int(200 * mm_to_inch * dpi), int(250 * mm_to_inch * dpi)),
        (int(300 * mm_to_inch * dpi), int(220 * mm_to_inch * dpi)),
    ]
    
    # Draw connections
    connections = [(0, 1), (2, 1), (2, 3)]
    for i, j in connections:
        x1, y1 = nodes[i]
        x2, y2 = nodes[j]
        # Create bond line with gradient effect
        bond_img = Image.new('RGBA', (width_px, height_px), (0, 0, 0, 0))
        bond_draw = ImageDraw.Draw(bond_img)
        bond_draw.line([(x1, y1), (x2, y2)], fill=bond_color, width=line_width)
        img = Image.alpha_composite(img.convert('RGBA'), bond_img).convert('RGB')
        draw = ImageDraw.Draw(img)
    
    # Draw nodes with radial gradient effect
    for x, y in nodes:
        # Create node with gradient
        node_img = Image.new('RGBA', (node_radius*4, node_radius*4), (0, 0, 0, 0))
        node_draw = ImageDraw.Draw(node_img)
        node_draw.ellipse([0, 0, node_radius*4, node_radius*4], 
                        fill=(255, 107, 53, 153))  # rgba(255, 107, 53, 0.6)
        node_draw.ellipse([node_radius*0.5, node_radius*0.5, node_radius*3.5, node_radius*3.5], 
                        fill=(255, 69, 0, 51))     # rgba(255, 69, 0, 0.2)
        node_img = node_img.filter(ImageFilter.GaussianBlur(radius=node_radius*0.3))
        img.paste(node_img, (x-node_radius*2, y-node_radius*2), node_img)
        draw = ImageDraw.Draw(img)
    
    # Mirror to other corners
    # Top right
    for x, y in nodes:
        mirror_x = width_px - x
        node_img = Image.new('RGBA', (node_radius*4, node_radius*4), (0, 0, 0, 0))
        node_draw = ImageDraw.Draw(node_img)
        node_draw.ellipse([0, 0, node_radius*4, node_radius*4], 
                        fill=(255, 107, 53, 153))
        node_draw.ellipse([node_radius*0.5, node_radius*0.5, node_radius*3.5, node_radius*3.5], 
                        fill=(255, 69, 0, 51))
        node_img = node_img.filter(ImageFilter.GaussianBlur(radius=node_radius*0.3))
        img.paste(node_img, (mirror_x-node_radius*2, y-node_radius*2), node_img)
    
    # Bottom left
    for x, y in nodes:
        mirror_y = height_px - y
        node_img = Image.new('RGBA', (node_radius*4, node_radius*4), (0, 0, 0, 0))
        node_draw = ImageDraw.Draw(node_img)
        node_draw.ellipse([0, 0, node_radius*4, node_radius*4], 
                        fill=(255, 107, 53, 153))
        node_draw.ellipse([node_radius*0.5, node_radius*0.5, node_radius*3.5, node_radius*3.5], 
                        fill=(255, 69, 0, 51))
        node_img = node_img.filter(ImageFilter.GaussianBlur(radius=node_radius*0.3))
        img.paste(node_img, (x-node_radius*2, mirror_y-node_radius*2), node_img)
    
    # Bottom right
    for x, y in nodes:
        mirror_x = width_px - x
        mirror_y = height_px - y
        node_img = Image.new('RGBA', (node_radius*4, node_radius*4), (0, 0, 0, 0))
        node_draw = ImageDraw.Draw(node_img)
        node_draw.ellipse([0, 0, node_radius*4, node_radius*4], 
                        fill=(255, 107, 53, 153))
        node_draw.ellipse([node_radius*0.5, node_radius*0.5, node_radius*3.5, node_radius*3.5], 
                        fill=(255, 69, 0, 51))
        node_img = node_img.filter(ImageFilter.GaussianBlur(radius=node_radius*0.3))
        img.paste(node_img, (mirror_x-node_radius*2, mirror_y-node_radius*2), node_img)
    
    draw = ImageDraw.Draw(img)
    
    # Draw subtle border with flame accent
    border_width = int(2 * mm_to_inch * dpi)
    border_margin = int(20 * mm_to_inch * dpi)
    border_img = Image.new('RGBA', (width_px, height_px), (0, 0, 0, 0))
    border_draw = ImageDraw.Draw(border_img)
    border_draw.rectangle([border_margin, border_margin, width_px-border_margin, height_px-border_margin],
                  outline=(255, 107, 53, 51), width=border_width)  # rgba(255, 107, 53, 0.2)
    img = Image.alpha_composite(img.convert('RGBA'), border_img).convert('RGB')
    draw = ImageDraw.Draw(img)
    
    # Save image
    if format.upper() == "PDF":
        # Convert to PDF
        rgb_img = img.convert('RGB')
        rgb_img.save(output_path.replace('.png', '.pdf'), 'PDF', resolution=dpi)
        print(f"Saved PDF: {output_path.replace('.png', '.pdf')}")
    else:
        img.save(output_path, 'PNG', dpi=(dpi, dpi))
        print(f"Saved PNG: {output_path}")
    
    print(f"File size: {width_px} × {height_px} pixels")
    return img

if __name__ == "__main__":
    import sys
    
    dpi = 300
    if len(sys.argv) > 1:
        dpi = int(sys.argv[1])
    
    format_type = "PNG"
    if len(sys.argv) > 2:
        format_type = sys.argv[2].upper()
    
    create_a0_background(dpi=dpi, format=format_type)
