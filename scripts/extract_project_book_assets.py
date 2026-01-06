"""
Extract readable text and key visual assets from the project book PDF.

Outputs (relative to repo root):
- exports/ignition_protocol_text.txt            (full text, one page after another)
- exports/ignition_protocol_pages/page_XXX.txt  (per-page text)
- exports/front_page_render.png                (rendered full page 1)
- exports/front_page_background.*              (largest embedded image on page 1)
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import fitz  # PyMuPDF


@dataclass(frozen=True)
class ExtractedImage:
    xref: int
    width: int
    height: int
    bpc: int
    colorspace: str
    ext: str
    size_bytes: int


def _safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _extract_text(doc: "fitz.Document", out_dir: Path) -> Path:
    pages_dir = out_dir / "ignition_protocol_pages"
    _safe_mkdir(pages_dir)

    full_text_path = out_dir / "ignition_protocol_text.txt"
    with full_text_path.open("w", encoding="utf-8", newline="\n") as f_all:
        for i in range(doc.page_count):
            page = doc.load_page(i)
            txt = page.get_text("text") or ""
            # Normalize: ensure each page is clearly separated.
            header = f"\n\n===== PAGE {i+1} / {doc.page_count} =====\n\n"
            f_all.write(header)
            f_all.write(txt)

            (pages_dir / f"page_{i+1:03d}.txt").write_text(
                txt, encoding="utf-8", newline="\n"
            )

    return full_text_path


def _render_first_page(doc: "fitz.Document", out_dir: Path) -> Path:
    page0 = doc.load_page(0)
    # Render at ~200 DPI equivalent (zoom ~ 2.78 for 72dpi base).
    mat = fitz.Matrix(2.8, 2.8)
    pix = page0.get_pixmap(matrix=mat, alpha=False)
    out_path = out_dir / "front_page_render.png"
    pix.save(str(out_path))
    return out_path


def _collect_page_images(doc: "fitz.Document", page_index: int) -> List[ExtractedImage]:
    page = doc.load_page(page_index)
    imgs = []
    for img in page.get_images(full=True):
        # img tuple: (xref, smask, width, height, bpc, colorspace, alt. colorspace, name, filter, ...)
        xref = int(img[0])
        width = int(img[2])
        height = int(img[3])
        bpc = int(img[4])
        cs = str(img[5])
        try:
            d = doc.extract_image(xref)
            ext = str(d.get("ext") or "bin")
            size_bytes = len(d.get("image") or b"")
        except Exception:
            ext = "bin"
            size_bytes = 0
        imgs.append(
            ExtractedImage(
                xref=xref,
                width=width,
                height=height,
                bpc=bpc,
                colorspace=cs,
                ext=ext,
                size_bytes=size_bytes,
            )
        )
    return imgs


def _export_largest_first_page_image(doc: "fitz.Document", out_dir: Path) -> Tuple[Optional[Path], List[ExtractedImage]]:
    imgs = _collect_page_images(doc, 0)
    if not imgs:
        return None, imgs

    # Heuristic: largest pixel area tends to be background.
    imgs_sorted = sorted(imgs, key=lambda x: (x.width * x.height, x.size_bytes), reverse=True)
    best = imgs_sorted[0]

    extracted = doc.extract_image(best.xref)
    raw = extracted.get("image") or b""
    ext = extracted.get("ext") or best.ext or "bin"

    raw_path = out_dir / f"front_page_background.{ext}"
    raw_path.write_bytes(raw)

    # Also try to save as PNG for consistent viewing.
    png_path = out_dir / "front_page_background.png"
    try:
        pix = fitz.Pixmap(doc, best.xref)
        if pix.n > 4:  # CMYK or other
            pix = fitz.Pixmap(fitz.csRGB, pix)
        pix.save(str(png_path))
        return png_path, imgs_sorted
    except Exception:
        return raw_path, imgs_sorted


def main() -> int:
    pdf_path = Path("Ignition Protocol SYSTE'26 (4).pdf")
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF not found: {pdf_path}")

    out_dir = Path("exports")
    _safe_mkdir(out_dir)

    doc = fitz.open(str(pdf_path))
    try:
        text_path = _extract_text(doc, out_dir)
        render_path = _render_first_page(doc, out_dir)
        bg_path, imgs = _export_largest_first_page_image(doc, out_dir)

        print(f"[OK] Pages: {doc.page_count}")
        print(f"[OK] Full text: {text_path}")
        print(f"[OK] First page render: {render_path}")
        print(f"[OK] First page images found: {len(imgs)}")
        if bg_path:
            print(f"[OK] Exported background candidate: {bg_path}")
            best = imgs[0]
            print(
                f"[IMG] xref={best.xref} {best.width}x{best.height} bpc={best.bpc} cs={best.colorspace} ext={best.ext} bytes={best.size_bytes}"
            )
        else:
            print("[WARN] No images found on the first page.")
    finally:
        doc.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


