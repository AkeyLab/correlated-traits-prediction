"""Recolour Figure 1 (Illustrator/PDF source) into the house palette and render at 600 dpi."""
import re, sys, numpy as np, fitz
from PIL import Image
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
src, out = sys.argv[1], sys.argv[2]
BLUE, ORANGE, GREEN = (0.122, 0.467, 0.706), (1.0, 0.498, 0.055), (0.173, 0.627, 0.173)
LIGHT_BLUE = (0.776, 0.859, 0.937)
# CMYK operands as written in the stream -> RGB (keep-colours use the values MuPDF renders them to)
CMYK = {(0.711, 0.672, 0, 0): BLUE,            # hexagons E1 E2 G1 G2
        (0.172, 0.172, 0, 0): LIGHT_BLUE,       # phenotype boxes
        (0, 0.914, 1, 0): ORANGE,               # "shared" dashed arrows
        (0.426, 0.348, 0.352, 0.012): (0.6, 0.6, 0.6),
        (0.746, 0.676, 0.668, 0.898): (0.012, 0.017, 0.018),
        (0.668, 0.59, 0.578, 0.41): (0.277, 0.282, 0.284),
        (0, 0, 0, 0): (1, 1, 1), (0, 0, 0, 1): (0, 0, 0)}
doc = fitz.open(src); page = doc[0]
stream = page.read_contents().decode('latin1')
def swap(m):
    vals = tuple(float(v) for v in m.group(1).split())
    key = tuple(round(v, 3) for v in vals)
    rgb = CMYK.get(key) or tuple((1 - c) * (1 - vals[3]) for c in vals[:3])
    op = 'rg' if m.group(2) == 'scn' else 'RG'
    return '%.3f %.3f %.3f %s' % (*rgb, op)
new, n = re.subn(r'((?:[-\d.]+\s+){4})(scn|SCN)\b', swap, stream)
print('colour operators rewritten:', n)
doc.update_stream(page.get_contents()[0], new.encode('latin1'))

# panel B is an embedded raster: remap its hues (teal -> blue, green -> green, light blue -> orange)
xref = [i[0] for i in page.get_images(full=True) if doc.extract_image(i[0])['width'] == 1222][0]
im = np.asarray(Image.open(fitz.Pixmap(doc, xref).tobytes('png') and __import__('io').BytesIO(fitz.Pixmap(doc, xref).tobytes('png'))).convert('RGBA')).astype(float) / 255
hsv = rgb_to_hsv(im[..., :3]); h, s, v = hsv[..., 0] * 360, hsv[..., 1], hsv[..., 2]
teal = (h > 175) & (h < 220) & (s > 0.3)
green = (h > 85) & (h < 165) & (s > 0.3)
light = teal & (v > 0.7)
def paint(mask, target, v_scale):
    t = rgb_to_hsv(np.array(target)[None, None, :])[0, 0]
    hsv[..., 0][mask] = t[0]; hsv[..., 1][mask] = t[1]
    hsv[..., 2][mask] = np.clip(v[mask] * v_scale, 0, 1)
paint(teal & ~light, BLUE, 0.706 / 0.50)
paint(green, GREEN, 0.627 / 0.41)
paint(light, ORANGE, 1.0 / 0.91)
rgba = np.concatenate([hsv_to_rgb(hsv), im[..., 3:]], axis=2)
Image.fromarray((rgba * 255).round().astype(np.uint8)).save(out + '_panelB.png')
page.replace_image(xref, filename=out + '_panelB.png')
# the chromosome is a second raster with a lavender band: recolour it to the light blue
xref0 = [i[0] for i in page.get_images(full=True) if doc.extract_image(i[0])['width'] > 10000][0]
print('chromosome SMask:', doc.xref_get_key(xref0, 'SMask'), doc.xref_get_key(xref0, 'ColorSpace'))
im0 = np.asarray(Image.open(__import__('io').BytesIO(doc.extract_image(xref0)['image'])).convert('RGBA')).astype(float) / 255
hsv0 = rgb_to_hsv(im0[..., :3]); h0, s0, v0 = hsv0[..., 0] * 360, hsv0[..., 1], hsv0[..., 2]
band = (h0 > 215) & (h0 < 290) & (s0 > 0.04) & (s0 < 0.6) & (v0 > 0.6)
t = rgb_to_hsv(np.array(LIGHT_BLUE)[None, None, :])[0, 0]
hsv0[..., 0][band] = t[0]; hsv0[..., 1][band] = t[1]
Image.fromarray((np.concatenate([hsv_to_rgb(hsv0), im0[..., 3:]], 2) * 255).round().astype(np.uint8)).save(out + '_chrom.png')
page.replace_image(xref0, filename=out + '_chrom.png')
print('chromosome band pixels recoloured:', int(band.sum()))
doc.save(out + '.pdf', garbage=3, deflate=True)
pix = fitz.open(out + '.pdf')[0].get_pixmap(dpi=600, alpha=False)
pix.save(out + '.png'); print('rendered', pix.width, pix.height)
