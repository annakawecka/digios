import os
from PIL import Image, ImageDraw, ImageFont

parent_dir = "/Users/kawecka/digios/analysis/corrected_working/plots_17O/minuit_Excorr/ang_dist"
output_pdf = "output.pdf"

dpi = 300
a4_width_in, a4_height_in = 11.7, 8.3  # inches
page_width = int(a4_width_in * dpi)
page_height = int(a4_height_in * dpi)

pages = []

try:
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Unicode.ttf", 100)
except OSError:
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Helvetica.ttc", 100)
except:
    font = ImageFont.load_default()

for folder_name in sorted(os.listdir(parent_dir)):
    folder_path = os.path.join(parent_dir, folder_name)
    if not os.path.isdir(folder_path):
        continue

    png_files = [f for f in os.listdir(folder_path) if f.lower().endswith(".png")]
    png_files.sort()

    if not png_files:
        continue

    page = Image.new("RGB", (page_width, page_height), "white")
    draw = ImageDraw.Draw(page)

    text = folder_name
    text_w, text_h = draw.textsize(text, font=font)
    draw.text(((page_width - text_w) // 2, 80), text, fill="black", font=font)

    margin = 200
    available_height = page_height - (text_h + 200)
    available_width = page_width - 2 * margin
    img_width = available_width // 3
    img_height = available_height

    x = margin
    y = text_h + 150
    
    for png in png_files[:3]:
        img_path = os.path.join(folder_path, png)

        with Image.open(img_path) as im:
            im = im.convert("RGB")
            im.thumbnail((img_width, img_height))
            page.paste(im, (x, y))
        x += img_width

    pages.append(page)

if pages:
    pages[0].save(output_pdf, save_all=True, append_images=pages[1:])
    print(f"PDF created: {output_pdf}")
else:
    print("No PNGs found.")
