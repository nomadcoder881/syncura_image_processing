# Image Segmentation & Bounding Box Detection

This repository contains a Python script for performing basic image binning, edge detection, and bounding box extraction from an image. The code is especially aimed at document images, designed to separate regions (such as lines or characters) into bounding boxes using naive statistical techniques, binning, and simplified edge detection.

## Features

- **Image binning**: Converts an image into a matrix with reduced pixel complexity using value binning
- **Naive edge detection**: Uses simple adjacent pixel checking instead of classical kernels
- **Statistical analysis**: Calculates average whitespace (horizontal and vertical) for line detection
- **Bounding box detection**: Groups contiguous "edge" pixels into bounding boxes, filtered for text or general groups
- **Collision resolution**: Merges overlapping bounding boxes until all collisions are resolved
- **Debug output**: Annotates and saves intermediate images with drawn boxes for inspection

## Requirements

- Python 3.x
- `numpy`
- `Pillow`
- `scikit-image`
- `pypng`

Install requirements with:

```bash
pip install numpy Pillow scikit-image pypng
```

## Usage

```bash
python logic.py path/to/input_image.png
```

This will produce PNG files in a `tmp/` directory, visualizing various detected bounding boxes over your input image.

**Arguments:**

- **path/to/input_image.png** – The input image you want to analyze (PNG/JPG/BMP supported via Pillow and skimage).

**Outputs:**

- `tmp/initial_blind_boxes_<uuid>.png` – Image with initial detected bounding boxes for probable text regions
- `tmp/initial_blind_boxes_raw_<uuid>.png` – Image with all detected bounding boxes (no text constraint)
- `tmp/collided_blind_boxes_<uuid>.png` – After merging overlapping boxes (fewer, larger boxes)

> **NOTE:** Default number of bins is **5**; you can alter this in the code (`gain_image_pixel_simple_matrix`).

## Main Functions

### gain_image_pixel_simple_matrix(image_path, num_bins=10)
Reads an image, bins pixel values, generates a simplified matrix indicating edge locations, computes horizontal/vertical tolerances, and removes vertical lines.

### gain_bounding_boxes(image_pixel_simple_matrix, horizontal_tolerance, vertical_tolerance, type=None)
Segments the simplified matrix into bounding boxes, optionally constraining box growth for text-like regions.

### collide_boxes(bounding_boxes)
Naively merges bounding boxes that overlap by more than 5% area, reducing duplication.

### draw_entities_on_image(image_path, entities, image_name)
Draws rectangles for detected entities directly over the original image for debugging/visualization.

## How it Works (Summary)

1. The script bins the image's pixel values to reduce noise and color complexity.
2. A naive edge is detected by differences between adjacent binned pixel averages.
3. Statistical analysis determines "tolerances" for how far to look when grouping pixels into boxes.
4. Line/grouping algorithm sweeps the matrix, gathering contiguous or nearby "edge" pixels into rectangles.
5. Overlapping/adjacent bounding boxes are recursively merged to remove duplicates.
6. Box locations are drawn over the original image and output.

## Limitations

- This is a **naive implementation**; it is much less robust than OpenCV or advanced ML-based approaches.
- Only horizontal/vertical "lines" are split—curves, rotated text, or noisy backgrounds will not perform well.
- Mostly effective on **monochrome or low-complexity document images** (e.g., printed text).

## File Structure

```
logic.py         # The main script (as provided)
README.md              # This file
tmp/                   # Output images are stored here
```

## Example

You can test the script with a simple scanned or generated document image:

```bash
python logic.py sample_document.png
```

And check the output images in the `tmp/` directory.

---
**Happy hacking!**