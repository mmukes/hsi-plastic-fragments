# Hyperspectral images of post-consumer plastic fragments (PET, HDPE, PP, PS)

This folder holds ten near-infrared hyperspectral images (cubes) of shredded post-consumer plastic fragments and their class labels. It also holds the same data cut into one small hyperspectral image per fragment: 200 images, 50 per class.

The cubes are the first ten of the 54-cube dataset used in:

> M. K. Singh, S. Hait and A. Thakur, "Hyperspectral imaging-based classification of post-consumer thermoplastics for plastics recycling using artificial neural network", *Process Safety and Environmental Protection* 179 (2023) 593–602.

## Contents

| Path | What it is |
|---|---|
| `HSI_images/scan0001.bil` … `scan0010.bil`, each with its `.hdr` | 10 cubes, each with 20 fragments in a 4 × 5 grid (5 of each class) |
| `HSI_images/1hotcode.csv` | The class of every fragment: 200 rows × 4 columns (PET, HDPE, PP, PS), one-hot, no header |
| `fragments/PET/`, `fragments/HDPE/`, `fragments/PP/`, `fragments/PS/` | One hyperspectral image per fragment, 50 in each folder |
| `fragments/fragment_index.csv` | One line per fragment: file, class, source cube, fragment number, scan number, crop position and size, number of fragment pixels |
| `extract_fragments.m` | The MATLAB code that made `fragments/` from `HSI_images/` |

## Cubes

- **Camera:** Resonon Pika NIR-320 (push-broom), 168 bands, 883.9–1719.3 nm. The band centres are listed in each `.hdr`.
- **Format:** ENVI, band-interleaved-by-line (BIL), unsigned 16-bit, little-endian. Each cube is 400 lines × 320 samples × 168 bands.
- **Values:** reflectance × 16383 (`reflectance scale factor = 16383` in the `.hdr`).
- **Pixel size on the belt:** about 0.66 mm across the scan and 0.69 mm along it.

## Labels

Row (N − 1) × 20 + k of `1hotcode.csv` is fragment k of `scanN`.

The fragments of a cube are numbered in the order that MATLAB's `bwlabel` finds them after segmentation, scanning the image column by column from the left. `extract_fragments.m` uses exactly this segmentation and numbering.

## Fragment images

Fragment k of `scanNNNN` is saved in the folder of its class as:

- **`scanNNNN_fragKK.bil` and `.hdr`:** the rectangle around the fragment plus a 2-pixel margin, with all 168 bands. The format and values are the same as the cube's. The `.hdr` records where the rectangle was cut from.
- **`scanNNNN_fragKK_mask.png`:** shows which pixels of the rectangle belong to the fragment (white) and which are belt (black).

The size of each rectangle (lines and samples) is in its `.hdr` and in `fragment_index.csv`.

MATLAB:

```matlab
f = 'fragments/PET/scan0001_frag01';            % 56 lines x 38 samples (see the .hdr)
X = multibandread([f '.bil'], [56 38 168], 'uint16', 0, 'bil', 'ieee-le');
M = imread([f '_mask.png']);                    % true = fragment
R = double(X) / 16383;                          % reflectance
P = reshape(R, [], 168);                        % pixels x bands
s = mean(P(M(:), :), 1);                        % mean spectrum of the fragment
```

Python:

```python
import numpy as np
from PIL import Image
f = 'fragments/PET/scan0001_frag01'
lines, samples = 56, 38                         # from the .hdr
X = np.fromfile(f + '.bil', dtype='<u2').reshape(lines, 168, samples).transpose(0, 2, 1)
M = np.array(Image.open(f + '_mask.png')) > 0   # True = fragment
R = X / 16383.0                                 # reflectance, lines x samples x bands
s = R[M].mean(0)                                # mean spectrum of the fragment
```

## Re-making the fragments

1. Open `extract_fragments.m` in MATLAB and click **Run**. It needs the Image Processing Toolbox.
2. Select the `HSI_images` folder.

The `fragments` folder is written next to the script, and any existing files in it are overwritten.

## Note for train/test splits

The fragments were re-imaged in different positions and orientations from one cube to the next, so the same physical piece of plastic can appear in more than one cube. Splitting by cube (for example, training on eight cubes and testing on two) gives a stricter test than splitting the 200 fragment images at random.
