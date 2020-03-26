# Laplacian Coordinates for Seeded Image Segmentation

This repository contains implementations of [Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation](https://ieeexplore.ieee.org/document/9000902) by Casaca et al. If you plan to use the code in this repository or any variant of it, please cite the paper.

---

## MATLAB implementation

- Author: Wallace Casaca ([wallace.coc@gmail.com](mailto:wallace.coc@gmail.com))
- Version: Preliminary prototype
- License: GPLv3

### Description

This code is a preliminary version of the hard-constrained Laplacian Coordinates framework (LCH) for seeded image segmentation.

The code is very simple to run and it has been implemented and tested in [MATLAB](https://www.mathworks.com/products/matlab.html) 9.2 under Windows 10 (64-bit).
No extra toolboxes or mex-C compilations are required to run this prototype, making it easy to use and less sensitive to OS.

### Usage

- `run_me`: Runs the LC segmentation for the sample images in the `Example_*.mat` files.
- `main_interactive('sample_image.png')`: Provides an interactive interface for segmentation.
       
### Note

- For any question, suggestion or bug reports, please contact [wallace.coc@gmail.com](mailto:wallace.coc@gmail.com).

---

## C++ implementation

- Authors: Harlen Batagelo ([hbatagelo@gmail.com](mailto:hbatagelo@gmail.com)) and João Paulo Gois ([jpgois@gmail.com](mailto:jpgois@gmail.com))
- License: GPLv3

### Description

This is a prototype C++ implementation of the Laplacian Coordinates segmentation framework with support to:
- Soft-constrained, pixel-based Laplacian Coordinates (LC).
- Hard-constrained, pixel-based Laplacian Coordinates (LCH).
- Soft-constrained, superpixel-based Laplacian Coordinates (SPLC).
- Hard-constrained, superpixel-based Laplacian Coordinates (SPLCH).

### Build instructions

- Install the dependencies: [Eigen3](https://eigen.tuxfamily.org/dox/) and [OpenCV 3.4](https://opencv.org/releases/).
- Run the qmake tool from [Qt 5](https://qt.io) in the directory that contains the `lcseg.pro` file. A C++17 compliant compiler is required.

The build will generate a command-line application called `lcseg`.

We tested the build process on Linux with GCC and MAC OS with Clang. It should work on Windows but may require some tweaking of the `lcseg.pro` file, e.g., for updating the OpenCV headers and libraries paths.

### Usage

`lcseg input seeds output [options]`

where

- `input` is the name of the input image file.
- `seeds` is the name of the image file containing the foreground and background seeds. By default, red pixels (`#ff0000`) correspond to foreground seeds and blue pixels (`#0000ff`) correspond to background seeds.
- `output` is the name of the output file that will contain the segmented image.

and `[options]` can be any combination of the following:

- `--fg`: Sets the color of the foreground seeds (default is `#ff0000`).
- `--bg`: Sets the color of the background seeds (default is `#0000ff`).
- `--hard`: Uses seeds as hard labeling constraints (LCH and SPLCH).
- `--superpixel`: Uses SLIC superpixels (SPLC and SPLCH).
- `--size`: Sets superpixel size (default is 100).
- `--compactness`: Sets superpixel compactness (default is 10).
- `-b` or `--binary`: Writes output as a binary image.
- `-q` or `--quiet`: Runs in silent mode.
