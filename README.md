# Laplacian Coordinates for Seeded Image Segmentation

This repository contains implementations of *Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation* by Casaca et al.

---

## MATLAB implementation

- Author: Wallace Casaca ([wallace.coc@gmail.com](mailto:wallace.coc@gmail.com))
- Version: Preliminary prototype
- License: See matlab/LICENSE
- Repository: https://github.com/hbatagelo/laplacianseg

### Description

This code is a preliminary version of the hard-constrained Laplacian Coordinates framework (LCH) for seeded image segmentation.

The code is very simple to run and it has been implemented and tested in [MATLAB](https://www.mathworks.com/products/matlab.html) 9.2 under Windows 10 (64-bit).
No extra toolboxes or mex-C compilations are required to run this prototype, making it easy to use and less sensitive to OS.

### Usage

- `run_me`: Runs the LC segmentation for the sample images in the `Example_*.mat` files.
- `main_interactive('sample_image.png')`: Provides an interactive interface for segmentation.
       
### Note

1. If you plan to use this code or any variant of it, please, cite the published LC papers.
2. For any question, suggestion or bug reports, please, contact [wallace.coc@gmail.com](mailto:wallace.coc@gmail.com).
