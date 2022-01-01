# Laplacian Coordinates for Seeded Image Segmentation

This repository contains implementations of [Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation](https://ieeexplore.ieee.org/document/9000902) by Casaca et al., [DOI 10.1109/TPAMI.2020.2974475](https://doi.org/10.1109/TPAMI.2020.2974475). If you plan to use the code in this repository or any variant of it, please cite the paper.

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

This is a cross-platform C++ implementation of the Laplacian Coordinates segmentation framework with support to:
- Soft-constrained, pixel-based Laplacian Coordinates (LC).
- Hard-constrained, pixel-based Laplacian Coordinates (LCH).
- Soft-constrained, superpixel-based Laplacian Coordinates (SPLC).
- Hard-constrained, superpixel-based Laplacian Coordinates (SPLCH).

### Build instructions

Install the build tools and dependencies:
- [Qt 6](https://qt.io).
- A C++20-compliant compiler.
- [CMake](https://cmake.org/), [Eigen 3.3](https://eigen.tuxfamily.org/), [OpenCV 4.2](https://opencv.org/) and [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html). On Linux these can be installed using the distribution's package manager.

  Ubuntu 20.04 and later / Debian 11 and later:
  ```
  sudo apt-get install cmake libeigen3-dev libopencv-dev libsuitesparse-dev
  ```
  Fedora 34 and later:
  ```
  sudo dnf install cmake eigen3-devel opencv-devel suitesparse-devel blas-devel lapack-devel
  ```
  macOS with Homebrew:
  ```
  brew install cmake eigen opencv suite-sparse
  ```
  For instructions on how to build SuiteSparse on Windows with MSVC, refer to the [suitesparse-metis-for-windows](https://github.com/jlblancoc/suitesparse-metis-for-windows) project.

After installing the dependencies, build the project in Qt Creator:

- Open `cpp/CMakeLists.txt` in Qt Creator and configure/build the project from there.

As an alternative, build in the command line:
- Set `CMAKE_PREFIX_PATH` to the location where Qt 6 is installed (e.g. `$HOME/Qt/6.2.2/gcc_64`).
- From the `cpp` directory, run the following commands (assuming `Release` build type):

  ```
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release ..
  cmake --build . --config "Release"
  ```
If a dependency is not found during the CMake configuration, manually set the missing path using the CMake variables `Eigen3_DIR`, `OpenCV_DIR` and `SuiteSparse_DIR`.

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
