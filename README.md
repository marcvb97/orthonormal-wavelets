# Orthonormal Wavelets — Image & Tensor Compression Toolbox

A MATLAB toolbox for image and video compression using custom orthonormal wavelets. The toolbox supports two wavelet families (**OW** and **VP**), operates on both 2D colour images and 3D tensors (e.g. grayscale video), and includes baseline comparisons against MATLAB's built-in `wavedec2`/`wavedec3`.

---

## Table of Contents

- [Overview](#overview)
- [Wavelet Families](#wavelet-families)
- [File Structure](#file-structure)
- [Quick Start](#quick-start)
- [Function Reference](#function-reference)
  - [1-D Building Blocks](#1-d-building-blocks)
  - [2-D Forward & Inverse Transforms](#2-d-forward--inverse-transforms)
  - [3-D Forward & Inverse Transforms](#3-d-forward--inverse-transforms)
  - [Decomposition Wrappers](#decomposition-wrappers)
  - [Image Compression — 2D Colour](#image-compression--2d-colour)
  - [Tensor Compression — 3D](#tensor-compression--3d)
  - [Utility & Analysis](#utility--analysis)
- [Key Parameters](#key-parameters)
- [Output Structure](#output-structure)
- [Requirements](#requirements)

---

## Overview

The toolbox implements a custom wavelet transform built on a **base-3 multiresolution analysis**: each decomposition step splits a signal of length `N = 3n` into `n` scaling coefficients and `2n` wavelet coefficients. The localization degree is controlled by the parameter `theta`, which is optimised automatically by searching over a grid of values.

Compression is achieved by **hard thresholding**: a fraction `keep` of the largest-magnitude coefficients is retained and the rest are set to zero. Quality is evaluated with PSNR and SSIM.

---

## Wavelet Families

| Family | Description | Key files |
|--------|-------------|-----------|
| **OW** (Orthonormal Wavelets) | Custom orthonormal wavelets with optional initial spectral transformation (`TR2D`/`TR3D`) | `FT1step_m.m`, `IFT1step_m.m`, `FWTmatrix1_m.m`, `IFWTmatrix_m.m` |
| **VP** (VP Wavelets) | Custom wavelets with explicit normalization factors `fsca` and `fwav` | `FWT1step_m_VP.m`, `IFWT1step_m_VP.m`, `FWTmatrix1_m_VP.m`, `IFWTmatrix_m_VP.m` |
| **wavedec3** (baseline) | MATLAB built-in 3-D wavelet via `wavedec3`/`waverec3` | `imagecompression_tensor_wavedec3.m` |

---

## File Structure

```
orthonormal_wavelets/
│
├── 1-D transforms
│   ├── FT1step_m.m                   OW  forward 1-D one-step transform
│   ├── FT1step_m_fast.m              OW  forward 1-D one-step transform (optimised)
│   ├── IFT1step_m.m                  OW  inverse 1-D one-step transform
│   ├── FWT1step_m_VP.m               VP  forward 1-D one-step transform
│   ├── FWT1step_m_VP_fast.m          VP  forward 1-D one-step transform (optimised)
│   └── IFWT1step_m_VP.m              VP  inverse 1-D one-step transform
│
├── 2-D transforms
│   ├── FWTmatrix1_m.m                OW  multi-level 2-D forward transform
│   ├── FWTmatrix1_m_VP.m             VP  multi-level 2-D forward transform
│   ├── IFWTmatrix_m.m                OW  multi-level 2-D inverse transform
│   ├── IFWTmatrix_m_VP.m             VP  multi-level 2-D inverse transform
│   └── IFWTmatrix_m_old.m            (legacy version)
│
├── 3-D transforms
│   ├── FWTmatrix3D_m.m               OW  multi-level 3-D forward transform
│   ├── FWTtensor1_m_VP.m             VP  multi-level 3-D forward transform
│   ├── IFWTmatrix3D_m.m              OW  multi-level 3-D inverse transform
│   └── IFWTtensor_m_VP.m             VP  multi-level 3-D inverse transform
│
├── Decomposition wrappers (pad + transform + split coeff.)
│   ├── DEC.m                         OW  2-D decomposition wrapper
│   ├── DEC_new.m                     OW  2-D decomposition wrapper (with TR2D option)
│   ├── DEC_VP.m                      VP  2-D decomposition wrapper
│   ├── DEC_3D.m                      OW  3-D decomposition wrapper
│   ├── DEC3D.m                       OW  3-D decomposition wrapper (clean version)
│   └── DEC3D_VP.m                    VP  3-D decomposition wrapper
│
├── Initial spectral transformations
│   ├── ITR1D.m                       Inverse spectral transform (1-D)
│   ├── ITR2D.m                       Inverse spectral transform (2-D)
│   ├── ITR3D.m                       Inverse spectral transform (3-D)
│   └── compute_nu.m                  Compute nu weights for spectral transform
│
├── Image compression — 2-D colour
│   ├── imagecompression.m            Basic grayscale compression script
│   ├── imagecompression_color_OW.m   OW  colour image compression
│   ├── imagecompression_color_OW_par.m     OW  (parallelised R/G/B)
│   ├── imagecompression_color_OW_new.m     OW  with initial TR2D transformation
│   ├── imagecompression_color_OW_new_par.m OW  new + parallelised
│   ├── imagecompression_color_OW_PSNR_SSIM.m  OW  dual optimisation (MSE + SSIM)
│   ├── imagecompression_color_VP.m   VP  colour image compression
│   ├── imagecompression_color_VP_par.m     VP  (parallelised R/G/B)
│   ├── imagecompression_color_other.m      Baseline: db2 / bior3.5 via wavedec2
│   └── imagecompression_color_other_PSNR_SSIM.m  Baseline with dual optimisation
│
├── Tensor compression — 3-D
│   ├── imagecompression_tensor_OW.m        OW  tensor compression
│   ├── imagecompression_tensor_VP.m        VP  tensor compression
│   └── imagecompression_tensor_wavedec3.m  Baseline: MATLAB wavedec3
│
├── Batch processing & analysis
│   ├── compression_of_images_of_a_folder.m      Compress all images in a folder
│   ├── compression_of_images_of_a_folder_par.m  Same, parallelised
│   ├── analyze.m                     Average RESULTS_best over multiple images
│   ├── analyze_PSNR_SSIM.m           Average results (PSNR+SSIM variant)
│   ├── analyze_results.m             Script: load .mat files and compare methods
│   ├── analyze_timings.m             Compute average timing by image size
│   ├── main_analyze_timings.m        Script: run timing analysis
│   └── check_RESULTS_GLOBAL.m        Inspect theta/level distribution in results
│
└── Data utilities
    ├── convert_video_to_tensor.m     Read a video file into a 3-D grayscale tensor
    └── generate_grayscale_video.m    Generate a synthetic uncompressed grayscale video
    └── generate_random_images.m      Generate random PNG images of various sizes
```

---

## Quick Start

### Compress a single colour image (VP wavelets)

```matlab
[RESULTS_GLOBAL] = imagecompression_color_VP('path/to/image.png');
```

### Compress a single colour image (OW wavelets)

```matlab
[RESULTS_GLOBAL] = imagecompression_color_OW('path/to/image.png');
```

### Compress a 3-D tensor (VP wavelets)

The tensor can be a `.mat` file containing a 3-D double array, or a folder of grayscale PNG/BMP/JPG frames.

```matlab
[RESULTS_GLOBAL, Ifin] = imagecompression_tensor_VP('path/to/tensor.mat');
% or
[RESULTS_GLOBAL, Ifin] = imagecompression_tensor_VP('path/to/frames_folder/');
```

### Compress a 3-D tensor (OW wavelets)

```matlab
[RESULTS_GLOBAL, Ifin] = imagecompression_tensor_OW('path/to/tensor.mat');
```

### Compress a 3-D tensor (MATLAB wavedec3 baseline)

```matlab
[RESULTS_GLOBAL, Ifin] = imagecompression_tensor_wavedec3('path/to/tensor.mat');
```

### Convert a video to a tensor

```matlab
% Edit convert_video_to_tensor.m to set the video path, then run:
convert_video_to_tensor
% Produces a uint8 tensor of size (H x W x nFrames)
```

### Compress all images in a folder

```matlab
% Edit compression_of_images_of_a_folder.m to set image_dir, then run:
compression_of_images_of_a_folder
```

---

## Function Reference

### 1-D Building Blocks

| Function | Description |
|----------|-------------|
| `FT1step_m(a, theta)` | OW one-step forward transform. Returns scaling (`anew`), wavelet (`bnew`) coefficients and padding count (`nplus`). |
| `FT1step_m_fast(a, theta)` | Optimised version of `FT1step_m` using vectorised padding. |
| `IFT1step_m(a, b, theta)` | OW one-step inverse transform. |
| `FWT1step_m_VP(a, theta, fsca, fwav)` | VP one-step forward transform with explicit normalization factors. |
| `FWT1step_m_VP_fast(a, theta, fsca, fwav)` | Optimised VP forward transform. |
| `IFWT1step_m_VP(a, b, theta, fsca, fwav)` | VP one-step inverse transform. |

All 1-D functions expect **column vector** input. The input length `N1` is automatically padded to the next multiple of 3.

### 2-D Forward & Inverse Transforms

| Function | Description |
|----------|-------------|
| `FWTmatrix1_m(A, theta, kmax, dmin)` | OW multi-level 2-D forward transform. Returns cell array `Adec`, `lrow`, `lcol`. |
| `FWTmatrix1_m_VP(A, theta, kmax, dmin, fsca, fwav)` | VP multi-level 2-D forward transform. |
| `IFWTmatrix_m(Adec, lrow, lcol, theta)` | OW multi-level 2-D inverse transform. |
| `IFWTmatrix_m_VP(Adec, lrow, lcol, theta, fsca, fwav)` | VP multi-level 2-D inverse transform. |

### 3-D Forward & Inverse Transforms

| Function | Description |
|----------|-------------|
| `FWTmatrix3D_m(A, theta, kmax, dmin)` | OW multi-level 3-D forward transform. Returns `Adec`, `lrow`, `lcol`, `lslice`. |
| `FWTtensor1_m_VP(A, theta, kmax, dmin, fsca, fwav)` | VP multi-level 3-D forward transform. |
| `IFWTmatrix3D_m(Adec, lrow, lcol, lslice, theta)` | OW multi-level 3-D inverse transform. |
| `IFWTtensor_m_VP(Adec, lrow, lcol, lslice, theta, fsca, fwav)` | VP multi-level 3-D inverse transform. |

The 3-D transforms apply 1-D transforms sequentially along **columns (dim 1) → rows (dim 2) → slices (dim 3)**. The inverse applies them in the reverse order.

### Decomposition Wrappers

These functions pad the input, call the forward transform, and split the result into a scaling vector `Iscal` and wavelet vector `Iwav`.

| Function | Dims | Family |
|----------|------|--------|
| `DEC(I, theta, kstep, dmin)` | 2-D | OW |
| `DEC_new(I, theta, kstep, dmin, initial_transformation)` | 2-D | OW + TR2D option |
| `DEC_VP(I, theta, kstep, dmin, fsca, fwav)` | 2-D | VP |
| `DEC_3D(I, theta, kstep, dmin)` | 3-D | OW |
| `DEC3D(I, theta, kstep, dmin)` | 3-D | OW (clean) |
| `DEC3D_VP(I, theta, kstep, dmin, fsca, fwav)` | 3-D | VP |

**Subband structure per decomposition step:**

- **2-D:** 1 scaling block (LL) + 3 wavelet subbands (LH, HL, HH)
- **3-D:** 1 scaling block (LLL) + 7 wavelet subbands (LLH, LHL, LHH, HLL, HLH, HHL, HHH)

### Image Compression — 2D Colour

All functions take a single `filepath` argument (path to a colour image) and return `RESULTS_GLOBAL`.

| Function | Method | Notes |
|----------|--------|-------|
| `imagecompression_color_OW(filepath)` | OW | Standard |
| `imagecompression_color_OW_par(filepath)` | OW | R/G/B channels in `parfor` |
| `imagecompression_color_OW_new(filepath)` | OW | Includes initial TR2D spectral transform |
| `imagecompression_color_OW_new_par(filepath)` | OW | New + parallelised |
| `imagecompression_color_OW_PSNR_SSIM(filepath)` | OW | Dual optimisation: best MSE and best SSIM reported separately |
| `imagecompression_color_VP(filepath)` | VP | Standard |
| `imagecompression_color_VP_par(filepath)` | VP | R/G/B channels in `parfor` |
| `imagecompression_color_other(filepath)` | db2 / bior3.5 | Baseline using MATLAB `wavedec2` |
| `imagecompression_color_other_PSNR_SSIM(filepath)` | db2 | Baseline with dual optimisation |

### Tensor Compression — 3D

| Function | Method | Input |
|----------|--------|-------|
| `imagecompression_tensor_OW(filepath)` | OW | `.mat` file or folder of frames |
| `imagecompression_tensor_VP(filepath)` | VP | `.mat` file or folder of frames |
| `imagecompression_tensor_wavedec3(filepath)` | MATLAB wavedec3 | `.mat` file or folder of frames |

All tensor functions return `[RESULTS_GLOBAL, Ifin]` where `Ifin` is the reconstructed tensor after compression with the best `theta`.

### Utility & Analysis

| Function | Description |
|----------|-------------|
| `compute_nu(n, m)` | Compute the `nu` weight vector used in spectral transforms. |
| `ITR1D(a, nu)` / `ITR2D(A, theta)` / `ITR3D(A, theta)` | Inverse spectral (initial) transform in 1-D, 2-D, and 3-D. |
| `analyze(res, factor)` | Average `RESULTS_best` over a cell array of `RESULTS_GLOBAL` structs. |
| `analyze_PSNR_SSIM(res, factor)` | Same as `analyze` but for the PSNR+SSIM variant. |
| `analyze_timings(RESULTS_GLOBAL, factor, kmin)` | Compute average timing per image size and decomposition level. |
| `check_RESULTS_GLOBAL` | Script: inspect the range of theta values and decomposition levels across a batch. |
| `convert_video_to_tensor` | Script: read a video file and save as a `uint8` tensor `(H × W × nFrames)`. |
| `generate_grayscale_video` | Script: generate a synthetic uncompressed AVI and read it back as a tensor. |
| `generate_random_images(max_power)` | Generate random RGB PNG images of sizes `2^4` up to `2^max_power`. |

---

## Key Parameters

| Parameter | Description | Typical values |
|-----------|-------------|----------------|
| `theta` | Localization degree of the wavelet. Controls the overlap between scaling and wavelet spaces. Searched over a grid to minimise MSE. | `0.1 : 0.1 : 0.9` |
| `keep` | Fraction of wavelet coefficients retained after thresholding. `keep=1` = lossless; `keep=0.005` = high compression. | `[0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005]` |
| `kstep` | Maximum number of decomposition levels. | `3` to `Kmax = floor(log(min_dim)/log(3))` |
| `dmin` | Minimum dimension of the scaling sub-block; decomposition stops when any axis reaches this size. | `3` |
| `flag` | `0` = threshold all coefficients (scaling + wavelet); `1` = threshold wavelet only (preserve scaling). | `1` |
| `fsca` | Normalization factor for scaling coefficients (VP family only). | `sqrt(3)` |
| `fwav` | Normalization factor for wavelet coefficients (VP family only). | `sqrt(3/2)` |

---

## Output Structure

All compression functions return a `RESULTS_GLOBAL` struct with the following fields:

| Field | Description |
|-------|-------------|
| `name` | Input filepath |
| `RESULTS` | Matrix of all (keep, kstep, theta) combinations: `[keep, PSNR, SSIM, %retained, best_theta, decomp_levels]` |
| `RESULTS_best` | Best row per `keep` value (highest PSNR) |
| `TIMES` | Timing matrix: `[kstep, keep, theta, t_decomp, t_threshold, t_reconstruct, MSE]` |
| `size` | Original input dimensions `[N, M]` (2-D) or `[N, M, P]` (3-D) |

---

## Requirements

- **MATLAB R2019b or later** (for `squeeze`, `dct` Type-3, `ssim`)
- **Image Processing Toolbox** (for `ssim`, `imread`, `imshow`)
- **Parallel Computing Toolbox** (optional — only required for `_par` variants)
- **Wavelet Toolbox** (optional — only required for `imagecompression_tensor_wavedec3` and the `other` baselines)
