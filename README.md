# Orthonormal Wavelets — MATLAB Toolbox

This toolbox implements a family of custom orthonormal wavelets and applies them to image and video compression. The core idea is a **base-3 subsampling** scheme: each decomposition step splits N coefficients into N/3 scaling (approximation) coefficients and 2N/3 wavelet (detail) coefficients. This contrasts with classical dyadic wavelets (Daubechies, biorthogonal), which use a 1:2 split.

The localization and smoothness of the wavelets are tuned by a parameter `theta ∈ (0, 1)`, which controls the overlap degree `m`. Two main transform variants are provided throughout the toolbox:

- **OW** — the base orthonormal wavelet transform
- **VP** — a normalized variant with explicit scaling factors `fsca = sqrt(3)` and `fwav = sqrt(3/2)`

Both variants are available for 1D signals, 2D images, and 3D tensors (e.g., video or image stacks). Compression is achieved by thresholding the smallest wavelet coefficients and retaining only a fraction `keep` of the total.

---

## Quickstart

### Grayscale image compression

```matlab
% Load and prepare image
I = double(im2gray(imread('./Images/lena512.bmp')));

theta = 0.5;   % localization parameter in (0, 1)
kstep = 4;     % number of decomposition levels
dmin  = 3;     % minimum subband size
keep  = 0.05;  % retain 5% of coefficients
flag  = 1;     % 1 = threshold wavelet only, 0 = threshold everything

% Forward transform + threshold + inverse
[Idec, Iscal, Iwav, lrow, lcol] = DEC(I, theta, kstep, dmin);
[~, ~, Icomp]                   = THR(Idec, Iscal, Iwav, keep, lrow, lcol, flag);
Irec = IFWTmatrix_m(Icomp, lrow, lcol, theta);
Irec = Irec(1:size(I,1), 1:size(I,2));

mse = mean((I(:) - Irec(:)).^2);
fprintf('MSE = %.4f,  PSNR = %.2f dB\n', mse, 10*log10(255^2/mse));
```

### Color image compression (sweep over theta)

```matlab
% Run the full OW color pipeline on any image file;
% returns MSE/PSNR/SSIM for all (theta, keep, level) combinations
RESULTS = imagecompression_color_OW('path/to/image.png');
```

### 3D tensor / video compression (VP pipeline)

```matlab
% Accepts a .mat file (variable I of size N x M x P)
% or a folder of grayscale PNG/BMP frames
[RESULTS, Ifin] = imagecompression_tensor_VP('path/to/movie.mat');
```

### Batch processing over a folder of images

```matlab
% Sequential
compression_of_images_of_a_folder;   % edit image_dir inside the script

% Parallel (uses parfor)
compression_of_images_of_a_folder_par;
```

---

## Repository Structure

```
orthonormal_wavelets/
├── Images/                  Test images (lena512, peppers512, baboon512)
├── FT1step_m.m              Core 1D forward transform step
├── IFT1step_m.m             Core 1D inverse transform step
├── FWTmatrix1_m.m           Multi-level 2D forward transform
├── IFWTmatrix_m.m           Multi-level 2D inverse transform
├── DEC.m                    High-level 2D decomposition
├── THR.m                    Wavelet coefficient thresholding
├── imagecompression*.m      Compression pipelines (color, tensor, baselines)
└── ...                      (see full file list below)
```

---

## File Reference

### Core 1D Transform Steps

The entire library builds on two primitive 1D functions. Each handles a single decomposition or reconstruction step on a vector of arbitrary length, padding to the nearest multiple of 3 as needed.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `FT1step_m.m` | Forward step: maps a length-N vector to N/3 scaling + 2N/3 wavelet coefficients. | — |
| `FT1step_m_fast.m` | Vectorized, optimized version of `FT1step_m`. | — |
| `IFT1step_m.m` | Inverse step: reconstructs the original vector from scaling and wavelet coefficients. | — |
| `FWT1step_m_VP.m` | Forward step, VP variant — applies normalization factors `fsca` and `fwav`. | — |
| `FWT1step_m_VP_fast.m` | Vectorized version of `FWT1step_m_VP`. | — |
| `FWT1step_m_VP_old.m` | Legacy version of `FWT1step_m_VP`, kept for reference. | `FWT1step_m_VP` |
| `IFWT1step_m_VP.m` | Inverse VP step corresponding to `FWT1step_m_VP`. | — |

### Multi-Level 2D Transforms

These iterate the 1D step along rows and columns of a matrix up to `kmax` times, stopping when the approximation submatrix shrinks below `dmin`.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `FWTmatrix1_m.m` | Multi-level 2D forward transform (OW). Returns cell array `Adec`. | `FT1step_m` |
| `FWTmatrix1_m_VP.m` | Multi-level 2D forward transform (VP). | `FWT1step_m_VP` |
| `IFWTmatrix_m.m` | Multi-level 2D inverse transform (OW). | `IFT1step_m` |
| `IFWTmatrix_m_VP.m` | Multi-level 2D inverse transform (VP). | `IFWT1step_m_VP` |
| `IFWTmatrix_m_old.m` | Legacy version of `IFWTmatrix_m`. | `IFT1step_m`, `IFWTmatrix_m` |

### 2D Decomposition / Reconstruction (High-Level)

These wrap the multi-level transforms with padding logic and extract flat coefficient vectors `Iscal` and `Iwav` for easy thresholding.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `DEC.m` | Pads image to multiple of 3, decomposes, returns `Iscal` and `Iwav`. | `FWTmatrix1_m` |
| `DEC_VP.m` | VP-normalized variant of `DEC`. | `FWTmatrix1_m_VP` |
| `DEC_new.m` | Like `DEC`, but optionally applies a `TR2D` pre-transform first. | `FWTmatrix1_m`, `TR2D` |
| `REC_new.m` | Inverse of `DEC_new`; optionally applies `ITR2D` post-transform. | `IFT1step_m`, `ITR2D` |

### 3D (Tensor) Transforms

All 3D functions extend their 2D counterparts by applying the 1D step along a third axis (slices). They support video tensors stored as `.mat` files or folders of frames.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `FWTmatrix3D_m.m` | Multi-level 3D forward transform (OW). | `FT1step_m` |
| `FWTtensor1_m_VP.m` | Multi-level 3D forward transform (VP). | `FWT1step_m_VP` |
| `IFWTmatrix3D_m.m` | Multi-level 3D inverse transform (OW). | `IFT1step_m` |
| `IFWTtensor_m_VP.m` | Multi-level 3D inverse transform (VP). | `IFWT1step_m_VP` |
| `DEC3D.m` | High-level 3D decomposition (OW). | `FWTmatrix3D_m` |
| `DEC3D_VP.m` | High-level 3D decomposition (VP). | `FWTtensor1_m_VP` |
| `DEC_3D.m` | Alternative 3D decomposition treating slices independently. | `FWTmatrix1_m` |

### Thresholding

Thresholding is how compression is achieved: coefficients below a threshold are zeroed, retaining only a `keep` fraction of the total.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `THR.m` | Thresholds 2D wavelet coefficients. `flag=0`: all coefficients; `flag=1`: wavelet only. | — |
| `THR3D.m` | 3D version of `THR`. | — |
| `THR3D_VP.m` | VP-normalized 3D thresholding. | `THR` |
| `THR_other.m` | Thresholding for MATLAB built-in `wavedec`/`wavedec2` output (baseline comparisons). | — |

### Auxiliary DCT-based Transform (TR / ITR)

Some compression variants apply a DCT-based change of basis before the wavelet transform to improve energy compaction. `TR2D` and `ITR2D` are invertible and self-contained.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `compute_nu.m` | Computes the weight vector `nu(r)` controlling the TR reweighting. | — |
| `TR1D.m` | 1D forward TR transform (DCT-domain reweighting). | — |
| `ITR1D.m` | Inverse of `TR1D`. | — |
| `TR2D.m` | 2D forward TR transform (rows then columns). | `TR1D`, `compute_nu` |
| `ITR2D.m` | Inverse of `TR2D`. | `ITR1D`, `compute_nu` |
| `TR3D.m` | 3D forward TR transform. | `TR1D`, `compute_nu` |
| `ITR3D.m` | Inverse of `TR3D`. | `ITR1D`, `compute_nu` |

### Compression Pipelines

#### Grayscale

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `imagecompression.m` | Script: sweeps `theta` and `keep`, reports MSE for each combination. | `DEC`, `IFWTmatrix_m`, `ITR2D`, `THR`, `TR2D` |
| `our_imagecompression.m` | Like `imagecompression.m` but uses the VP pipeline. | `DEC`, `IFWTmatrix_m`, `THR` |

#### Color (R/G/B channels processed independently)

All functions below accept a file path and return a `RESULTS_GLOBAL` struct with MSE, PSNR, SSIM, and timings for each `(theta, keep, level)` combination.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `imagecompression_color_OW.m` | OW pipeline, RGB. | `DEC`, `IFWTmatrix_m`, `ITR2D`, `THR`, `TR2D` |
| `imagecompression_color_OW_par.m` | OW pipeline, RGB, parallelized over channels. | `DEC`, `IFWTmatrix_m`, `ITR2D`, `THR`, `TR2D` |
| `imagecompression_color_OW_new.m` | OW pipeline with TR2D pre-transform. | `DEC_new`, `REC_new`, `THR` |
| `imagecompression_color_OW_new_par.m` | Parallel version of `imagecompression_color_OW_new`. | `DEC_new`, `REC_new`, `THR` |
| `imagecompression_color_OW_PSNR_SSIM.m` | OW pipeline, optimizes independently for MSE and SSIM. | `DEC`, `IFWTmatrix_m`, `ITR2D`, `THR`, `TR2D` |
| `imagecompression_color_OW_YCbCr.m` | OW pipeline, YCbCr color space. | `DEC`, `IFWTmatrix_m`, `ITR2D`, `THR`, `TR2D` |
| `imagecompression_color_VP.m` | VP pipeline, RGB. | `DEC_VP`, `IFWTmatrix_m_VP`, `THR` |
| `imagecompression_color_VP_par.m` | VP pipeline, RGB, parallelized. | `DEC_VP`, `IFWTmatrix_m_VP`, `THR` |
| `imagecompression_color_VP_YCbCr.m` | VP pipeline, YCbCr color space. | `DEC_VP`, `IFWTmatrix_m_VP`, `THR` |
| `imagecompression_color_other.m` | Baseline using MATLAB `wavedec2` (e.g., `db2`, `bior3.5`). | `THR_other` |
| `imagecompression_color_other_PSNR_SSIM.m` | Baseline with independent MSE/SSIM optimization. | `THR_other` |
| `imagecompression_color_other_YCbCr.m` | Baseline in YCbCr color space. | `THR_other` |

#### Tensor / Video

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `imagecompression_tensor_OW.m` | OW pipeline on a 3D tensor. Accepts `.mat` or folder of frames. | `DEC3D`, `IFWTmatrix3D_m`, `ITR3D`, `THR3D`, `TR3D` |
| `imagecompression_tensor_VP.m` | VP pipeline on a 3D tensor. | `DEC3D_VP`, `IFWTtensor_m_VP`, `THR3D_VP` |
| `imagecompression_tensor_wavedec3.m` | Baseline using MATLAB `wavedec3`/`waverec3`. | — |

### Batch Processing

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `compression_of_images_of_a_folder.m` | Runs `imagecompression_color_OW_YCbCr` on every image in a directory. | `imagecompression_color_OW_YCbCr` |
| `compression_of_images_of_a_folder_par.m` | Same, parallelized over images with `parfor`. | `imagecompression_color_OW_par` |

### Analysis and Visualization

After running compression experiments, results are saved to `.mat` files and analyzed with these functions.

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `analyze.m` | Averages MSE results across a set of images. | — |
| `analyze_PSNR_SSIM.m` | Same but for PSNR/SSIM-optimized results. | — |
| `analyze_timings.m` | Extracts and averages timing data by image size. | — |
| `main_analyze_timings.m` | Loads multiple result files and plots timing comparisons. | `analyze_timings` |
| `check_RESULTS_GLOBAL.m` | Diagnostic: reports level range and `theta` distribution in a results struct. | — |
| `run_analysis.m` | Loads OW, VP, and baseline results and plots comparative MSE curves. | `analyze` |
| `run_analysis_YCbCr.m` | Same as `run_analysis.m` for YCbCr results. | `analyze` |
| `PSNR_in_function_of_theta_and_level.m` | Plots PSNR vs. `theta` and decomposition level. | — |

### Utilities and Data Generation

| File | Description | Local Dependencies |
|------|-------------|--------------------|
| `generate_random_images.m` | Generates random RGB test images of sizes `2^4` to `2^max_power`. | — |
| `generate_grayscale_video.m` | Generates a synthetic grayscale video and reads it back as a 3D tensor. | — |
| `convert_video_to_tensor.m` | Reads an MP4, converts to grayscale, saves as a `.mat` tensor. | — |
| `write_tensor_into_video.m` | Writes a 3D uint8 tensor to an uncompressed AVI file at 30 fps. | — |
| `timings_1D.m` | Benchmarks `FT1step_m`, `FWT1step_m_VP`, and `wavedec` across signal lengths. | `FT1step_m`, `FWT1step_m_VP` |

### Unit Tests

| File | What it tests | Local Dependencies |
|------|---------------|--------------------|
| `test_FT1step_IFT1step.m` | Perfect reconstruction: `IFT1step_m( FT1step_m(a) ) ≈ a` | `FT1step_m`, `IFT1step_m` |
| `test_FWTmatrix1_IFWTmatrix.m` | Perfect reconstruction on a 1000×1000 matrix. | `FWTmatrix1_m`, `IFWTmatrix_m` |
| `test_FWT_IFWTtensor.m` | Six round-trip tests on VP 3D tensors (aligned, padded, asymmetric, large, constant, quasi-2D). | `FWTtensor1_m_VP`, `IFWTtensor_m_VP` |
| `test_TR1D_ITR1D.m` | `ITR1D( TR1D(f) ) ≈ f` | `TR1D`, `ITR1D`, `compute_nu` |
| `test_TR2D_ITR2D.m` | `ITR2D( TR2D(F) ) ≈ F` on a 1000×900 matrix. | `TR2D`, `ITR2D` |
| `test_imagecompression_tensor_VP.m` | Full pipeline test: `DEC3D_VP → THR3D_VP → IFWTtensor_m_VP`, covering lossless, quality sweeps, and both threshold modes. | `DEC3D_VP`, `IFWTtensor_m_VP`, `THR3D_VP` |
| `test_roundtrip3D.m` | `DEC3D → IFWTmatrix3D_m` round-trip on five tensor sizes including non-multiples of 3. | `DEC3D`, `IFWTmatrix3D_m` |
| `test01.m` | Diagnostic: checks `TR1D` output dimensions for row vs. column inputs. | `TR1D`, `compute_nu` |

---

## Test Images

Three standard benchmark images are included in `Images/`:

| File | Size | Color |
|------|------|-------|
| `lena512.bmp` | 512×512 | Grayscale |
| `peppers512.tiff` | 512×512 | Color |
| `baboon512.tiff` | 512×512 | Color |

---

## Dependency Graph

The diagram below shows how project files call each other (MATLAB built-ins omitted). Leaf nodes at the top have no local dependencies and form the foundation of the library; pipelines at the bottom assemble everything together.

```
Leaf functions (no local dependencies)
───────────────────────────────────────
FT1step_m           IFT1step_m
FT1step_m_fast      IFWT1step_m_VP
FWT1step_m_VP       TR1D / ITR1D
FWT1step_m_VP_fast  compute_nu
THR / THR3D         THR_other

1D → Multi-level 2D transforms
────────────────────────────────
FT1step_m       ──► FWTmatrix1_m    ──► FWTmatrix3D_m
FWT1step_m_VP   ──► FWTmatrix1_m_VP ──► FWTtensor1_m_VP
IFT1step_m      ──► IFWTmatrix_m    ──► IFWTmatrix3D_m
IFWT1step_m_VP  ──► IFWTmatrix_m_VP ──► IFWTtensor_m_VP

TR / ITR chain
──────────────
compute_nu, TR1D  ──► TR2D, TR3D
compute_nu, ITR1D ──► ITR2D, ITR3D

THR chain
─────────
THR ──► THR3D_VP

High-level decomposition / reconstruction
───────────────────────────────────────────
FWTmatrix1_m              ──► DEC, DEC_3D
FWTmatrix1_m_VP           ──► DEC_VP
FWTmatrix1_m + TR2D       ──► DEC_new
FWTmatrix3D_m             ──► DEC3D
FWTtensor1_m_VP           ──► DEC3D_VP
IFT1step_m + ITR2D        ──► REC_new

2D OW compression pipelines
─────────────────────────────
DEC + IFWTmatrix_m + TR2D + ITR2D + THR
  ──► imagecompression
  ──► imagecompression_color_OW  (+_par, +_PSNR_SSIM, +_YCbCr)

DEC_new + REC_new + THR
  ──► imagecompression_color_OW_new  (+_par)

2D VP compression pipeline
───────────────────────────
DEC_VP + IFWTmatrix_m_VP + THR
  ──► imagecompression_color_VP  (+_par, +_YCbCr)

2D baseline (built-in wavelets)
────────────────────────────────
THR_other
  ──► imagecompression_color_other  (+_PSNR_SSIM, +_YCbCr)

3D OW tensor pipeline
──────────────────────
DEC3D + IFWTmatrix3D_m + TR3D + ITR3D + THR3D
  ──► imagecompression_tensor_OW

3D VP tensor pipeline
──────────────────────
DEC3D_VP + IFWTtensor_m_VP + THR3D_VP
  ──► imagecompression_tensor_VP  (+test)

Batch processing
─────────────────
imagecompression_color_OW_YCbCr ──► compression_of_images_of_a_folder
imagecompression_color_OW_par   ──► compression_of_images_of_a_folder_par

Analysis
─────────
analyze         ──► run_analysis, run_analysis_YCbCr
analyze_timings ──► main_analyze_timings
```

---

## MATLAB Toolbox Requirements

| Toolbox | Functions used | Required by |
|---------|---------------|-------------|
| **Signal Processing** | `dct` | `FT1step_m`, `FT1step_m_fast`, `FWT1step_m_VP`, `FWT1step_m_VP_fast`, `IFT1step_m`, `IFWT1step_m_VP`, `TR1D`, `ITR1D` |
| **Wavelet** | `wavedec2`, `waverec2`, `wavedec3`, `waverec3` | `imagecompression_color_other*`, `imagecompression_tensor_wavedec3`, `timings_1D` |
| **Image Processing** | `imread`, `imshow`, `im2gray`, `rgb2ycbcr`, `ycbcr2rgb`, `ssim` | All `imagecompression_color_*` and `imagecompression_tensor_*` functions |
| **Parallel Computing** | `parfor` | `imagecompression_color_OW_par`, `imagecompression_color_OW_new_par`, `imagecompression_color_VP_par`, `compression_of_images_of_a_folder_par` |

The Parallel Computing Toolbox is optional — all `_par` functions have sequential equivalents.
