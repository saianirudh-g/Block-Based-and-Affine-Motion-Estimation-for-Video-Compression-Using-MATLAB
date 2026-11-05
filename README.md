# Block-Based and Affine Motion Estimation for Video Compression

**MATLAB implementation and comparison of integer-pixel, half-pixel, diamond-search, and four-parameter affine motion estimation**

## Project Overview

This project explores motion estimation and motion compensation for inter-frame video prediction. It compares the quality of predicted frames and the computation required by four approaches:

1. **Integer-pixel exhaustive block matching** using 8×8 and 4×4 blocks.
2. **Half-pixel motion estimation** using local subpixel refinement and bilinear interpolation.
3. **Fast diamond-search block matching** using large and small diamond search patterns.
4. **Four-parameter affine motion estimation** for translation, rotation, and scaling (extra credit).

The experiments use the **Beauty**, **Bosphorus**, and **Jockey** video test sequences, converted to **1280×720, 8-bit luma (Y)** images. The first five frames provide four consecutive reference/current-frame pairs. The diamond-search and affine comparisons cover the first three pairs.

**Domain:** Digital Signal Processing (DSP) / Video Coding and Multimedia Communications.

## Objectives

- Implement block-based P-frame prediction with a ±16-pixel displacement search range.
- Compare residual mean squared error (MSE) for 8×8 and 4×4 blocks.
- Evaluate the quality improvement from half-pixel motion vectors.
- Compare the average number of SSD evaluations for exhaustive search and diamond search.
- Compare translational motion estimation with a four-parameter affine motion model.
- Visualize original frames, predicted frames, residual images, and motion vectors.

## Technologies

- MATLAB (Image Processing Toolbox required for frame resizing)
- FFmpeg (to extract luma PNG frames from raw YUV video)
- Video test sequences: Beauty, Bosphorus, Jockey

## Project Structure

```text
Block-Based-and-Affine-Motion-Estimation/
├── HW3_complete_motion_estimation.m
├── README.md
├── raw_videos/
│   ├── Beauty.yuv
│   ├── Bosphorus.yuv
│   └── Jockey.yuv
├── frames/
│   ├── Beauty/
│   │   ├── Beauty-01.png
│   │   └── ... Beauty-05.png
│   ├── Bosphorus/
│   │   ├── Bosphorus-01.png
│   │   └── ... Bosphorus-05.png
│   └── Jockey/
│       ├── Jockey-01.png
│       └── ... Jockey-05.png
└── hw3_results/              # Generated automatically by MATLAB
    ├── Q1_Q2_MSE.csv
    ├── Q3_search_counts.csv
    ├── Q4_affine_MSE.csv
    └── *.png
```

The large raw-video files and generated results do not need to be committed to GitHub. See `.gitignore` suggestion below.

## Data Preparation (Windows)

1. Download the 8-bit YUV sequences from [Ultra Video Group test sequences](https://ultravideo.fi/#testsequences).
2. Install [FFmpeg](https://ffmpeg.org/download.html) and ensure `ffmpeg` is available in Command Prompt (or use the full path to `ffmpeg.exe`).
3. Create the folders `raw_videos`, `frames\Beauty`, `frames\Bosphorus`, and `frames\Jockey`.
4. Place the downloaded raw videos in `raw_videos` and update the commands below to match their actual filenames.

**For 1920×1080, 8-bit planar YUV 4:2:0 input**, open Command Prompt in the project directory and run:

```bat
mkdir frames\Beauty frames\Bosphorus frames\Jockey

ffmpeg -f rawvideo -pixel_format yuv420p -video_size 1920x1080 -i "raw_videos\Beauty.yuv" -vf "scale=1280:720,format=gray" -frames:v 5 "frames\Beauty\Beauty-%02d.png"
ffmpeg -f rawvideo -pixel_format yuv420p -video_size 1920x1080 -i "raw_videos\Bosphorus.yuv" -vf "scale=1280:720,format=gray" -frames:v 5 "frames\Bosphorus\Bosphorus-%02d.png"
ffmpeg -f rawvideo -pixel_format yuv420p -video_size 1920x1080 -i "raw_videos\Jockey.yuv" -vf "scale=1280:720,format=gray" -frames:v 5 "frames\Jockey\Jockey-%02d.png"
```

**Important Windows note:** In an *interactive Command Prompt* use `%02d.png`, not `%02d.png` (the doubled percent sign is only for a `.bat` file). For direct Command Prompt usage, replace each `%%02d` above with `%02d`.

If the input is **3840×2160**, change the `-video_size` setting to `3840x2160`. If it is not **YUV420p**, use the correct input pixel format; do not infer the chroma layout from the `.yuv` extension alone. The resulting files must be named `Beauty-01.png` ... `Beauty-05.png`, etc.

## Running the MATLAB Code

1. Open MATLAB and navigate to the project root (containing `HW3_complete_motion_estimation.m`).
2. Verify the input PNG files have been extracted.
3. Run:

```matlab
HW3_complete_motion_estimation
```

The program is configured for:

```matlab
cfg.names = {'Beauty','Bosphorus','Jockey'};
cfg.blockSizes = [8 4];
cfg.search = 16;
cfg.nFrames = 5;
cfg.outDir = fullfile(pwd,'hw3_results');
cfg.saveFigures = true;
cfg.runDiamond = true;
cfg.runAffine = true;
```

**Quick smoke test:** Temporarily set `cfg.names = {'Beauty'}`, `cfg.blockSizes = 8`, `cfg.nFrames = 2`, `cfg.runDiamond = false`, and `cfg.runAffine = false`; restore the original values for the full experiment.

**Performance warning:** Exhaustive search over a full 720p frame with ±16 motion range, especially for 4×4 blocks and multiple frame pairs, is computationally demanding.

## Motion Estimation Methods

### 1. Integer-Pixel Exhaustive Search

For each current-frame block, test candidate displacements \((d_x,d_y)\) in the range \([-16,16]\) in both axes and choose the candidate with minimum sum of squared differences (SSD):

\[
SSD(d_x,d_y)=\sum_{i=0}^{B-1}\sum_{j=0}^{B-1}[C(x+i,y+j)-R(x+i+d_x,y+j+d_y)]^2.
\]

The exhaustive search has **33×33 = 1,089 displacement candidates per block** in this implementation (including border positions using edge replication).

### 2. Half-Pixel Motion Estimation

Refine the best integer motion vector using half-pixel offsets and bilinear interpolation. Measure the residual MSE improvement:

\[
\Delta MSE = MSE_{integer} - MSE_{half}.
\]

### 3. Diamond Search

Use a **large diamond search pattern (LDSP)** to move toward lower SSD and a **small diamond search pattern (SDSP)** for final refinement. Record the number of unique candidate SSD evaluations and compare it to the exhaustive search:

\[
Reduction(\%) = 100\left(1-\frac{N_{diamond}}{N_{exhaustive}}\right).
\]

Diamond search generally evaluates fewer locations but may return a suboptimal motion vector.

### 4. Four-Parameter Affine Motion Estimation

Model per-pixel motion inside each block as:

\[
MV_x(x,y)=ax+by+c,\qquad MV_y(x,y)=-bx+ay+f.
\]

The parameters represent local scale/rotation and horizontal/vertical translation. Compare the affine prediction residual MSE with the translational prediction residual MSE. The supplied MATLAB script implements a simplified iterative estimator; it is **not** a complete reproduction of the HEVC encoder described in the reference paper.

## Outputs

The MATLAB script creates `hw3_results/` containing:

| Output file | Description |
| --- | --- |
| `Q1_Q2_MSE.csv` | Integer/half-pixel MSE and improvement percentage by sequence, frame pair, and block size |
| `Q3_search_counts.csv` | Mean exhaustive/diamond SSD counts per block, computational reduction, and MSE |
| `Q4_affine_MSE.csv` | Translational vs. affine MSE (8×8 blocks) |
| `*.png` | Prediction, residual, and motion-vector visualizations |

## Results and Discussion

> **Results pending:** The test video files must be processed in MATLAB before inserting numerical results. Do not treat illustrative values as measured values.

After running the code, summarize:

| Experiment | Metric to report | Findings |
| --- | --- | --- |
| Integer search, 8×8 vs. 4×4 | Residual MSE | *Add values from `Q1_Q2_MSE.csv`* |
| Integer vs. half-pixel | Residual MSE and % improvement | *Add values from `Q1_Q2_MSE.csv`* |
| Exhaustive vs. diamond | Average SSD evaluations/block, reduction %, MSE | *Add values from `Q3_search_counts.csv`* |
| Translational vs. affine | Residual MSE | *Add values from `Q4_affine_MSE.csv`* |

Interpret the results using the actual sequence characteristics, including motion magnitude, local movement, rotation/zoom, and any loss in quality from a faster search method.

## Limitations

- Uses only the luma channel, not chroma.
- Evaluates the first five frames rather than full video sequences.
- Uses a simplified half-pixel interpolation scheme and affine estimator.
- Measures distortion and search evaluations, not encoded bitrate or full rate-distortion performance.
- Results depend on raw-video input format, spatial scaling, and interpolation choices.

## References

1. ECE5578 Multimedia Communication, *Homework 3: Motion Estimation* (assignment specification).
2. Supplied course sample, `motion_estimation.m` (block matching and motion-vector plotting reference).
3. S. Zhu and K.-K. Ma, “A New Diamond Search Algorithm for Fast Block-Matching Motion Estimation,” *IEEE Transactions on Image Processing*, vol. 9, no. 2, pp. 287–290, 2000.
4. L. Li et al., “An Efficient Four-Parameter Affine Motion Model for Video Coding,” *IEEE Transactions on Circuits and Systems for Video Technology*, vol. 28, no. 8, pp. 1934–1948, Aug. 2018, doi: [10.1109/TCSVT.2017.2699919](https://doi.org/10.1109/TCSVT.2017.2699919).
5. [Ultra Video Group – Test Sequences](https://ultravideo.fi/#testsequences).

## Suggested `.gitignore`

```gitignore
raw_videos/
hw3_results/
*.yuv
*.mp4
*.hevc
```

## License and Dataset Notes

The source video sequences belong to their respective providers; consult the dataset website for permitted uses and redistribution terms. This repository does not imply any license for the source video data or the cited research paper.
