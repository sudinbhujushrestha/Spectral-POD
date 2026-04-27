# Spectral Proper Orthogonal Decomposition of Wind-over-Waves DNS Data

A comprehensive analysis of coherent turbulent structures in a wind-over-waves boundary layer using Spectral Proper Orthogonal Decomposition (SPOD).

## Overview

This project applies **Spectral Proper Orthogonal Decomposition (SPOD)** to Direct Numerical Simulation (DNS) data of a wind-over-waves boundary layer. SPOD is a modal decomposition technique that extracts spatially and temporally coherent structures at specific frequencies, making it ideal for identifying physically meaningful phenomena in turbulent flows.

### What is SPOD?

Standard POD identifies energetic spatial structures but ignores temporal evolution, mixing multiple physical mechanisms operating at different frequencies. SPOD solves this by:
- Applying a temporal Fourier transform to the space-time correlation tensor
- Producing modes that are coherent in **both space and time**
- Isolating structures at specific frequencies (e.g., vortex shedding, traveling waves)

## Dataset

- **Grid resolution**: 192 × 192 × 65 points
- **Snapshots**: 500 temporal samples
- **Velocity components**: U (streamwise), V (spanwise), W (wall-normal)
- **Frequency bands**: 26 resolved frequencies from Δf = 2.667 Hz
- **Total TKE**: 1.512 × 10⁻¹

## Key Results

### Spectral Energy Distribution

- **Global peak** at f = 2.67 Hz where the leading mode captures **33.9%** of local spectral energy
- **Coherent–incoherent transition** identified at f ≈ 18.67 Hz
- **Three frequency regimes**:
  - Low frequencies (f ≲ 5.33): Strong rank-1 dominance
  - Intermediate (5.33 ≲ f ≲ 16.0): Meaningful energy advantage
  - High frequencies (f ≳ 18.67): Broadband incoherent fluctuations

### Spatial Structure

**Mode 1 progression across frequencies:**
- **f = 2.67 Hz**: Large, smooth structure with strong streamwise dominance (U:V:W = 68:22:10)
- **f = 8.0 Hz**: Structure contracts, forward-leaning inclination
- **f = 16.0 Hz**: Surface-locked anomaly with W = 35.5% (only frequency where W exceeds V)
- **f = 53.33 Hz**: Near-isotropic partition (U:V:W ≈ 27:40:33)

### Component-Energy Cascade

A systematic anisotropy cascade is observed:
- **Coherent regime** (f ≲ 5.33): Strongly U-dominated (70:21:10)
- **Intermediate coherent range**: U share drops sharply, W peaks
- **Incoherent plateau** (f ≥ 18.67): Near-isotropic, mildly V-biased small-scale fluctuations
- U crosses equipartition (1/3) at f = 26.67 Hz

### Wave Character

Real-part energy fraction ≈ **50%** for all 12 mode–frequency combinations, confirming that **all structures are downstream-propagating waves** (not standing waves).

## Project Structure

```
SPOD_project/
├── README.md                          # This file
├── Report/
│   ├── SPOD_Report.pdf               # Full technical report (12 pages)
│   ├── SPOD_Report.tex               # LaTeX source
│   └── figures/
│       ├── spod_spectrum.png         # Eigenvalue spectrum (2 panels)
│       ├── spod_mode1_multifreq_xzslices.png
│       ├── spod_mode2_multifreq_xzslices.png
│       ├── spod_mode3_multifreq_xzslices.png
│       └── spod_zprofile_uvw.png     # z-profile & component cascade
└── Code/
    └── projectv2_final.ipynb         # Complete analysis pipeline
```

## Methodology

### Welch's Method for CSD Estimation

The full DNS record is partitioned into **Nblk = 19 overlapping blocks** with **50% overlap**. Each block is:
1. Windowed with a **Hann window** to mitigate spectral leakage
2. Fourier-transformed independently along the temporal dimension
3. Assembled into frequency-domain snapshot matrices X(fⱼ)

This provides a controlled trade-off between frequency resolution and statistical noise reduction.

### Weighted Method of Snapshots

The SPOD eigenvalue problem is solved via the **reduced Method of Snapshots**:

```
R Ψ = Ψ Λ
```

where:
- R = (1/Nblk) X* W X  (reduced correlation matrix)
- W = spatial weight matrix (accounts for non-uniform grid clustering)
- Ψ = reduced eigenvectors (orthonormal)
- Λ = modal energies

Spatial modes are reconstructed as:
```
Φ = (Nblk Λ)^(-1/2) X Ψ
```

### Spatial Weighting

Critical for non-uniform grids (as in this DNS):
- **Horizontal**: uniform spacing (Δx, Δy)
- **Vertical**: clustered near surface (non-uniform z)
- Weight matrix uses **trapezoidal rule**: Δzₖ = ½(zₖ₊₁ - zₖ₋₁)
- Verification: sum of all weights = 9π² ≈ 88.83 ✓

### Phase Alignment

Complex SPOD modes are rotated to canonical phase by maximizing real-part energy:

```
φ* = (1/2) arg(Σₙ Wₙ Φₙ²)
```

Real-part energy fraction characterizes wave type:
- **~100%**: Standing wave
- **~50%**: Purely traveling (propagating) wave
- **~30%**: Mixed standing + traveling

## How to Use

### Requirements

```
Python 3.8+
numpy
scipy
matplotlib
jupyter
```

### Running the Analysis

1. **Open the notebook**:
   ```bash
   jupyter notebook Code/projectv2_final.ipynb
   ```

2. **Data Setup** (Step 1-2):
   - Load DNS binary data from Fortran (`ACCESS='STREAM'`)
   - Rearrange into time-domain snapshot matrix Q

3. **Preprocessing** (Step 3-4):
   - Compute weight matrix W from non-uniform grid
   - Apply Hann windowing, block segmentation, fluctuation extraction

4. **Spectral Estimation** (Step 5-7):
   - FFT each overlapping block
   - Assemble frequency-domain matrices X(fⱼ)
   - Solve reduced eigenvalue problem one frequency at a time

5. **Mode Reconstruction** (Step 8-9):
   - Reconstruct 3D spatial modes Φ(fⱼ)
   - Verify orthonormality ║Φ* W Φ - I║ ≈ 10⁻¹⁵

6. **Diagnostics** (Step 10-13):
   - Eigenvalue spectrum (Figure 2)
   - Vertical slices at 4 target frequencies (Figures 3-5)
   - Wall-normal profiles & component cascade (Figure 6)

### Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ntime` | 500 | Number of snapshots |
| `nblk` | 19 | Number of overlapping blocks |
| `block_size` | 50 | Snapshots per block |
| `overlap` | 25 | Snapshots of overlap (50%) |
| `nfreq` | 26 | Resolved frequency bins |
| `Δf` | 2.667 | Frequency resolution (Hz) |
| `f_Ny` | 66.67 | Nyquist frequency (Hz) |

## Figures

### Figure 2: SPOD Eigenvalue Spectrum
- **Left panel**: Log-energy vs. linear frequency with global peak at f = 2.67 Hz
- **Right panel**: Log-log axes showing f⁻⁵/³ Kolmogorov reference (inertial subrange)
- All 19 modes ranked by energy

### Figures 3-5: Vertical Structure of Leading Modes
Phase-aligned real parts of modes 1, 2, 3 at four target frequencies:
- f = 2.67 Hz (primary coherent peak)
- f = 8.0 Hz (intermediate coherent range)
- f = 16.0 Hz (highest rank separation anomaly)
- f = 53.33 Hz (incoherent turbulent reference)

Each figure shows U, V, W components stacked vertically.

### Figure 6: Wall-Normal Distribution & Component Cascade
- **(a) Left**: Energy profiles E(z) at four target frequencies
- **(b) Right**: U, V, W component-energy shares across all 26 frequencies with coherent–incoherent transition marked

## References

The analysis is based on the following foundational work:

1. **Schmidt, O.T. and Colonius, T.** (2020)  
   "Guide to Spectral Proper Orthogonal Decomposition"  
   *AIAA Journal*, Vol. 58, No. 3, pp. 1023–1033  
   https://doi.org/10.2514/1.J058809

2. **Towne, A., Schmidt, O.T., and Colonius, T.** (2018)  
   "Spectral Proper Orthogonal Decomposition and its Relationship to Dynamic Mode Decomposition and Resolvent Analysis"  
   *Journal of Fluid Mechanics*, Vol. 847, pp. 821–867  
   https://doi.org/10.1017/jfm.2018.283

3. **Lumley, J.L.** (1970)  
   *Stochastic Tools in Turbulence*  
   Academic Press, New York

4. **Welch, P.D.** (1967)  
   "The Use of Fast Fourier Transform for the Estimation of Power Spectra: A Method Based on Time Averaging Over Short, Modified Periodograms"  
   *IEEE Transactions on Audio and Electroacoustics*, Vol. 15, No. 2, pp. 70–73  
   https://doi.org/10.1109/TAU.1967.1161901

5. **Bartlett, M.S.** (1948)  
   "Smoothing Periodograms from Time Series with Continuous Spectra"  
   *Nature*, Vol. 161, No. 4096, pp. 686–687  
   https://doi.org/10.1038/161686a0

6. **Harris, F.J.** (1978)  
   "On the Use of Windows for Harmonic Analysis with the Discrete Fourier Transform"  
   *Proceedings of the IEEE*, Vol. 66, No. 1, pp. 51–83  
   https://doi.org/10.1109/PROC.1978.10837

7. **Taira, K., Brunton, S.L., Dawson, S.T.M., Rowley, C.W., Colonius, T., McKeon, B.J., Schmidt, O.T., Gordeyev, S., Theofilis, V., and Ukeiley, L.S.** (2017)  
   "Modal Analysis of Fluid Flows: An Overview"  
   *AIAA Journal*, Vol. 55, No. 12, pp. 4013–4041  
   https://doi.org/10.2514/1.J056060

8. **Sirovich, L.** (1987)  
   "Turbulence and the Dynamics of Coherent Structures. Part I: Coherent Structures"  
   *Quarterly of Applied Mathematics*, Vol. 45, No. 3, pp. 561–571

## Report

A complete **12-page technical report** is included in `Report/SPOD_Report.pdf`. The report contains:

- **Theoretical Formulation** (Section 2): Full mathematical derivation of weighted SPOD with flowchart
- **Results** (Section 3): Spectral analysis, mode structure, wall-normal distribution, component cascade
- **Implementation** (Section 4): Memory-efficient discrete pipeline for large DNS datasets
- **Appendix**: Complete eigenvalue table for all 26 frequencies

## Author

**Sudin Bhuju Shrestha**  
(sbhujush@cougarnet.uh.edu, sudinbhujushrestha95@gmail.com)
PhD Researcher, Direct Numerical Simulation of Wind-over-Waves Boundary Layers 
University of Houston,
Department of Mechanical and Aerospace Engineering

April 2026

## License

This project is provided as-is for research and educational purposes.

## Citation

If you use this code or results in your research, please cite the technical report and the foundational SPOD papers (Schmidt & Colonius 2020, Towne et al. 2018).

---

**Questions or issues?** Please open an issue or contact the author.
