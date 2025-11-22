## A method to quantify 3D variation in photosynthetic ability in maize canopies
The pipeline includes:
- 3D canopy structure from point clouds  
- RGB-based SPAD estimation  
- C4 leaf photosynthesis parameterization  
- Canopy-scale photosynthesis simulation
  <img width="3341" height="2154" alt="Fig1" src="https://github.com/user-attachments/assets/8658969d-1957-4aec-88d0-1027c324d3bd" />

## Refrence 
Liu, F. et al. A method to quantify 3D variation in photosynthetic ability in maize canopies. Plant Physiol. 199, kiaf557 (2025). https://doi.org/10.1093/plphys/kiaf557

The methods can be freely used for academic purposes. For commercial purposes, a special license is required.
## Installation Note

Before running the code, please make sure to install the required dependency:

- [MATLAB Bayesian Estimation Toolbox](https://github.com/NilsWinter/matlab-bayesian-estimation)

You can install it by cloning the repository and adding it to your MATLAB path:

# C4-photosynthesis modules

## 1. `point_clouds_color_correction_alignment/`

**Goal:** Pre-process and align 3D point clouds so that geometry and RGB colors are consistent across views and plants.

### Key functions / tasks

- **Point cloud I/O**
  - Read raw point clouds from NeRF or SfM  exports (`.ply` etc.).
  - Save cleaned / aligned clouds for downstream modules.

- **Color correction & normalization**
  - Apply color correction based on the reference panel.
  - Normalize RGB values across views / sessions.

- **Geometric alignment & registration**
  - Register multiple views of the same plant into a common coordinate system.
  - Optionally align multiple plants at different time points for comparison.

- **Filtering & cropping**
  - Remove background points and noise.
  - Crop to regions of interest (i.e. a single plant ).

**Output:**  
Clean, aligned, color-corrected point clouds used by:

- `../compute_dgci_and_spad_distribution/`
- `../canopy_photosynthesis/`

---

## 2. `compute_dgci_and_spad_distribution/`

**Goal:** Convert RGB information into **SPAD distribution** on leaves or canopy elements.

### Key functions / tasks

- **Color index computation**
  - Compute **DGCI (Dark Green Color Index)** and other RGB-based indices from point colors or images.
  - Aggregate color statistics per leaf / mesh.

- **Calibration to SPAD**
  - Read measured SPAD data (from SPAD meter) and matching DGCI values.
  - Fit regression models (e.g. linear / nonlinear) to relate DGCI (or similar indices) to SPAD.

- **SPAD prediction & mapping**
  - Apply the calibrated DGCI→SPAD relationship to all leaf pixels/points.
  - Produce SPAD values for each leaf, triangle, or spatial layer in the canopy.

- **Export & visualization**
  - Export SPAD maps / distributions (e.g. per triangle of a mesh or per leaf).


**Output:**  
SPAD  maps used as input in:

- `../canopy_photosynthesis/`

---

## 3. `C4_photosynthesis_model_parameterization/`

**Goal:** Estimate biochemical parameters of the **C4 leaf photosynthesis model** from gas-exchange measurements (typically using the MATLAB Bayesian Estimation Toolbox).

### Key functions / tasks

- **Data handling**
  - Import A–Ci / A–Q gas-exchange data from Excel files.
  - Clean and organize data by leaf, treatment, or measurement sequence.

- **C4 leaf model implementation**
  - Implement the steady-state C4 biochemical model (Rubisco, PEPCase, electron transport, bundle sheath leakage, etc.).
  - Optionally include temperature response of parameters.

- **Bayesian / optimization wrapper**
  - Define prior distributions and likelihood for model parameters.
  - Run MCMC or optimization to fit parameters to measured curves.
  - Monitor convergence / goodness-of-fit.

- **Diagnostics & export**
  - Plot fitted vs. observed gas-exchange curves.
  - Visualize posterior distributions or confidence intervals.
  - Save estimated parameter sets (e.g. `Vpmax`, `Vcmax`, etc.) to `.mat` / table files.

**Output:**  
Calibrated C4 parameter sets used by:

- `../canopy_photosynthesis/`

---

## 4. `canopy_photosynthesis/`

**Goal:** Simulate **canopy-scale C4 photosynthesis** using 3D structure, SPAD distribution, and calibrated biochemical parameters.

### Key functions / tasks

- **Input & scene setup**
  - Load 3D canopy geometry (e.g. meshes / point clouds with leaf IDs).
  - Load SPAD / nitrogen distribution from `../compute_dgci_and_spad_distribution/`.
  - Load C4 biochemical parameters from `../C4_photosynthesis_model_parameterization/`.
  - Assign per-leaf / per-element photosynthetic capacity based on SPAD or N.

- **Light & energy balance**
  - Compute local light interception (ray tracing or geometric shading).
  - Optionally compute leaf temperature via energy balance (shortwave, longwave, convection).

- **Leaf-level photosynthesis**
  - Evaluate C4 leaf photosynthesis for each leaf element given local light, temperature, and biochemical parameters.
  - Distinguish between sunlit/shaded fractions or explicit triangle-level resolution.

- **Scaling to canopy**
  - Sum leaf-level assimilation to plant / canopy level.
  - Integrate over time to obtain diurnal / daily canopy photosynthesis.

- **Analysis & plotting**
  - Compute and visualize vertical profiles of light, temperature, and photosynthesis.
  

**Output:**  
Canopy-level photosynthesis time series, spatial patterns, and figures for analysis.

---

## Workflow: how these modules connect

A typical end-to-end workflow for the repository:

1. **Pre-process point clouds**  
   - Use scripts in `point_clouds_color_correction_alignment/`.  
   - Input: raw point clouds from your imaging system.  
   - Output: clean, aligned, color-corrected point clouds.

2. **Estimate SPAD distribution**  
   - Use `compute_dgci_and_spad_distribution/`.  
   - Calibrate DGCI ↔ SPAD with measured SPAD values.  
   - Output: SPAD maps / distributions.

3. **Calibrate the C4 leaf model**  
   - Use `C4_photosynthesis_model_parameterization/`.  
   - Fit biochemical parameters from gas-exchange data.  
   - Output: C4 parameter sets (`Vpmax`, `Vcmax`, etc.).

4. **Simulate canopy photosynthesis**  
   - Use `canopy_photosynthesis/`.  
   - Combine 3D canopy, SPAD distribution, and C4 parameters.  
   - Output: canopy-scale photosynthesis and related analyses.
