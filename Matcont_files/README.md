# MatCont files and figures for the NMM manuscript

This directory contains the MatCont system definition, saved continuation diagrams, MATLAB figure scripts, and precomputed frequency-map data used in *Fast Frequency Estimation via Bifurcation Analysis in Neural Mass Models*. The files support the figures in the revised manuscript, including the Hopf-to-limit-cycle frequency comparison added during revision.

> **Important:** To use the saved continuation diagrams in MatCont, install **MatCont 7p6** separately from [SourceForge](https://sourceforge.net/projects/matcont/). MatCont is **not** included in this repository. Each user must adjust the local MatCont path in the MATLAB scripts.

## Requirements

- MATLAB (the version used for the final computations should be recorded separately if exact environment reproduction is needed).
- MatCont 7p6, installed locally (the installation folder is referred to below as `MatCont7p6`).
- The MATLAB scripts and `.mat` data files provided here.

## Files in this directory

```text
Matcont_files/
├── Dopamine_NMM_fast4.m             # MatCont system definition
├── Dopamine_NMM_fast4.mat           # MatCont system settings
├── diagram_DeltaMd/                 # Saved equilibria, Hopf curves, LPCs, and limit cycles
├── script_NMM_figures.m             # Original figure-generation script
├── script_NMM_figures_revised.m     # Revised figure-generation script
├── R2_8_Hopf_LC_validation.m        # Hopf-to-limit-cycle frequency validation (Fig. 12)
├── freq_data_EP_20.mat              # Equilibrium eigenvalue-based estimates
├── freq_data_FP_20.mat              # Peak-based frequency estimates
├── freq_data_LS_20.mat              # Lomb–Scargle frequency estimates
├── freq_data_LC_20.mat              # Limit-cycle frequency data
└── README.md
```

The `diagram_DeltaMd/` directory contains the saved MatCont bifurcation and orbit data. Keep its filenames unchanged.

## MatCont folder layout (expected by the plotting scripts)

The original `script_NMM_figures.m` loads continuation data using the following layout **inside the local MatCont installation**:

```text
MatCont7p6/
└── Systems/
    ├── Dopamine_NMM_fast4.m
    ├── Dopamine_NMM_fast4.mat
    └── Dopamine_NMM_fast4/
        └── diagram_DeltaMd/
            ├── H_DeltaMd_1a.mat
            ├── H_DeltaMd_1b.mat
            ├── LC_Md_Delta_0.5_left_dense.mat
            ├── P_bursting_Md_0.668.mat
            ├── P_spiking_Md_5.mat
            └── ...
```

The directory name `MatCont7p6` is illustrative; the installation may be anywhere on the computer. The path above reflects what the supplied plotting code expects, rather than a requirement to install MatCont on a particular disk.

## Setup

1. **Install MatCont 7p6** from [SourceForge](https://sourceforge.net/projects/matcont/) and make its launcher accessible in MATLAB. You can start the graphical interface with:

   ```matlab
   matcont
   ```

2. **Copy the MatCont system files** from this directory to the `Systems` directory of your local MatCont installation:

   ```text
   MatCont7p6/Systems/Dopamine_NMM_fast4.m
   MatCont7p6/Systems/Dopamine_NMM_fast4.mat
   ```

3. **Copy the entire `diagram_DeltaMd` folder** into:

   ```text
   MatCont7p6/Systems/Dopamine_NMM_fast4/diagram_DeltaMd/
   ```

   Preserve its subdirectory structure and the saved `.mat` files. These are precomputed MatCont results; regenerating the continuations is not necessary merely to redraw the figures.

4. **Set the MatCont path in the plotting script(s).** At the start of `script_NMM_figures.m` (and in the revised script if it uses the same convention), edit the variable `Path` to point to the `Systems` directory **on your computer**. For example:

   ```matlab
   Path = 'C:\prace\MatCont7p6\Systems';
   ```

   This is intentionally a user-specific absolute path. Replace the example path with your own. The scripts then find stored diagrams under `[Path, '\Dopamine_NMM_fast4\diagram_DeltaMd\...']`.

5. **Keep the frequency-data files accessible.** The original plotting script loads `freq_data_EP_20.mat`, `freq_data_FP_20.mat`, `freq_data_LS_20.mat`, and `freq_data_LC_20.mat` by filename from MATLAB's current working directory. The simplest approach is to keep these files beside the plotting scripts, and set MATLAB's **Current Folder** to this `Matcont_files/` directory before executing the relevant sections.

## Reproducing the manuscript figures

Open `script_NMM_figures_revised.m` in MATLAB, edit the local MatCont path, and execute the required figure sections (including any preceding setup section). The original `script_NMM_figures.m` is retained for comparison with the figures before revision.

| Manuscript figure | Code and input data |
| --- | --- |
| **Figure 1 – collective bursting** | Figure section in `script_NMM_figures_revised.m`; saved orbit `diagram_DeltaMd/P_bursting_Md_0.668.mat`. |
| **Figure 2 – collective oscillations** | Figure section in `script_NMM_figures_revised.m`; saved orbit `diagram_DeltaMd/P_spiking_Md_5.mat`. |
| **Figure 3 – FP, LS, EP maps** | Revised figure script, using `freq_data_FP_20.mat`, `freq_data_LS_20.mat`, `freq_data_EP_20.mat`, and saved bifurcation curves. |
| **Figure 4 – method differences and one-parameter cuts** | Revised figure script, using the frequency-map files and saved curves. The revised relative differences exclude points with a zero EP denominator. `freq_data_LC_20.mat` is also supplied for comparisons using continued cycles. |
| **Figure 8 – two-parameter bifurcation diagrams** | Revised figure script; saved Hopf and LPC continuations in `diagram_DeltaMd/`. |
| **Figure 9 – bursting limit-cycle diagrams** | Revised figure script; saved left-branch equilibrium and limit-cycle continuations in `diagram_DeltaMd/`. |
| **Figure 10 – oscillatory limit-cycle diagrams** | Revised figure script; saved right-branch equilibrium and limit-cycle continuations in `diagram_DeltaMd/`. |
| **Figure 12 – local Hopf-frequency validation** | `R2_8_Hopf_LC_validation.m`; comparison of critical Hopf eigenfrequencies with frequencies from continued periodic orbits at `Delta = 0.5`. The required saved limit-cycle diagrams must be available. |

**Figures 5, 6, and 11** are covered by the separate `coarse_continuation/` workflow at the repository root. **Figure 7** is covered by `manifold_grid_sampling/`. Refer to their own README files for setup and scripts.

## Notes

- The original plotting script uses **MATLAB code sections** (`%% FIG. ...`). Execute the setup section first and then the figure sections you need; there is no need to recompute all numerical continuations for plotting.
- The original script contains a call to `bluewhitered(256)` for a diverging colormap. This function is not part of every standard MATLAB installation; if it is absent, supply the corresponding helper or use an appropriate available colormap. Check whether this dependency is still present in the revised script.
- Model time is expressed in milliseconds. The frequency conversions used in the manuscript are `f_H = 1000*abs(imag(lambda_H))/(2*pi)` and `f_LC = 1000/T`, with `T` in milliseconds.
- The validation near the Hopf points assesses **local onset-frequency accuracy**; it does not establish accuracy of the equilibrium-based estimate for arbitrary finite-amplitude oscillations.
- The original plotting script and saved diagram files were checked for the paths and filenames above. The newly revised plotting and validation scripts were supplied by the authors for this repository update; if their final versions use additional paths or data files, document those dependencies here before release.