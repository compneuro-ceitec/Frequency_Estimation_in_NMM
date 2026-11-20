# Fast Frequency Estimation via Bifurcation Analysis in Neural Mass Models

This repository contains the code used for the numerical experiments in the paper *"Fast Frequency Estimation via Bifurcation Analysis in Neural Mass Models."*

The repository is organized into five complementary components:

- **`eigenvalue-based_frequency_mapping/`** — MATLAB script (`NMM_Gast_interburst_freq_estim_GRID.m`) that maps oscillation frequencies inside the bursting and spiking bubbles of the `$(M_{D1}, \Delta_{v_\theta})$` plane. It implements and compares four methods: fast **equilibrium-continuation/eigenvalue (EP)** estimation, **limit-cycle continuation (LC)**, simulation with **Lomb–Scargle periodogram (LS)**, and simulation with **peak detection (FP)**.

- **`coarse_continuation/`** — MATLAB scripts that use **MatCont** to trace equilibrium and Hopf bifurcation curves for a Gast-type neural mass model (NMM). The scripts compute dominant and secondary Hopf frequencies across 2-D parameter slices, coarse 3-D grids, and target-band searches.

- **`manifold_grid_sampling/`** — Python scripts that perform a **direct symbolic/numerical grid search** over the Hopf manifold in the `(Md, Vt)` plane. The workflow produces 1-D cuts, 2-D grids, static 3-D frequency figures, and interactive 3-D visualisations.

- **`analytic_hopf_bifurcation_manifold/`** — **Maple worksheet** (`Gast_Hopf_biproduct_Delta0.mw`) and **PDF notes** (`Gast_Delta_Hopf_condition.pdf`) containing the analytic derivation of the Hopf bifurcation manifold for the Gast-type NMM in the $\Delta_{v_\theta}=0$ limit. The Maple worksheet builds the bialternate-product Hopf condition, eliminates the state variables, and produces the symbolic resultant manifold used to initialize continuation and grid-search computations.

- **`Matcont_files/`** — MATLAB figure scripts for the original and revised manuscript, MatCont system files and saved bifurcation diagrams, frequency-map data, and the Hopf-to-limit-cycle validation script added for the revision (Figure 12). See **`Matcont_files/README.md`** for installation, file placement, path configuration, and figure reproduction. MatCont 7p6 must be installed separately; the full software package is **not** included in this repository.

## Quick start

- For MatCont system and saved diagrams, figure reproduction (Figures 1–4, 8–10, 12), and MATLAB setup, see **`Matcont_files/README.md`**.
- For continuation-based frequency estimation in MATLAB, see **`coarse_continuation/README.md`**.
- For the Python grid-search and plotting pipeline, see **`manifold_grid_sampling/README.md`**.
- For the analytic Maple derivation of the Hopf manifold, open **`analytic_hopf_bifurcation_manifold/Gast_Hopf_biproduct_Delta0.mw`** in Maple (or read the accompanying PDF notes).
- For the comparison of eigenvalue-based and simulation-based frequency maps in the `$(M_{D1}, \Delta_{v_\theta})$` plane, see **`eigenvalue-based_frequency_mapping/README.md`**.

For setup instructions and script descriptions, consult the README file in each relevant subdirectory. The analytic derivation is documented in the Maple worksheet and accompanying PDF.
