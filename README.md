# Graphene-Nanoribbon-NEGF

MATLAB scripts for simulating electron transport in graphene nanoribbons using the Non-Equilibrium Green's Function (NEGF) framework.

## Features

- NEGF-based transport workflow for graphene nanoribbon models
- MATLAB implementations for Hamiltonian construction and transport analysis
- Example scripts for comparative studies and plotting
- Repository structure suitable for further extension and reproducible studies

## Requirements

- MATLAB (recommended: a recent version)
- Required MATLAB toolboxes/dependencies: **to be specified by the repository owner**

## Installation

1. Clone or download this repository.
2. Open MATLAB.
3. Set the working directory to the repository root:
   - `/home/runner/work/Graphene-Nanoribbon-NEGF/Graphene-Nanoribbon-NEGF`
4. Add the repository path in MATLAB if needed:
   - `addpath(genpath(pwd))`

## Usage

> The exact run sequence may vary based on your workflow.  
> Replace placeholders below with your confirmed script order and parameters.

1. Identify the main analysis script(s) to run.
2. Configure simulation parameters (e.g., geometry, bias/energy grid, temperature) in the relevant `.m` files.
3. Run the selected script(s) from MATLAB.
4. Review generated figures/data and adjust settings as needed.

Suggested placeholders to fill later:

- **Primary entry script(s):** `<main_script_name.m>`
- **Optional preprocessing script(s):** `<preprocess_script_name.m>`
- **Postprocessing/plotting script(s):** `<plot_script_name.m>`

## Repository Structure

Current repository contents are MATLAB scripts, including:

- `AGNR_Hamiltonian.m`
- `ZGNR_Hamiltonian.m`
- `GNR_Curren.m`
- `COMPARE_BANDSTRUCTURES.m`
- `Compare_GNR_IV_publication_fast.m`
- `H_7AGNR.m`
- `H_11AGNR.m`
- `H_8ZGNR.m`

> Add or update this section as the project structure evolves.

## Outputs

Expected outputs are project-dependent and should be confirmed by the repository owner.  
Typical outputs for NEGF transport studies may include:

- Transmission-related data/curves
- Current–voltage (I–V) characteristics
- Band structure plots
- Saved MATLAB figures and/or `.mat` data files

## References / Citation

If this repository supports a publication, add the citation details here.

Example placeholder:

```text
Author(s), "Title," Journal/Conference, Year, DOI.
```

If you use this code in academic work, please cite the corresponding paper/project documentation once available.
