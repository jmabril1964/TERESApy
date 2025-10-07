# TERESApy: Python Package for ²¹⁰Pb-based Sediment Dating

**Author:** J.M. Abril-Hernández  
**Affiliation:** University of Seville, Department of Applied Physics I  
**Contact:** jmabril@us.es  
**ORCID:** [0000-0003-2540-5576](https://orcid.org/0000-0003-2540-5576)

## Overview

TERESApy is a Python implementation of the TERESA model (Time Estimates from Random Entries of Sediments and Activities) for radiometric dating of recent sediments using unsupported ²¹⁰Pb. The model generates large sets of solvers with normally distributed initial activity concentrations and sedimentation rates, evaluates their fit to empirical data, and identifies the best chronology through minimization of the adjusted chi-square function.

## Features

- Solver generation with normal distributions
- Chronology estimation and sedimentation rate profiles
- Confidence intervals (68.3%) for model parameters
- Optional inclusion of time marks as attractors
- Modular code structure for sequential execution
- Output files compatible with Gnuplot and other plotting tools

## Included Scripts

- `Random_generator.py`: Generates canonical normal distributions and random pairings
- `TERESA_map.py`: Maps the chi-square function across solver space
- `TERESA_clouds.py`: Extracts confidence regions
- `TERESA_clouds_ages.py`: Focused confidence region for chronology
- `TERESA_plots.py`: Generates theoretical profiles and chronology
- `TERESA_plot_ages.py`: Focused chronology plots
- `Core_C1.txt`: Example dataset
- `aleat_N.txt`: Randomized distribution library

## Requirements

- Python ≥ 3.6
- NumPy
- Gnuplot (optional for plotting)

## Usage

1. Place all scripts in a working folder.
2. Run `Random_generator.py` to create the distribution library.
3. Prepare your empirical data file (e.g., `Core_C1.txt`) with mass depth, ²¹⁰Pbexc, and uncertainty.
4. Execute `TERESA_map.py` with initial parameter estimates.
5. Use `TERESA_clouds.py` and `TERESA_plots.py` for full model analysis.
6. For focused chronology, use `TERESA_clouds_ages.py` and `TERESA_plot_ages.py`.

## Citation

If you use TERESApy in your research, please cite:

> Abril-Hernández, J.M. (2025). TERESApy: Python Software Package for 210Pb-based Radiometric Dating of Recent Sediments. Zenodo. DOI: [to be assigned]

## License

This software is distributed under the MIT License. See the `LICENSE` file for details.

## Acknowledgements

This package builds upon the theoretical framework developed in:
- Abril & Brunskill (2014), *Journal of Paleolimnology*
- Abril (2016, 2020, 2023, 2025), *Journal of Environmental Radioactivity* and *Quaternary Geochronology*

