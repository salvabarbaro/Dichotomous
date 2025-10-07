# Dichotomous Preferences: Concepts, Measurement, and Evidence

This repository contains the replication files accompanying the paper  
**"Dichotomous Preferences: Concepts, Measurement, and Evidence"**  
by [Salvatore Barbaro](mailto:sbarbaro@uni-mainz.de) and [Anna-Sophie Kurella](mailto:a.kurella@ipw.uni-hannover.de).

📄 Link to the paper: *[to be inserted]*

---

## Repository structure

- **RF_France22.R** – Replication file for the France 2022 dataset  
- **RF_Graz.R** – Replication file for the Graz dataset  
- **RF_Grenoble.R** – Replication file for the Grenoble dataset  
- **AuxFunctions.RData** – Contains auxiliary functions used by the replication files (reproducible via **AuxFunctions.R**)  
- **OLSNEW.R** – Runs the OLS regressions reported in Section 7 of the paper  
- **LogRegNEW.R** – Runs the logistic regressions reported in Section 8 of the paper  
- **CSES_cluster.R** – Replication file for Subsection 6.2  
- **CSESSocioDem.R** – Replication file for Section 8  
- **DATA/** – Directory containing all raw datasets required to run the replication files (freely accessible)

---

## Software environment

- Tested on **TuxedoOS (Debian/Ubuntu)** with  
  **R version 4.3.3 (2024-02-29) — "Angel Food Cake"**

### Required R packages

```r
library(dplyr)
library(cluster)
library(factoextra)
library(ggplot2)
library(ineq)
library(dineq)
library(rstatix)
library(tidyr)
library(IC2)        # install via remotes::install_version("IC2")
library(parallel)
library(modelsummary)
library(rlang)
```

---

## Reproduction instructions

1. **Clone this repository**
   ```bash
   git clone https://github.com/salvabarbaro/Dichotomous.git
   ```

2. **Open R (≥ 4.3.3)**  
   We recommend Linux (tested on TuxedoOS), as some functions rely on OS-specific parallelisation.

3. **Install the required packages**  
   See list above. For `IC2` (no longer maintained), install via:
   ```r
   remotes::install_version("IC2")
   ```

4. **Ensure the `DATA/` folder is present**  
   All replication files depend on the raw datasets in this directory.

5. **Run the replication scripts** according to the section of the paper:
   - `RF_France22.R`, `RF_Graz.R`, `RF_Grenoble.R` – dataset analyses  
   - `OLSNEW.R` – OLS regressions (Section 7)  
   - `LogRegNEW.R` – Logistic regressions (Section 8)  
   - `CSES_cluster.R` – Subsection 6.2  
   - `CSESSocioDem.R` – Section 8  

---

## Data sources

- Grenoble data: [Zenodo DOI: 10.5281/zenodo.3548303](https://doi.org/10.5281/zenodo.3548303)  
- Graz data: provided by Andreas Klamler (see Darmann and Klamler 2023, *Public Choice*)  
- France22 data: [Zenodo DOI: 10.5281/zenodo.10568798](https://doi.org/10.5281/zenodo.10568798)  
- CSES data: freely available at [cses.org](https://cses.org) (not included in `DATA/`)

---

## Notes and warnings

- ⚠️ **Parallelisation**  
  Some replication files use `mclapply()`, which only works on Linux/Unix.  
  On Windows/macOS, replace with `lapply()` or adapt the code.

- 🕒 **Computation time**  
  Running all replication scripts takes about **1.5 hours** on a standard workstation.

- 📂 **Auxiliary functions**  
  `AuxFunctions.RData` is required by the replication scripts.  
  If missing, regenerate by running `AuxFunctions.R`.

---

## Citation

If you use this code or data, please cite the paper:

*Barbaro, S., & Kurella, A.-S. (2025). Dichotomous Preferences: Concepts, Measurement, and Evidence. [link to be inserted]*

---

## License

This project is licensed under the **Creative Commons Attribution (CC-BY)** license.  
You are free to share and adapt the material, provided appropriate credit is given.

## Acknowledgements

The data from the Graz experiment were provided by **Andreas Klamler (University of Graz)** and are published in this repository with his kind permission.  
The authors are grateful to Andreas Klamler for his generous support.
