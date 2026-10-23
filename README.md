# local_wccm — barotropic Lagrangian filter (MATLAB)

Public MATLAB working copy of the **barotropic Lagrangian data-assimilation** solver used during the MSc.

`lagrangianfilter_baro.m` integrates Lagrangian tracers and a filter on a barotropic spectral model (`q` / streamfunction / vorticity). Drivers generate synthetic barotropic observations, coarsen vorticity, run the filter, and plot omnidirectional spectra.

This is research code, not a GUI. The browser lab is [qg-lada-lab](https://github.com/ahmaddroobi99/qg-lada-lab). Private thesis sources stay private.

## How to run

Needs MATLAB (string arrays, so R2017a or newer). `server_run.m` and `analysis_run.m` also call `parpool`, which needs the Parallel Computing Toolbox. Not run in CI. There is no test suite.

```matlab
addpaths                 % puts this folder and drivers/ on the path
driver_baro_obs          % reads data/<case>/diagnostics.mat and t.bin
driver_coarsenqh         % coarsens q for every case in loadcases()
server_run               % runs the filter (opens a parallel pool)
```

Case names come from `drivers/loadcases.m`. The `data/` folder they point to is not in this repository; it has to be supplied. Root-level `.m` files are working copies; several differ from their `drivers/` namesakes (`analysis.m`, `coarsenqh.m`, `driver_baro_obs.m`, `generate_baro_obs.m`, `lagrangianfilter_baro.m`, `runfilter.m`), and when you run from the repo root MATLAB uses the root copy, since the current folder takes precedence over the path.

---
MSc-period work (University of Calgary, 2024). Not the public lab — see [qg-lada-lab](https://github.com/ahmaddroobi99/qg-lada-lab) for the public QG DA demo.

