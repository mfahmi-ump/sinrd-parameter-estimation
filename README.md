# SI4RD Parameter Estimation

Bayesian (MCMC) parameter estimation for an **SI4RD compartmental epidemic model** in MATLAB. The model tracks Susceptible (S), four infectious classes (I1-I4), Recovered (R), and Deceased (D) populations, and is fitted to Malaysian case data using adaptive MCMC on automatically detected time segments.

## What it does

1. **Loads** daily compartment data from `si4rd_data.csv`.
2. **Detects change points** in the epidemic curve from peaks in the smoothed second derivative of I1 and I3, splitting the series into segments.
3. **Builds starting values** for all 18 parameters from the data in each segment.
4. **Refines** them with least squares (`fminsearch`).
5. **Samples the posterior** with DRAM (Delayed Rejection Adaptive Metropolis) using Marko Laine's `mcmcrun`.
6. **Simulates** the model (RK4) with the posterior means and plots it against the observed data.
7. **Stores** the time-varying parameter estimates in a matrix (`paramestimate`, 18 x T).

## Estimated parameters (18 per segment)

| Parameters | Meaning |
|---|---|
| `p1`-`p4` | Share of `pvax` carried by each infectious class I1-I4 |
| `beta1`-`beta4` | Transmission rates |
| `r1`-`r4` | Recovery rates |
| `d1`-`d4` | Death rates |
| `gamma` | Small rate parameter (bounded above by 1e-5) |
| `pvax` | Fraction derived from `ES ./ S` in the data |

The model equations are implemented in `si4rdsrk4.m` and the fitting objective in `si4rdss.m`.

## Repository structure

```
.
├── si4rd_parameter_estimation.m   # main script
├── si4rdss.m                      # sum-of-squares function (objective)
├── si4rdsrk4.m                    # RK4 solver for the SI4RD model
├── si4rd_data.csv                 # input data
└── README.md
```

## Requirements

- MATLAB R2016b or later
- Signal Processing Toolbox (`findpeaks`)
- [MCMC toolbox](https://mjlaine.github.io/mcmcstat/) by Marko Laine (`mcmcrun`, `chainstats`)
- [Panel](https://www.mathworks.com/matlabcentral/fileexchange/20003-panel) by Ben Mitch (for the plots)

Add these to the MATLAB path before running.

## Input data

`si4rd_data.csv` is read with `readmatrix`. The script uses these columns:

| Column | Content |
|---|---|
| 1 | Time (day index) |
| 2 | I1 |
| 3 | R |
| 4 | D |
| 5 | I2 |
| 6 | I3 |
| 7 | I4 |
| 8 | S (as a fraction of the total population) |
| 11 | ES (used with S to derive `pvax`) |

Columns 9-10 are not used by this script. Population size is set in the script (`Ptotal = 3.153e7`) and the first row corresponds to 25 February 2021 (`day0 = 24 Feb 2021`).

## Usage

```matlab
% from the repository folder, with dependencies on the path
si4rd_parameter_estimation
```

Key settings you may want to change at the top of the script:

| Setting | Default | Purpose |
|---|---|---|
| `Ptotal` | `3.153e7` | Total population |
| `day0` | `2021-02-24` | Date offset for plot axes |
| `smoothWin1`, `smoothWin2` | `32`, `20` | Smoothing windows for change-point detection |
| `options.nsimu` | `5000` | MCMC samples per run (three runs are chained) |

## Output

- A 7-panel figure (S, I1-I4, R, D) showing observed data against the fitted model for each segment.
- `paramestimate`: an 18 x T matrix of posterior-mean parameters per day.

## Notes

- Zeros in the infectious series are replaced by a tiny random positive value (`0.01*rand()`), so results vary slightly between runs unless you set the random seed (`rng(1)`).
- Segments are fitted from the end of the series backwards.

## License

MIT (see `LICENSE`).

## Citation

If you use this code, please cite the associated paper: *[add reference / DOI here]*.
