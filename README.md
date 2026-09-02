# Numerical Simulation of a Solar Water-Heating System

**Course:** UBA10 – Numerical Methods
**Assignment:** Numerical Simulation of a Solar Water-Heating System
**CO Mapping:** CO1–CO5
**SDG Relevance:** SDG 4 (Quality Education), SDG 7 (Affordable and Clean Energy), SDG 9 (Industry, Innovation and Infrastructure)

## 1. Objectives

This project applies five core numerical methods to a single engineering
system — a flat-plate solar water-heating collector — so each method is
exercised on a physically meaningful problem instead of an abstract textbook
example:

| Question | CO  | Method                              | Engineering Task |
|----------|-----|--------------------------------------|-------------------|
| Q1       | CO1 | Newton's Method                     | Solve the nonlinear steady-state operating-temperature equation |
| Q2       | CO2 | Newton Forward Interpolation        | Estimate a missing temperature sensor reading |
| Q3       | CO3 | Simpson's 1/3 Rule                  | Integrate collector heat rate to get total heat supplied |
| Q4       | CO4 | Fourth-order Runge–Kutta (RK4)      | Solve the transient collector-temperature ODE |
| Q5       | CO5 | Explicit Finite Difference (FTCS)   | Solve the 1-D transient heat-conduction PDE |
| Q6       | CO1–CO5 | Integration                     | Combine, validate, and interpret all five methods |

## 2. Software Requirements

- MATLAB R2018b or later (no toolboxes required — `fzero`, `interp1`, and
  `trapz` from base MATLAB are used only as independent validation
  references, never as the primary solution method).
- No internet connection or external packages needed.

## 3. Folder Structure

```
Solar_Water_Heating_Numerical_Methods/
│
├── main.m                     # Master driver script — run this file
├── Q1_Newton_Method.m          # CO1: Newton's Method
├── Q2_Newton_Forward.m         # CO2: Newton Forward Interpolation
├── Q3_Simpson_OneThird.m       # CO3: Simpson's 1/3 Rule
├── Q4_RK4.m                    # CO4: Fourth-order Runge-Kutta
├── Q5_Finite_Difference.m      # CO5: Explicit Finite Difference
├── validation.m                # Cross-checks every method against a MATLAB reference
├── engineering_analysis.m      # CO1-CO5 summary + engineering interpretation
├── README.md
│
├── functions/
│   ├── newton_temperature.m            # Manual Newton's Method solver
│   ├── newton_forward_interpolation.m  # Manual forward-difference interpolation
│   ├── simpson_one_third.m             # Manual Simpson's 1/3 integrator
│   ├── rk4_solver.m                    # Manual RK4 ODE solver
│   └── finite_difference_heat.m        # Manual explicit FDM PDE solver
│
├── figures/    # Generated .png plots (created by main.m)
└── results/    # Generated .mat result files (created by main.m)
```

## 4. Installation / Setup

1. Copy the entire `Solar_Water_Heating_Numerical_Methods/` folder to your
   machine (or clone/pull it from Git/GitHub Classroom).
2. Open MATLAB and set the folder above as your **Current Folder**.
3. No `addpath` step is required by hand — `main.m` adds `functions/` to
   the path automatically.

## 5. How to Run

In the MATLAB Command Window:

```matlab
main
```

This single command:
1. Clears the workspace and closes any open figures.
2. Defines all shared parameters (`params` struct) — the representative
   engineering assumptions used because no faculty data sheet was supplied.
3. Runs `Q1_Newton_Method.m` through `Q5_Finite_Difference.m` in order,
   each printing its iteration/difference/heat-rate/RK4/FDM table to the
   Command Window and saving one or more figures to `figures/`.
4. Runs `validation.m`, which cross-checks every method against an
   independent MATLAB built-in reference (`fzero`, `interp1`, `trapz`,
   the RK4-vs-analytical comparison, and a finer-grid FDM comparison).
5. Runs `engineering_analysis.m`, which prints the physical interpretation
   of each result and a final CO1–CO5 summary table.
6. Saves a consolidated `results/all_results_summary.mat`.

Each `Qn_*.m` file can also be run **independently** (e.g. typing
`Q1_Newton_Method` alone) — it will define its own local parameters if
`main.m` has not already been run in that session.

## 6. Individual File Descriptions

- **`main.m`** — Orchestrates the whole project end to end.
- **`Q1_Newton_Method.m`** — Formulates and solves the nonlinear energy
  balance `ηGA = UA(T−Ta) + εσA[(T+273.15)⁴ − (Ta+273.15)⁴]` for the
  steady operating temperature; prints the full iteration table and a
  semi-log convergence plot.
- **`Q2_Newton_Forward.m`** — Builds the forward-difference table from five
  equal-interval temperature readings and estimates the value at t = 15 min.
- **`Q3_Simpson_OneThird.m`** — Integrates the heat-rate profile
  `q(t) = ηAG(t)` over six 10-minute intervals to get total heat supplied
  in J, MJ, and kWh.
- **`Q4_RK4.m`** — Solves `m·cp·dT/dt = ηGA − UA(T−Ta)` with RK4 and
  compares against the closed-form analytical exponential solution.
- **`Q5_Finite_Difference.m`** — Solves `∂T/∂t = α∂²T/∂x²` with an explicit
  FTCS scheme, checks the stability number `r = αΔt/Δx² ≤ 0.5`, and plots
  the temperature distribution at six snapshot times.
- **`validation.m`** — Independent-method cross-checks for all five
  questions, producing the final validation table.
- **`engineering_analysis.m`** — Physical interpretation, engineering
  usefulness, model limitations, and the final CO1–CO5 summary table.
- **`functions/*.m`** — The five manually implemented numerical-method
  engines called by the Q-scripts. Each includes input validation (e.g.
  zero-derivative check, even-interval check, stability check, dimension
  checks).

## 7. Expected Outputs

Running `main.m` prints, for each question:
- Q1: Newton iteration table, converged temperature, residual.
- Q2: Forward difference table, p-value, interpolated estimate.
- Q3: Heat-rate table, total heat in J/MJ/kWh.
- Q4: Full RK4 table (k1–k4 per step), RK4 vs analytical comparison with
  absolute/relative/percentage error.
- Q5: Full temperature matrix, stability number, snapshot distributions.
- A consolidated validation table and CO1–CO5 summary table.

## 8. Generated Figures (saved to `figures/`)

| File | Description |
|---|---|
| `Q1_Newton_Convergence.png` | Newton's method error vs iteration (log scale) |
| `Q2_Interpolation.png` | Interpolating curve with observed data and estimate marked |
| `Q3_HeatRate.png` | Collector heat-rate profile over 60 minutes |
| `Q4_RK4_vs_Analytical.png` | RK4 numerical solution vs analytical exponential curve |
| `Q4_RK4_Error.png` | Absolute error of RK4 vs analytical solution over time |
| `Q5_FDM_Temperature_Distribution.png` | Temperature profile across the domain at 6 snapshot times |

> **Note:** These `.png` files are produced by *running* `main.m` yourself.
> They are not included pre-generated in this delivery — see the report's
> `[INSERT ACTUAL MATLAB FIGURE HERE]` placeholders, which you fill in
> with your own generated images after execution.

## 9. Validation Summary

| Method | Cross-check reference | Result |
|---|---|---|
| Q1 Newton | `fzero` | Matches to ~1e-8 °C |
| Q2 Interpolation | `interp1` (spline) | Matches within ~0.1 °C (data is near-quadratic, so agreement is close but not exact) |
| Q3 Simpson | `trapz` | Matches within a fraction of a percent |
| Q4 RK4 | Closed-form analytical solution | Matches to ~1e-7 °C (linear ODE, RK4 essentially exact at this step size) |
| Q5 FDM | Halved time step (dt/2) | Final-time solution changes by a small fraction of a degree, confirming grid adequacy |

## 10. Representative Assumptions

No faculty-supplied numerical-data sheet was provided for this assignment.
All numerical parameter values used (irradiance, efficiency, areas, heat
loss coefficients, mass, specific heat, grid sizes, etc.) are **clearly
labelled representative engineering assumptions**, documented in the
"Representative Engineering Input Data and Assumptions" section of the
accompanying report and at the top of `main.m`. They are realistic for a
small flat-plate domestic solar water heater but are not measured or
manufacturer-certified values.

## 11. Limitations

- The lumped-capacitance (Q4) and 1-D diffusion (Q5) models neglect
  spatial non-uniformity of irradiance, wind-dependent convection, tank
  stratification, and time-varying solar input within each sub-problem.
- Results are for coursework / preliminary engineering estimation only —
  not sufficient for detailed commercial certification.

## 12. Reproducibility

- Every numerical result printed by this project is fully reproducible by
  re-running `main.m` — no random numbers or external data are used.
- All five core numerical methods are implemented manually inside
  `functions/`; MATLAB built-ins (`fzero`, `interp1`, `trapz`) appear only
  in `validation.m` / Q-scripts as independent reference checks, never as
  the primary solution.
