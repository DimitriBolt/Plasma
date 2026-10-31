# Mini-Research 4: independent scientific review handoff

The published dispatch message supplies the full SHA and immutable GitHub/PDF links for the commit containing this file. In a checkout, its identity is obtained with `git rev-parse HEAD`. The scientific base was `1a12fac295dc904569ab2a57ce4f68a8ea75f0c4`; MR1--MR3 and all 611 tracked inputs remain byte-identical. Only Mini-Research 4 is implemented.

Report: `Mini_Research/04_Temperature_Barrier/report/04_Temperature_Barrier.pdf` (10 pages, seven main-text pages). Source, all runs, trajectories, parameters, tests and source records are in the same package. The report is by Dimitri Bolt, academic supervisor Prof. Dr. Ildar R. Gabitov.

## Main MR3/MR4 comparison, eta = 1

OW = one-way on the same saved mechanics used by MR3; CO = fully coupled. M1 uses variable dry-gas capacity with the legacy constant vapor-pressure convention. M2 uses a fixed initially saturated H2O inventory in pressure and energy, with no separate constant pv. A censored observed maximum is **not** a global five-cycle Tmax.

| Pa (bar) | Gas | Mode | MR3 Tmax (K) | M1 observed (K) | M2 observed (K) | M2 scope | M2 global delta vs MR3 |
|---:|---|---|---:|---:|---:|---|---:|
| 1.00 | air | OW | 441.26 | 439.84 | 440.81 | five-cycle global | -0.10% |
| 1.00 | air | CO | 707.13 | 693.36 | 552.87 | five-cycle global | -21.82% |
| 1.00 | argon | OW | 575.96 | 575.96 | 574.00 | five-cycle global | -0.34% |
| 1.00 | argon | CO | 886.82 | 886.82 | 661.45 | five-cycle global | -25.41% |
| 1.15 | air | OW | 3898.73 | 2667.17 | 2647.05 | five-cycle global | -32.10% |
| 1.15 | air | CO | 4423.44 | 3416.99 | 2798.39 | five-cycle global | -36.74% |
| 1.15 | argon | OW | 22142.15 | 22142.11 | 6000.00 | censored | undefined |
| 1.15 | argon | CO | 7893.02 | 7893.02 | 5674.86 | five-cycle global | -28.10% |
| 1.30 | air | OW | 15174.00 | 6777.47 | 6000.00 | censored | undefined |
| 1.30 | air | CO | 10555.02 | 6741.21 | 6000.00 | censored | undefined |
| 1.30 | argon | OW | 213181.59 | 213180.80 | 6000.00 | censored | undefined |
| 1.30 | argon | CO | 26421.63 | 26421.63 | 6000.00 | censored | undefined |

## Energy budgets at the M2 observed thermal maximum

All energy columns below are nJ. Work is by the gas; heat loss is signed. The source and chemical-energy change are exactly zero in these frozen-mixture ODEs. The residual identity is delta-U + W-by + Q-loss = 0. The residual column is the maximum scaled residual over the whole observed interval.

| Pa (bar) | Gas | Mode | delta-U | W-by | signed Q-loss | max scaled residual |
|---:|---|---|---:|---:|---:|---:|
| 1.00 | air | OW | 0.145338 | -0.181059 | 0.035721 | 8.721e-10 |
| 1.00 | air | CO | 0.258552 | -1.154035 | 0.895484 | 2.012e-11 |
| 1.00 | argon | OW | 0.166652 | -0.256323 | 0.089671 | 1.707e-10 |
| 1.00 | argon | CO | 0.218683 | -0.219551 | 0.000868 | 2.134e-12 |
| 1.15 | air | OW | 2.876257 | -2.775492 | -0.100766 | 4.128e-07 |
| 1.15 | air | CO | 3.083595 | -2.958436 | -0.125159 | 1.251e-07 |
| 1.15 | argon | OW | 3.482130 | -3.459308 | -0.022822 | 5.967e-09 |
| 1.15 | argon | CO | 3.281382 | -4.705742 | 1.424360 | 7.269e-08 |
| 1.30 | air | OW | 7.664122 | -7.461882 | -0.202240 | 8.353e-09 |
| 1.30 | air | CO | 7.664122 | -7.452968 | -0.211154 | 5.742e-09 |
| 1.30 | argon | OW | 3.482130 | -3.365093 | -0.117037 | 4.177e-09 |
| 1.30 | argon | CO | 3.482130 | -3.359789 | -0.122341 | 3.511e-09 |

End-state budgets, gross compression input, the last expansion-maximum-to-thermal-peak budget, signed heat histories and raw residuals are saved separately. This is a closed gas-subsystem balance, not conservation of gas-liquid total energy.

## Numerical evidence

- production_runs: 92.
- complete_runs: 76.
- censored_runs: 16.
- mathematical_and_artifact_tests: 16.
- tests_passed: True.
- convergence_cases: 56.
- max_MR3_regression_relative_peak_difference: 3.708497731569338e-06.
- max_strict_peak_difference: 3.478365850240195e-06.
- max_Radau_peak_difference: 1.8205184204589386e-06.
- max_input_peak_difference: 4.906907042245834e-06.
- max_production_energy_residual_scaled: 4.376246103923962e-07.
- max_independent_16point_quadrature_residual: 1.2350301686429456e-06.
- valid_global_sensitivity_contrasts: 45.
- source_energy_zero: True.

Weak 1.0 bar argmax times can move by one period when successive peaks nearly tie. `peak_time_conditioning.csv` retains this explicitly; the actual selected maxima/times are not averaged or silently relabeled. The near-peak band is 1e-8 relative, not an error bound. All 16 incomplete runs stop at the H2O coefficient upper boundary, 6000 K, during their observed interval. No table boundary is called a physical temperature cap.

## Implemented and assessed mechanisms

Implemented: NASA9 N2/O2/H2O sensible energy with continuous energy joins; positive variable capacity; fixed 79/21 N2/O2 or ground-state Ar; fixed vapor in both EOS and energy; exact full-pressure KM derivative and viscous acceleration term; signed MR3 heat exchange; separate one-way/coupled histories; pressure mismatch; refined thermal/radius events; energy and applicability diagnostics.

Controls: exact MR3 equation limit; same-component capacities frozen at 293.15 K; saturated and enhanced trapped inventories; dry mixture pressure-convention control; eta=0.5,1,2 at 1.3 bar and complete air comparisons at 1.15 bar. Vapor-fraction contrasts hold initial total pressure/moles fixed and displace dry gas. The vapor-free mixture control consequently contains more dry gas than legacy M1. In wet Ar, variable/frozen capacity isolates H2O excitation because Ar capacity is unchanged.

Energy assessment only: N2/O2 dissociation, first-step/full H2O dissociation, and first ionization of Ar/H/O/N. NIST reaction costs at 298.15 K and isolated-atom thresholds are compared to frozen thermal energy and compression input. At the 1.3 bar saturated coupled boundary, first/full water dissociation demand 0.431/0.801 nJ, N2 dissociation 34.074 nJ, and complete first Ar ionization 69.559 nJ. These are separate hypothetical demands, not inferred reaction fractions, kinetic rates, equilibrium abundances or new chemically corrected temperatures. Overlapping dissociation channels must not be added.

## Physical and data limits

N2/O2 NASA data: 200--20000 K; H2O: 200--6000 K. Trial-stage polynomial continuation is used only to bracket the first terminal event; accepted output is not extrapolated. Analytic ground-state Ar has constant capacity at all positive T mathematically, while electronic excitation/ionization are absent. Frozen-capacity controls are separate formal calorically perfect gases, not NASA polynomial extrapolations.

Enhanced vapor fractions (5/10 times the initially saturated fraction) are mechanically balanced hypothetical trapped states, supersaturated at 293.15 K; they do not model capture or evaporation. Evaporation/condensation, latent heat, chemical kinetics/equilibrium, spatial temperature gradients, nonideal gas/liquid EOS, shocks, radiation and shape instability remain unresolved. Gas density, pressure, wall Mach, liquid pressure/bulk-modulus proxy and relaxation/collapse time diagnostics are saved. Conductivities remain the MR3 300 K reference constants for controlled comparison, not validated high-temperature mixture transport. Numerical convergence is not physical validation. No extreme regime is declared impossible and no finite scan proves an absolute maximum.

## Independent review request

Please independently audit the thermodynamic energy definition and interval joins; fixed vapor accounting/static balance; full KM pressure derivative; exact MR3 limit; same-composition controls; one-way mismatch and mechanical feedback; censored maxima and thermal-event conditioning; signed/integrated energy budgets and quadrature; coefficient ranges; chemical inventory energies and the limits of energy-only inference. Please identify necessary corrections before scientific acceptance. The four conclusion categories remain numerically reproducible, conditional model-based estimate, physically supported mechanism, and requires a more complete physical model. Mini-Research 4 is not declared accepted. No MR5 work has been started.

Additional clean-copy verification: all 16 tests passed, four fresh-process control runs reproduced the saved maxima exactly at floating-point precision, all completion/censoring statuses matched, and all 611 upstream input hashes matched. See `clean_reproduction.csv` and `clean_snapshot_verification.json`.
