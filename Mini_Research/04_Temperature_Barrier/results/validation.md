# Numerical and artifact validation

The package is prepared for independent scientific review, not declared scientifically accepted.

| Check | Result |
|---|---|
| production_runs | 92 |
| complete_runs | 76 |
| censored_runs | 16 |
| mathematical_and_artifact_tests | 16 |
| tests_passed | True |
| convergence_cases | 56 |
| max_MR3_regression_relative_peak_difference | 3.708497731569338e-06 |
| max_strict_peak_difference | 3.478365850240195e-06 |
| max_Radau_peak_difference | 1.8205184204589386e-06 |
| max_input_peak_difference | 4.906907042245834e-06 |
| max_production_energy_residual_scaled | 4.376246103923962e-07 |
| max_independent_16point_quadrature_residual | 1.2350301686429456e-06 |
| valid_global_sensitivity_contrasts | 45 |
| all_boundary_stop_reasons | ['coefficient_upper_boundary'] |
| source_energy_zero | True |
| physical_validation | False |
| scientific_acceptance | awaiting independent web review |

The exact legacy equation limit is checked pointwise against MR3. Peak regression differences measure the new integration, not a fitted target. The 56 numerical checks contain 36 strict DOP853 runs (24 new-model principal cases and 12 M0 controls), 12 independent Radau mixture runs, and eight input-resolution tests. The 144 quadrature records use 8/16-point Gauss-Legendre rules at the peak and observed end for all strict runs.

The largest argmax-time differences in weak 1.0 bar cases correspond to adjacent nearly tied cycles. Actual maxima/times remain unmodified; the declared 1e-8 relative near-peak band is a conditioning convention, not an error bound. Boundary-time differences and the uncomputed tail remain separate. Censored maxima are never global five-cycle values.

All NASA heat-capacity evaluations used for accepted states remain in their saved domains. Independent quadrature checks the closed gas subsystem, not the energy of the whole gas-liquid system. Coefficient uncertainty, fixed transport, ideal-gas/nonreacting composition and KM inadequacy are not reduced by tighter integration.

The 10-page English PDF has seven pages of main text, one reference page and two appendix pages. Equations, all tables, five JPG figures, captions and every rendered page were visually inspected; no overfull boxes, unresolved references or clipped content remained. The scientific report contains no commit SHA or assistant workflow instructions.

A separate clean copy made from the accepted Git archive plus the MR4 publication file list passed all 16 tests. Four bounded fresh-process runs reproduced the saved temperature maxima exactly at floating-point precision and retained both complete/censored statuses. All 611 upstream hashes matched. `clean_reproduction.csv` and `clean_snapshot_verification.json` retain the outcomes and hashes of the tested summary, coefficients and PDF. The verification journal itself was added after that clean calculation.
