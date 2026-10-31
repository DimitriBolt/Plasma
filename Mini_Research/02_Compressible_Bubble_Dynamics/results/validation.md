# Numerical and physical validation

Reference: accepted RP commit d59bc0d0119b514152f763e2537a0477c7c63663. The model file SHA256 is checked by an automated test. Only the two authorized MR1 documentation files change. The old MR1 checksum list remains the historical snapshot.

## Coverage and source checks

All five named primary papers and Part B of Plasma-plan.tex were read. Equation convention and full-source locations are documented in source_audit.md. The printed minus prefactor in Sadighi-Bonabi is not combined with the adopted plus-prefactor KM. Source hashes, solver/library versions and settings are saved in run_metadata.json.

56 amplitudes, 112 total RP/KM runs, 7 incomplete runs. All 56 KM runs complete. No Mach filter, pressure clipping, velocity clipping or denominator regularization is used. All incomplete records retain their interval, state and reason; their full-interval comparison differences are undefined and their plotted lines have gaps.

## Automated and independent checks

Mandatory tests cover static gas balance, forcing phase, dimensional base-unit rescaling, the viscous pressure derivative along an independent quadratic path, the implicit residual, nonzero KM initial acceleration, algebraic RP limit, unchanged reference hash, positive radius/domain guards, event zero crossings, table/trajectory consistency, saved-parameter reproduction, censoring and cycle correspondence. The accepted MR1 tests also pass independently (14 tests).

Independent clean environment: Python 3.12.14; package versions and differing numerical settings are in clean_environment.json. All six original amplitudes are repeated for both equations. Separate four-configuration convergence uses DOP853/Radau at 1.0, 1.10, 1.15 and 1.30 bar. Convergence bounds are empirical observed differences, not rigorous error certificates. Other scan points are labeled production-only rather than independently converged.

Largest production/tight/Radau relative differences over these four anchors:

| Model | Rmin | Rmax | max abs(U) | collapse time | paired duration |
|---|---:|---:|---:|---:|---:|
| RP | 1.44512e-05 | 5.21981e-06 | 3.03898e-05 | 1.92445e-07 | 3.43379e-06 |
| KM | 2.72172e-10 | 7.83842e-12 | 7.30302e-10 | 7.58988e-12 | 4.79493e-11 |

The RP 1.10 bar production run has greater sensitivity than the other anchors; this is retained explicitly, with a numerical agreement criterion of 1e-4 relative across the checked metrics. KM uses 2e-5; its actual differences above are much smaller. Root-time tolerance is not physical timing accuracy. The strongest complete RP scan case receives an additional independent audit in strongest_RP_validation.csv.

## Mathematical incompressible limit

The exact explicit expression reduces to the imported RP right-hand side as c tends to infinity. The moderate five-cycle trajectory converges rapidly. The extreme RP core requires much larger c because tiny phase errors amplify later rebounds. At the last multiplier 1e+15, the event-aware maximum radius-trajectory error/R0 is 0.000394336, and speed-trajectory error/reference peak is 0.0616896. This is an observed numerical-resolution limitation of the sharp RP trajectory, not a physical sound speed. Event extrema approach RP much earlier. All factors, nonmonotonic intermediate errors and time shifts are retained; no convergence curve is silently trimmed.

## Extreme KM and physical status

Strongest saved KM: 1.3 bar, Rmin=0.550808707 um, peak inward speed=1072.9808 m/s, Mach=0.724008635. D/R remains positive; the minimum is 0.738452137. The pressure/bulk-modulus proxy is 1.27581699; this is not a computed density change.

Steepest measured KM log-speed interval is [1.115, 1.120] bar on the selected grid. This does not establish a bifurcation or a unique threshold. Global extrema may change cycle; the cycle_comparison file is required for event delays.

Numerical convergence does not establish physical validity. Shape instability, gas-pressure uniformity, thermal/mass transfer, constant kappa, vapor pressure, excluded volume, continuum and EOS validity remain unresolved. Tensile liquid pressure does not prove impossibility; cavitation nucleation statistics are absent. The strongest KM extrema require a validated Gilmore/EOS comparison and gas-closure checks. The follow-up cases are 1.0, 1.15 and 1.3 bar; no unverified Gilmore model or CFD is included. No temperature or reaction is evaluated.

## Tests in final snapshot

```text
.....................                                                    [100%]
21 passed in 33.90s
```

Every report page and all 13 PDF/PNG figure pairs are visually checked before publication. Publication and clean-snapshot verification are recorded separately; file availability and remote SHA are verified after pushing.
