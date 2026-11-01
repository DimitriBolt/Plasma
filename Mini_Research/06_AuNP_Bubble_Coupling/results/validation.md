# Numerical and artifact verification

The checks establish reproducibility of the stated conditional models; they do
not establish experimental validity or independent scientific acceptance.

- Starting local and remote main: `0204ce7a8429406b7ec235b02cec8581274fc8c6`.
- All 1083 tracked starting files were compared byte-for-byte to accepted Git
  blobs and SHA-256 inventoried. No prior stage, plan or unrelated user file changed.
- 24 tests pass in the production environment. Tests include saved-field storage
  and maximum-principle checks across five serialized physical geometries, rather
  than only checking the internally integrated scalar energy ledger.
- All CSV and all NPZ numeric arrays are readable, finite where defined, SI
  documented and retained with validity/censoring flags. Formal extremes are not
  clipped. Tiny signed far-tail FEM undershoots are retained.

## Selected quantitative checks

| Check | Result |
|---|---:|
| spatial_cases | 15 |
| time_method_cases | 3 |
| domain_cases | 6 |
| production_gas_energy_spatial_relative_max | 0.001956097602 |
| production_gold_spatial_relative_max | 0.01205506181 |
| production_water_spatial_relative_max | 0.01191987823 |
| isolated_gold_peak_finest_relative_max | 0.0008923868612 |
| independent_Radau_gold_relative | 1.411301934e-09 |
| independent_Radau_gas_relative | -1.905364755e-12 |
| energy_balance_relative_max | 2.657475174e-14 |
| gaussian_normalization_relative_max | 2.220446049e-16 |
| mechanics_input_radius_difference_over_R0_max | 6.251416331e-06 |
| baseline_static_matching_error | 0.04214559387 |

Spatial differences are empirical selected-case discrepancies, not rigorous
universal bounds. Independent sampled wall-flux integration differs from the
ODE ledger by at most 7.405e-6 relative in its selected control. The gas increment
2 ns/0.5 ns step-cap comparison differs by at most 8.538e-7 relative. The isolated
comparison uses exactly matched Au/water coefficients and does not validate the
new bubble geometry. The analytic finite-h point-sphere control is independently
derived; early vanishing far tails are explicitly not relatively resolved.

Domain/insulated controls distinguish gas delivery convergence from artificial
bath escape. Expanding 1 to 2 micrometre padding changes reference gas delivery
about 0.2%, while water-versus-bath storage changes materially. The 0.5 micrometre
domain is too close for the same result. Full sphere/planar comparisons and all
high-h limits remain available.

## Physical and model limits

All 45 event snapshots fail the declared frozen-motion screens. Their event-aligned
outputs remain formal diagnostics, including the conditional mechanical baselines.
The 1.3 bar flow has additional inherited compressibility limitations. Passive
particles are not asserted to remain at a prescribed wall gap. Near-interface
optics, capture forces, hot transport, mass transfer, phase change and chemistry
remain unresolved. The 4.2146% static MR4-boundary matching difference exceeds the
tiny one-particle gas-temperature effects and is not an unsteady error bound.

The gas calculation evolves a common nonlinear baseline plus its first-order
laser perturbation to resolve tiny increments without subtractive noise. Its own
perturbed work/heat/source identity is checked. Large pressure mismatch with the
original mechanical gas is retained. The actual first NASA boundary censors each
strong history; no 6000 K physical ceiling, missing tail or global censored Tmax
comparison is claimed. Late sources have no effect on preceding peaks.

## Report and package

The final English scientific PDF has eight pages: seven main pages and one
appendix including references. All eight pages were rendered to JPG and visually
inspected. An initially orphaned paragraph and an overlong figure axis label were
corrected; the latest pages are readable without clipping or overlap. The
serialized water field was independently checked against water storage and the
maximum-temperature bound. Equations, legends, units and references are readable.
The final LaTeX log has no overfull box, undefined reference or warning.

Eight figures are JPG only, 300 dpi, rebuilt from saved data. The numbered
notebook has Markdown result text and embedded saved figures for GitHub reading.
The scientific report contains no Git SHA, Codex/Web task instructions, internal
chat references or stage-acceptance status. The source interface documents thermal
memory, finite intervals, signed baseline/increment and frozen-radius admission.

## Clean full recalculation

A separate archive of the accepted base plus the MR6 sources was extracted without
.git. A newly created Python 3.12.14 environment was installed from the full lock;
its dependency check passed. All experiments, finite designs, event diagnostics,
mesh/time/domain/analytic/isolated controls, gas-step refinement, figures and PDF
were recomputed. Four additional full-field benchmarks were independently rebuilt
in that same fresh environment. The final tests, including their field/storage
checks, passed: 24 tests.

All 21 CSV tables, all 538 arrays in 53 NPZ archives, all eight JPG files, five
scientific JSON documents and extracted PDF text match exactly. The clean PDF
rebuilt successfully. This is same-machine fresh-environment reproduction, not
cross-platform byte identity or independent physical validation. `clean_reproduction.json`
and `clean_environment_tests.txt` retain exact outcomes.

## Final publication procedure

The explicit FILE_LIST and SHA256SUMS include the full package and exclude caches,
virtual environments and TeX intermediate files. A final staged Git archive is
checked for file-list equality, all checksums and all 24 tests before publication.
The accepted upstream fingerprints are rechecked, and staged paths must all lie
inside MR6. The remote main is read again and the push does not use force.

The publication dispatch records actual remote SHA, immutable artifact links,
remote-file comparisons and delivery acknowledgment after they occur. The current
Web destination was identified/read as Mini-Research 6 (Web); the actual proven
mechanism is send_message_to_thread followed by destination readback. No previous
stage's destination identity is reused. Publication does not close the stage or
begin MR7.
