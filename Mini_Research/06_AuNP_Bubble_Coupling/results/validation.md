# Final-correction verification: frozen-geometry MR6

The agreed final scope is stationary AuNP-liquid-bubble thermal coupling.
0.3258% is a conditional static-geometry result for the specified 2 microsecond
window, not a physically verified efficiency during acoustic collapse. Moving
boundaries, advection and self-consistent mechanics are deferred to MR7.
Publication is not independent acceptance or official stage closure.

## Original data and scope preservation

The correction starts from published MR6
`faa2828f554517a0e2e43d9ae42cff37d682fbd0`; its initial 117 checksums passed.
All 1083 accepted upstream fingerprints pass again. MR1-MR5, Plasma-plan.tex
and unrelated user files remain unchanged. No new parameter optimization,
large sweep or mechanical/gas trajectory calculation was performed.

`correction_input_manifest.json` verifies 55 original numerical files byte-for-byte
(21 CSVs, 29 non-source NPZs and five scientific JSONs). All original time, laser,
baseline and sum arrays in 24 source NPZs are numerically identical after explicit
reference-state renaming. Reference gas temperature and radius are now included.
The historical full-run record remains in `clean_reproduction_original.json`.

## Corrections checked

- Handoff late-window fraction corrected from 0.01817% to 0.00817%; the original
  CSV value is 0.008166454345095196%. No underlying result was changed.
- 59 rounded scientific claims and bounds pass their source-data audit, including
  the six newly reported feedback estimates/bounds. See `handoff_numeric_audit.csv`.
- The failed mesh test no longer regenerates a mesh to assert bitwise equality.
  Archived nodes/triangles/edge markers are explicit inputs. The reported 1405
  versus 1404 node difference is therefore immaterial to stored-field validation.
  A test disables the mesher entirely and still loads the event_1.00 field.
- Five archived fields pass independent degree-2 triangle quadrature of water
  energy, wall-flux integration and balance/temperature checks. Signed far-tail
  undershoots remain unmodified.
- Three bounded saved-mesh replays compare temperature, power and energy at a
  2e-5 max-norm tolerance, 100 times the 2e-7 time-solver tolerance. This is a
  tolerance for replay on known meshes, not an error certificate for new meshes.
- The v2 API exports reference baseline/sum and laser increment separately, with
  no q_wall_total_W. An ODE already evaluating current-state heat loss adds only
  laser_increment_W(). Tests check changed gas temperature, absence of ambiguous
  keys, source-state consistency, legacy loading and no double-counted baseline.
- Mechanism A is explicitly a transient early-heating diagnostic, not a new
  equilibrium R0/vapor-pressure/composition model.
- Figure 04's vertical label is shortened to fit; Figure 07's A legend and the
  saved-results notebook use early-heating terminology. The PDF is rebuilt.

## New bounded checks versus retained original convergence

Fresh replay worst observable relative difference: 1.4755496e-08 (admitted 2e-5).

Independent stored-water energy residual / absorption: 2.412255361e-15.

The original mesh/time/domain/isolated studies are retained rather than repeated:

| Original selected-case metric | Relative discrepancy |
|---|---:|
| production_gas_energy_spatial_relative_max | 0.001956097602 |
| production_gold_spatial_relative_max | 0.01205506181 |
| production_water_spatial_relative_max | 0.01191987823 |
| isolated_gold_peak_finest_relative_max | 0.0008923868612 |
| independent_Radau_gold_relative | 1.411301934e-09 |
| independent_Radau_gas_relative | -1.905364755e-12 |
| energy_balance_relative_max | 2.657475174e-14 |

These are empirical selected-case differences, not universal error bounds.
The static MR4-boundary matching discrepancy remains 4.2146%, with no bound
claimed for unsteady matching. All 45 existing event snapshots remain formal.

## Omitted gas-to-liquid feedback estimate

The estimate uses only the saved gas increments, an unobstructed stationary sphere
in infinite water and an analytic surface step response. It does not change the
original gas trajectories, liquid finite-particle calculations or mechanics.
The gas ODE already includes its own -4*pi*k_g*R(t)*delta_Tg loss derivative;
no second direct loss is added. Only the possible extra positive liquid-return
term is estimated. For B, return/source-energy ratios are approximately 0.2323%,
0.09284%,0.000225% at 1.0,1.15,1.3 bar. The 1.3 bar result covers its censored
observed interval only. Maximum-principle bounds are much looser: 74.3%,12.1%,0.127%.
Thus the small estimates do not prove universal negligibility or a physically
validated collapse efficiency. Finite-Au obstruction, flow and changing radius
are outside this analytic control. `gas_feedback_estimate.csv` retains nine
A/B/ideal-C controls, already-counted gas loss, bounds and two-grid differences.

## Clean environment and presentation

A new Python 3.12.14 environment was installed from the published dependency lock;
pip check passed. A separate repository archive without .git was overlaid with
correction sources while retaining the original v1 source files. The fresh run
therefore re-exported those sources through the new schema, ran all bounded
checks/replays and all 31 tests, rebuilt all eight JPGs and compiled the PDF.
No original sweep was repeated. The five field checks and the mesher-bypass test
verify the corrected portability contract; a second operating system was not run
in this session and is not claimed.

All 53 preserved/re-exported NPZ files (586 arrays, including metadata/state
arrays), source conversion results, 59 claim checks and bounded numerical checks
agree with the clean run. Replay comparison uses 2e-5, analytic feedback reproduction
uses 1e-6 relative tolerance; tiny numerical residuals are checked by absolute
physical bounds rather than required to have identical roundoff digits. All eight
JPGs happen to be byte-identical on this host; PDF text also matches exactly.
This is not a demand for cross-platform remesher or raster byte identity.
`clean_reproduction.json` and `clean_environment_tests.txt` contain actual results.

All eight final PDF pages were rendered and visually inspected. Figure 04's full
vertical label, equations, signs, units, legends, Appendix feedback formula and
references are readable without clipping or overlap. The final LaTeX log has no
warning, overfull box or undefined reference. The scientific PDF contains no Git
SHA, Codex/Web instructions, internal correspondence or acceptance status.

## Publication gate

FILE_LIST and SHA256SUMS cover the final package, excluding environments, caches
and intermediate TeX files. The staged Git archive is verified for manifest
completeness, checksums and all tests before commit. Only MR6 paths may be staged.
Remote main and immutable artifacts are checked after the non-force push. The
publication dispatch supplies the actual new SHA, links, counts and delivery
acknowledgment to the current Mini-Research 6 (Web) chat. MR6 is not declared
closed, RoadMap is untouched, and MR7 is not begun.
