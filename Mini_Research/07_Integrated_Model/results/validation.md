# Numerical validation journal

This journal distinguishes equation tests, empirical convergence, numerical
completion and physical interpretation. It does not certify experimental validity.

## Compact publication scope

The user authorized a selected public package and explicitly rejected publishing
the roughly 1 GB full archive. The full scientific snapshot and original files
remain local. The public history starts at accepted MR6, so the large unpublished
commits do not travel as ancestors. Source, report, aggregate evidence and all run
metadata are retained; selected exact references support continuation and figure
regeneration. The source-array hashes and exact selections are recorded.

The historical 115-run full constitutive replay and 31-test clean full-archive
check are preserved as records. Current compact checks re-evaluate physical powers
from the four shipped accepted-state baselines, optical histories where supplied,
integrity of every shipped array, all metadata and aggregate scientific assertions.
They do not claim to rerun omitted fields. The unchanged scientific PDF and every
figure/table are compared with the original package after compact regeneration.

## Development findings retained

- Initial automatic finite-difference Jacobian adaptation in SciPy overflowed
  inactive zero-state columns. An explicitly evaluated complex-step Jacobian
  removed this numerical issue; no physical state was clipped.
- An overly uniform 1e-15 absolute tolerance for the variational state made the
  initial pulse step smaller than floating-point spacing at nonzero absolute
  time. The incident-photon ledger has a much larger scale than the tiny gas
  response. Explicit component tolerances separate those scales. This is a
  numerical error-control choice, not a modification of heat flow or source energy.
- The first material-volume equality test failed at a 2e-11 relative threshold;
  direct subtraction of nearly equal physical r³ values at large R produced a
  4.8e-11 cancellation error. The independent check now allows 1e-10. Capacities
  themselves are computed from fixed xi differences and never reconstructed by
  that cancellation-prone expression.
- Piecewise constant particle/water sampling had derivative jumps as the moving
  Au surface crossed cell edges. The final source integrates its angular surface
  density against piecewise linear radial reconstruction and uses adjoint
  deposition. It retains normalized positive weights and the same thermal balance.
- Some interrupted development shells left child workers alive. Their exact
  MR7 working directories and command lines were verified; obsolete worker PIDs
  were stopped before final results were admitted. Final tables and raw metadata
  are checked for scenario consistency; corrupted development log fragments are
  not used as numerical evidence.
- Direct subtraction of finite laser and baseline Tmax can reverse sign when
  the true response is far below absolute integration error. Such values remain
  in pairs.csv and convergence.csv. The variational equation and separately
  evaluated nonlinear increment formulation supply cancellation-controlled checks.

## Final checks

The final validation summary, test log, convergence table, independent quadrature,
source hash verification and report-value audit are separate saved artifacts.
BDF/Radau and mesh/domain/time errors are assessed separately. The static MR6
check matches its parameters and window without fitting; local Au temperature
and integrated gas energy have different errors. Material, chemistry, geometry,
optical and compressibility failures are retained in applicability.csv.

Every final PDF page is rendered to JPG for visual review. PDF author/title,
fonts, tables, formula layout and raw-data-driven values are checked. Rebuilding
figures or PDF does not start numerical production. The file manifest includes
all admitted source, raw numerical, scientific and verification artifacts;
environments, caches, temporary rendering and LaTeX auxiliaries are excluded.

A targeted constitutive audit found that one early F=10 raw archive had been
written by an obsolete prototype worker after its replacement run. Its Au-water
power did not match the final normalized reconstruction at the same saved state.
That archive was rejected and rerun with the final kernel; all admitted dynamic
archives are independently checked against current physical powers at saved
states. This is distinct from ordinary integration tolerance error.

For a trajectory censored before any accepted radius minimum, the timing label
is an observed-interval endpoint reference. It is never presented as a completed
first-cycle minimum. The source remains at its archived absolute pulse center.

The same constitutive replay rejected a second obsolete prototype archive,
conv_production_tangent. The current production baseline and its variational
trajectory were freshly recomputed. Final admission requires every archived
power to match the current kernel; matching a summary number alone is insufficient.

At the zero gas-conductance limit, the reconstructed *liquid* wall temperature
is the insulated first-cell temperature, rather than the decoupled gas temperature.
The diagnostic was corrected and separately tested; active-conductance production
trajectories and powers were unaffected.

The first clean export contained MR4/MR5 but omitted their transitive MR3/MR2
dependencies, so import checks failed before any numerical work. A second
partial export passed the fresh numerical replay but its full input-preservation
test correctly rejected an omitted root file. The final check uses the entire
tracked Git tree, with the original 1209-file inventory test unchanged. These
were incomplete export selections, not source or numerical corrections.
