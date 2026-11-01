not_delivered

# Mini-Research 7 scientific review handoff

The user directly authorized compact publication and sending this review, while explicitly rejecting publication of the roughly 1 GB full archive. Compact publication is being verified; delivery is not yet claimed. All large local files and the full scientific snapshot remain available locally. The published history starts at accepted MR6 and excludes the full local commits.

Accepted input and last verified remote base: `493dbdc213f3789b7195d93760bd5c0c867cb98c`.
The local commit and final package inventory identify the reviewable artifact.
Official stage closure and MR8 are outside this work.

## Artifacts and independent reproduction

- Report: `report/07_Integrated_Model.pdf`; source: `report/07_Integrated_Model.tex`.
- Scientific navigation: `README.md`, `07_Integrated_Model.ipynb`.
- Exact package inventory: `FILE_LIST.txt`; SHA256 checks: `SHA256SUMS`.
- Inputs and equations: `results/parameters.json`, `model_decisions.md`,
  `source_audit.md`, `input_inventory.json`, `data_schema.md`.
- Numerical evidence: `results/references/`, `results/metadata/`, `all_runs.csv`, `pairs.csv`,
  `stable_pairs.csv`, `timing_sensitivity.csv`, `convergence.csv`,
  `extra_controls.csv`, `strong_collapse_resolution.csv`,
  `independent_Radau.csv`, `constitutive_replay.csv`, `paired_energy.csv`,
  `physical_diagnostics.csv`, `missing_channel_scales.csv`, `applicability.csv`.
- Validation: `results/tests.txt`, `replay_check.txt`, `validation.md`,
  `input_preservation.json`; final PDF/package checks are recorded separately.
- Saved-data presentation: `python -m src.presentation`; full expensive replay:
  `bash reproduce_all.sh`; bounded fresh replay: `python -m src.replay_check`.
  Python 3.12 and `requirements-lock.txt` match the archived environment.

## Model and source contract

The primary Level 3 state couples R, U, NASA9 fixed-mixture gas temperature,
material-water cell temperatures and Au storage. The material coordinate
ξ=(r³−R³)/3 removes incompressible radial advection, with dV=4πdξ and fixed
positive capacities. Current spherical conductances change with R. A finite
material bath boundary has an independently evaluated signed heat flux.

One full current-state Qlg=4πR²hg(Twall−Tgas) enters gas and liquid with opposite
signs. The MR4 clamped-water loss is replaced. No MR6 frozen baseline, sum or
constant efficiency is added. MR6 contract v2 is used as a reference boundary;
Level 1 is contextual rather than a moving-geometry result.

pg=nRuT/V, dpg=pg(dT/T−3U/R), Cv(T)dT=−pg dV+Qlg feed the full viscous KM
equation, including the external pressure derivative. H2O is already in the
fixed inventory; there is no extra constant vapor pressure. No reacting EOS or
phase-change temperature correction is claimed.

Actual isolated-host MR5 Mie absorption supplies Pabs=σabs I once. Au-to-water
power is removed from Au and deposited through normalized adjoint surface
weights, using a spherical mean water temperature. This omits the finite-Au
local spreading resistance. N is a local count; the beam photon allowance,
particle ring geometry, surface gap and contact/outer-escape events are checked.

## Principal finite nonlinear pairs

R0=6 µm, f=30 kHz, T0=293.15 K; one a=20 nm Au particle, nominal gap100 nm,
532 nm, G=10⁸ W/(m² K), 10 ns top-hat F=10 J/m², 100 µm beam radius. The
position is prescribed rnp=R+a+d, with no validated retaining force. Pulse center
is five microseconds before each no-laser selected radius event. Five cycles are
requested, I=[0,166.667] µs. Each heated pair keeps the same absolute pulse time.

|Pa (bar)|Observed baseline (K)|Observed heated (K)|Global ΔTmax (µK)|t baseline / heated (µs)|Scope|
|---|---:|---:|---:|---|---|
|1|553.0230367272|553.0230381360|1.4088163|147.13186439767 / 147.13186440028|complete five cycles|
|1.15|2798.2812622841|2798.2812644690|2.1848791|16.23774009043 / 16.23774009120|complete five cycles|
|1.3|6000.0000000010|5999.9998704094|undefined|18.54531415805 / 18.54531415805|coefficient/reference censored|

These pairs share a checked C1 no-laser reference and a nonlinear directional
Jacobian-integral formulation (three Gauss nodes, five-node check). Full nonlinear
independent Radau subtraction at 1.0 bar gives ΔTmax=1.40727445 µK, within 0.11%
of the stable value. Absolute baseline error is larger than the tiny increment;
the displayed baseline digits are reproducibility data, not physical accuracy.
Direct production BDF differences remain archived and unresolved at this scale.

At 1.0 bar, Level 2 → Level 3 is 0.857788 → 1.408816 µK; feedback adds
0.551029 µK (about 64%). At 1.15 bar it is 5.656307 → 2.184879 µK; feedback
reduces it by 3.471428 µK (about 61%). There is no general beneficial-feedback rule.

Initial net wall-flux increments peak about39/35 ns after the pulse center;
the first-cycle thermal peaks occur about2.46/5.00 µs later at 1.0/1.15 bar.
The pulse can heat gas before these peaks, but the response is tiny. At 1.0 bar
the selected first-cycle radius event follows the thermal peak by about2.54 µs.
The five-cycle global peak is later, so first-cycle causality and global memory
are reported separately. Net wall-flux increments include feedback and are not
laser-tagged delivery fractions.

## Censoring and formal extrema

The wet 1.3 bar model stops near 18.5453 µs at 6000 K. Its full maximum and global
paired difference are undefined. Common-prefix diagnostics are labeled separately,
including the stable negative endpoint change; an unknown tail is not inferred.
When no accepted radius minimum precedes censoring, the pulse reference is the
observed-interval endpoint, not a completed collapse event.

Dry 1.3 bar baselines reach 6263.579 K (air) and 22808.716 K (ground-state constant-Cv
Ar), with complete mathematical five-cycle trajectories. Their F=100 pulses
use a separate minus 1 µs protocol. These are fixed-composition formal sensitivities,
not temperatures of the wet bubble and not validated laser gains. Wet H2O ends
at 6000 K; N2/O2 coefficients end at 20000 K; Ar's constant-Cv law has no finite
coefficient ceiling. All 1.3 bar cases remain formal.

Retained solver/event statuses across the dynamic archive:

- `completed`: 96 dynamic archives.
- `gas_coefficient_upper_boundary`: 12 dynamic archives.
- `particle_contact`: 3 dynamic archives.
- `Au_moderate_rise_30K`: 2 dynamic archives.
- `Required step size is less than spacing between numbers.`: 1 dynamic archives.
- `reference_censored_end: gas_coefficient_upper_boundary`: 1 dynamic archives.

Static/analytic controls are separate files. Every run's requested interval,
actual end, event type, gas domain, completion and observed/global values remain
in its JSON and the aggregate tables. Failed strict pulse-edge checks remain
visible, with the repaired incident-ledger tolerance documented.

## Conservation, convergence and physical screens

Gas, water and Au balances are checked separately and together:
Δ(Ug+El+EAu)+Wbygas+Qbath=Eabsorbed. Every physical power is independently
computed, signed, and integrated; bath heat is not a balancing residual.
This is the thermal subsystem with mechanical work, not the total acoustic-fluid
energy budget. The stable pairs have maximum incremental total residual
2.583e-07 of absorbed energy. The historical full local validation admitted 115 dynamic archives with saved-field
constitutive power replay. The compact package independently replays four shipped full-state baseline references; the other original fields remain local. All 1209 previously tracked input files remain unchanged.

At 1.0 bar the 48→96-cell derivative changes 0.30%; at 1.15 bar the early derivative
changes 0.6% and independent Radau changes 0.21%. The extremely small late pulse
centered on collapse at 1.15 bar varies by nearly an order of magnitude under
mesh refinement and is explicitly **unresolved**. The pure outer-boundary test
preserves all original inner faces; the older outer50 test also changes inner
resolution and is not a pure domain-error estimate.

Clamped MR4 regression recovers both completed peaks to better than 10⁻⁶ relative;
the censored 1.3 case is checked by status and end time, not just its imposed 6000 K
boundary. The analytic transient Robin point control has at most 0.326% relative
error at sampled positive times on 320 cells. Matching MR6 static parameters gives
selected integrated gas-delivery differences of 0.04–1.68%, yet Au peaks are only
21–30% of MR6. Matching nearest bath clearance retains bath-shape differences.
These controls quantify a large local-Au model error rather than hide it.

High gas density/pressure, first-order KM, incompressible thermal advection,
constant gas conductivity, fixed composition, omitted evaporation/condensation,
disassociation/ionization, isolated-host optics and finite-Au hot spots remain
limitations. The wet 1.3 path approaches 2 GPa and 1157 kg/m³ in its ideal EOS.
Separate formation/phase/ionization reference energy scales do not evolve reaction
extents or correct temperatures. Material thresholds are diagnostic references,
not universal nanoparticle damage or nucleation limits. No matched experiment
validates the present combined acoustic/particle/optical configuration.

## Questions for independent scientific review

1. Are the moving-coordinate derivation, boundary resistances, source weights and
   signed balances correct, without MR4/MR6 double counting?
2. Does the checked reference / nonlinear increment / independent Radau evidence
   support the reported conditional microkelvin changes and the unresolved late tail?
3. Are radius-event timing, thermal-peak causality, absolute pulse time and
   censored common-prefix comparisons separated correctly?
4. Is the static finite-Au discrepancy sufficiently restrictive on moving-model
   interpretation, especially before using any local material-validity screen?
5. Are the formal dry extremes, gas domains, EOS/phase/chemistry and mechanical
   limitations stated with adequate force? No physical temperature ceiling,
   achievable maximum, plasma state or fusion yield is asserted.

Independent review and RoadMap closure remain outstanding. Do not begin MR8
from this handoff.

## Clean archived snapshot verification

The complete tracked snapshot `a809c8240328460d69c5415b042ca705a0846f40` passes all 31 tests, all 367 manifest checks and the bounded fresh numerical replay. See `results/clean_archive_check.json` and `results/clean_archive_tests.txt`. Subsequent local changes are verification records, whitespace cleanup and reproduction documentation; the scientific arrays and PDF are unchanged.

## Compact public data and audit coverage

The complete source, report, figures, parameters, aggregate results, convergence/applicability tables and all 115 dynamic metadata records are shipped. Four accepted-state references cover the three acoustic amplitudes and the separately used censored high-amplitude pair. Selected exact scalar traces cover the main finite nonlinear responses, geometry controls and formal dry extremes. The publication inventory records every retained array and its source hash; five water profiles are exact indexed snapshots. No interpolation, rounding, clipping or refitting was introduced.

Other high-dimensional sweep outputs and solver matrices remain local. The full local snapshot c467678d2f48db51ceafe09ea4cb5040b6ecc293 is not published or included as a public ancestor. Historical full-archive checks are preserved separately from the current compact tests and reference checks. These distinctions limit what an independent reviewer can replay from the public subset without regenerating or requesting additional local arrays.

Current compact verification: 34 tests passed; four shipped accepted-state references replay their physical powers; eight shipped optical histories pass pulse-energy checks; the fresh one-cycle replay agrees with the preserved baseline. All 10 JPGs, report tables, PDF and TeX match the full original scientific package byte-for-byte. See `results/compact_verification.json`, `results/reference_checks.csv` and `results/tests.txt`.
