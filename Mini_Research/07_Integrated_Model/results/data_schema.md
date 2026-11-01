# SI data contract

## Public subset versus full local archive

`publication_scope.json` and `public_data_inventory.json` define the compact
public profile. It ships `references/<run>.npz` and all 115 `metadata/<run>.json`
files; `raw/<run>` below describes the complete local/regenerated archive.
`src.data_access` selects the public profile by default, or the full archive only
when `MR7_DATA_MODE=full`. The four public baseline references retain every
accepted state. Principal nonlinear pairs retain exact full scalar traces and
their six incremental ledgers; the figure reference retains five exact liquid
profiles. Other selected geometry/static/dry traces contain only documented
arrays. A missing array is not zero and must not be reconstructed silently.
All retained arrays have individual hashes; selected snapshots record original
array hashes and exact index selections. No numerical data were interpolated,
rounded or changed to create the public subset.

`full_archive_tests.txt`, `constitutive_replay.csv` and `clean_archive_check.json`
record the preceding full local validation. `tests.txt` and `reference_checks.csv`
are current compact-package checks. Their coverage is deliberately distinct.

All physical inputs and archived arrays use SI. Bar display columns are explicit:
1 bar=100000 Pa. Solver time internally uses microseconds; archived t_s is seconds.
NaN (CSV) or null (JSON) means unavailable. Infinity denotes an unbounded analytic
coefficient domain, not physically unlimited temperature.

raw/<run>.npz has t_s, R_m, U_m_s, T_g_K, p_g_Pa, T_l_rise_K[cell,time],
T_Au_rise_K[time], reconstructed T_wall_K, r_np_m and surface gap_m. Water radii at
any time follow r=(R^3+3 xi_centers_m3)^(1/3). The field is spherically averaged;
no local Au hotspot is resolved. Stored material capacities are constant, J/K.
accepted_t_s and accepted_y retain solver nodes and scaled states: R/R0, U/c,
Tg/T0, water rises/T0, Au rise/T0, then W,Qlg,QAu,Qbath,Eabs,Eincident divided
by n Ru T0. Tangent archives append the unit-fluence variational states. These
are accepted nodes; they do not reconstruct the original BDF/Radau polynomials
exactly. Expensive replay reintegrates. Level 2 uses a separately checked C1
Hermite reconstruction of Level 3 no-laser R with stored U as its derivative.

Flows Q_l_to_g_W, Q_Au_to_l_W, Q_bath_W have signed physical directions. No
channel is defined from a residual. Integrated counterparts and U_gas_change_J,
E_water_J and E_Au_J are independent ODE/temperature-derived quantities.
The residuals are gas: DU+W-Qlg; water: El-QAu+Qlg+Qbath;
Au: EAu-Eabs+QAu; total: DU+El+EAu+W+Qbath-Eabs.
This is thermal-subsystem bookkeeping with prescribed mechanical work. It excludes
fluid kinetic energy, external acoustic drive, viscous heat and confinement work.
Positive/negative integral parts in summaries are independent 8-point Gauss rules
on accepted intervals; signs are retained. Signed gas heating can be negative.

summary.csv and per-run JSON include actual end, terminal status, requested end,
observed versus global extrema, separate temperature and radius event times,
component coefficient domain and physical screens. Main initial placement is at
t=0, not at pulse time. Event tc is the first local R minimum with R<R0/2;
otherwise it is cycle-1 deepest minimum, consistent with the MR6 fallback.
First strong minimum, first-cycle deepest, global R minimum and thermal maximum
are separate. Global extrema are defined only for a complete five-period run.
If no accepted local radius minimum precedes censoring, the explicitly named
`observed_interval_endpoint_reference` is a timing reference only. It is not a
completed collapse or first-cycle event.

pairs.csv compares the same Level 3 model or the same Level 2 model within each
pair. Pulse center is absolute and fixed by the no-laser tc; no recentering occurs.
If either run is censored, global_delta_Tmax_K is NaN and only an explicitly
labeled common-prefix value on [0,min(end0,endLaser)] is available. Observed maxima
on different intervals are individually reported and never subtracted as global.

Tangent dT_dF_K_m2_J and dTmax_dF_K_m2_J are derivatives at F=0, not finite
nonlinear temperature differences. The peak derivative is evaluated at the
selected no-laser thermal maximum; near ties require care. Censored endpoint
derivatives are observed diagnostics, not global sensitivities. dledgers_dF_J
retains all six differentiated energy channels. The total first-order residual
uses Cv(Tbase)*dT, not the nonlinear sensible energy of Tbase+dT. Finite-difference
checks are separately archived. Timing offsets are pulse center minus no-laser tc.

static_MR6_comparison.csv matches R,a,d,h,G,k,F,duration and the finite observation
window. Its 1D/2D differences include finite Au exclusion, angular thermal
resistance and finite external-domain shape. They are not universal error bounds.
The radial source is the exact angular surface average tested against linear
water reconstruction; its adjoint deposits energy once. Static raw archives use
water/gold temperature rises, actual fluxes and independent stored energies.

convergence.csv separates mesh, outer-domain, solver-step/tolerance and independent
method effects. It includes extrema, time, Rmin, tangent effect and material-field
norms. Local_sensitivity.csv is a finite one-at-a-time design, not a global optimum.
All scientific figures/tables are regenerated from saved files, without solving.

## Publication storage version 2

All scalar histories and their extrema remain at every original saved time.
Liquid arrays T_l_rise_K, dTl_dF and paired_delta_Tl_K use field_t_s, with
field_indices mapping exactly into t_s. These are original float64 snapshots,
not interpolated or rounded values. Selection retains two uniform diagnostic
grids, every liquid-cell minimum/maximum and pulse/thermal-event neighborhoods.
No-laser reference files retain every accepted_y state for replay. Other files
omit this redundant solver-state matrix. storage_compaction.json records the
original and retained field counts and verifies unchanged scalar arrays. A full
numerical regeneration produces larger intermediate archives before compaction.

For native prescribed-mechanics Level 2, the first two accepted_y variables are
inactive placeholders; physical R_m and U_m_s contain the imposed C1 path.
Stable paired Level 2 is explicitly tagged with its coupled no-laser reference.
The older dTmax_dF column is interpreted through response_scope: after censoring
it denotes only the formal fixed-reference endpoint response. The separate
 global_dTmax_dF_K_m2_J is null for those cases. No unavailable tail is inferred.

Stable nonlinear pairs additionally store paired_delta_T_K, paired_delta_R_m,
selected exact liquid increments, all six paired_delta_ledgers_J histories and
four independently formed paired subsystem/total residual histories in J.
The latter are formed before field compaction and remain complete scalar traces.
After compaction, scientific_audit rechecks their recorded maxima without
inventing missing liquid-field times. feedback_decomposition.csv samples the
changed work and heat at the no-laser thermal reference time; it is distinct
from the difference of independently maximized temperatures.

Physical diagnostics report pressure in Pa, density in kg/m³, R/cgas in seconds
and sampled half-height temperature widths in seconds. Standard-state missing
channel scales are J and dimensionless ratios to observed frozen sensible
energy. They overlap and are not additive reaction or phase-change predictions.
Boolean material/optical/Mach screens are diagnostic conditions, not evidence
that omitted local fields, retaining forces or an experimental state are valid.

The final report table presents two independent documentation columns:
`Numerical completion` and `Physical applicability`. This does not change any
CSV field or classification. `completed=True` records integration completion
only. In particular, the 1.15 bar completed pair remains a formal fixed-composition
scenario, the 1.00 bar completed pair remains conditional on unvalidated prescribed
confinement, and the wet 1.30 bar global five-cycle difference remains unavailable.
