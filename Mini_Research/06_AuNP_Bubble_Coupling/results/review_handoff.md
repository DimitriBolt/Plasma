# Mini-Research 6: independent scientific review

Please audit the mathematics, numerics, physical interpretation and reusable
source interface of **Mini-Research 6: Frozen-Geometry AuNP-Bubble Thermal Coupling**. Publication
does not constitute scientific acceptance or RoadMap closure. No MR7 calculation,
global coupled optimization or reaction-rate calculation is included. The final agreed
scope is stationary AuNP-to-liquid-to-bubble heat transfer. Moving boundaries,
advection and self-consistent mechanical interaction are deferred to MR7.

This correction follows independent review of the original MR6 snapshot
`faa2828f554517a0e2e43d9ae42cff37d682fbd0`.

The publication dispatch supplies the full remotely verified commit SHA and
immutable commit/PDF/package links. This file belongs to that containing snapshot;
`git rev-parse HEAD` identifies it in a clone. The accepted starting main was
`0204ce7a8429406b7ec235b02cec8581274fc8c6`. All 1083 starting tracked files are
byte-compared to their Git objects and SHA-256 inventoried. MR1-MR5, the plan
and unrelated user files are preserved. Only `Mini_Research/06_AuNP_Bubble_Coupling/`
is published.

The English report is `report/06_AuNP_Bubble_Coupling.pdf`, **8 pages: seven main
scientific pages and one parameter/numerics/references appendix**. Author:
Dimitri Bolt. Academic supervisor: Prof. Dr. Ildar R. Gabitov. The package includes
the shared SI module, numbered GitHub-readable saved-results notebook, eight JPG
figures, all raw CSV/NPZ results, tests, dependency lock, manifests and validation
records. Administrative evidence remains outside the scientific PDF.

## 1. Model to audit

The main branch is a finite-particle, axisymmetric linear heat equation on a
**frozen** geometry. Both the entire Au sphere and entire spherical bubble surface
are explicitly meshed in a meridian, with the exact axisymmetric 2*pi*rho weight.
Sphere boundaries are polygonal chords, with refinement checks. The liquid domain
is a finite cylinder with a fixed-temperature outer bath; insulation and domain
expansion are separate controls. Au is one equilibrated temperature node with
capacity 4*pi*a³*rho_Au*cp_Au/3 and finite Au-water interface G.

The water perturbation obeys rho*cp*theta_t=k*Laplacian(theta). Au input is
`P_abs=I_local*sigma_abs`, with the actual MR5 Mie API and n_medium=1.33. The
Au interface exchange is removed from gold and added to water through the same
boundary links. With n_g from gas into liquid, the positive laser increment is
`delta_q = integral(h*theta_water dA)`. Gas perturbation is held zero in this
laser-tagged calculation. All perturbations initially vanish.

The signed gas energy model uses **MR4 baseline plus only this increment**:

`C_V(T)*Tdot = -p(T)*Vdot - 4*pi*k_g*R*(T-T0) + delta_q`.

It does not add another full boundary flux. The absolute resolved hot-gas
liquid-side baseline during collapse is not solved. Nominal h=k_g/R is an inherited
effective gas-conduction assumption, not measured molecular interface conductance.
A uniform static exterior-water-resistance control changes the baseline conductance
by 1/(1+k_g/k_l)=0.9578544061: a **4.2146% matching discrepancy**, not a bound on
the unsteady mismatch. Please assess whether the controlled small-increment
interpretation is sufficiently explicit.

P1 FEM uses exact constant-coefficient diffusion integrals, positive row-sum mass,
and conservative boundary quadrature. Independent gas/bath flux states check
`E_abs = E_Au + E_water + Q_gas + Q_bath`; bath heat is not defined as a residual.
The stored physical temperature field also independently reproduces water storage.
Full Au/water/power histories and spatial fields for additional near-contact,
far-particle and 1.0/1.15 bar frozen event geometries are in `results/fields/`.
All signed far-tail undershoots are retained rather than clipped. These decrease
with refinement; exponentially small far-source early tails have no reliable
relative accuracy or physical-energy interpretation.

## 2. Main conditional numerical results

Reference: R=6 µm, a=20 nm, **surface gap d=100 nm**, local F=1 J/m², 532 nm,
10 ns top hat, G=100 MW/(m² K), h=4400 W/(m² K), 1 µm outer padding.

- Sigma_abs = 3.4337324474063743e-15 m²; E_abs = 3.4337324474063743 fJ.
- Gold peak rise about 2.3066 K; facing Au-surface liquid peak about 1.6322 K.
  These are rises above 293.15 K, not absolute gas temperatures.
- Gas power peaks near 53.94 ns from pulse start, or 48.94 ns after pulse center,
  at about 14.0702 pW. Gold peaks near pulse end, not at the gas-flux peak.
- By 2 µs: Q_gas about 11.1865 aJ, **0.3258% of local absorption**.
- At the power peak, max local flux about 161.64 W/m² versus whole-surface mean
  0.03110 W/m²: ratio about 5197. About 0.0481% of the sphere carries 90% of
  positive flux. Raw signed fields and nodal area weights permit independent checks.
- The largest reference late-window fraction on the declared timing lattice is
  about 0.00817% for [tc-10,tc+10] ns at laser-center offset -50 ns. This is
  a fixed-geometry timing result, not a global physical optimum.

The finite design has 216 geometry/duration/wavelength cases, continuously capped
by a 30 K Au-rise allowance, 1000 J/m² and 1e10 W/m². At d/a>=5 and the declared
equilibration screen, the largest gas energy by 2 µs is about **3.9517 fJ** from
about **3.48462 pJ** absorption, a=100 nm,d=500 nm,100 ns. The largest absorbed
energy is about 3.48476 pJ at another gap; this tiny gap ranking is below spatial
uncertainty. Maximum fraction instead uses a=20 nm,d=100 nm,1 ns. Thermal-cap
wavelength plateaus tie to roundoff; the largest cross section is the stated
tie-break, not proof of a unique optimum. `design_extrema.json` preserves each
objective and `finite_design.csv` preserves every candidate.

Including d/a=1 gives a formal all-grid gas-energy maximum about 9.43494 fJ,
subject to unvalidated interfacial optics. The near-contact 0.1 ps / high-fluence
scan retains formal sensible-equivalent Au rises in the tens of thousands K.
Those are not actual gold/water states. Time equilibration, continuum,
phase stability, optics and capture failures are enumerated in the data.

N=0,1,2,4,8 local particles are a separately screened, nonshadowing equatorial ring,
not a macroscopic absorption layer. Source energy is bounded by a stated local
beam allowance, separation exceeds optical and finite-time thermal screens, and
a diffusion-tail overlap proxy is retained. Conditional linear patch superposition
is not extrapolated to a dense random dispersion or repeated photon absorption.

## 3. Motion, timing and why event predictions remain formal

For passive incompressible flow, r_np³-R³ is constant. Differentiation gives
`gap_dot=((R/r_np)**2-1)*Rdot`; the exterior gap **increases during collapse**.
Fixed-center and prescribed capture laws are distinct controls. Capture has no
modeled retaining force. Nonpositive gap stops the liquid-only snapshot branch
with a contact/intersection status.

Selected complete mechanical inputs:

| Amplitude | Radius event | tc | Rc | Inherited use |
|---|---|---|---|---|
| 1.0 bar | cycle-1 deepest, no R<R0/2 event | 20.092624 µs | 4.39309 µm | conditional baseline |
| 1.15 bar | first strong R<R0/2 | 15.664319 µs | 1.22816 µm | conditional sensitivity |
| 1.3 bar | first strong R<R0/2 | 18.341064 µs | 0.550809 µm | formal extreme |

R<2Rc widths are about 21.53 ns and 2.278 ns at 1.15/1.3 bar. The separate MR4
coupled thermal half-height width near 17.85 ns is not called the width of the
same mechanical trajectory. New Hermite interpolation is distinguished from the
original dense solution, with decimation and accepted-node sensitivity checks.

Snapshots freeze R and gap at pulse center, then explicitly measure Pe, radius
change over diffusion time and response duration, and gap change. Event-aware
nodes capture the narrow velocity peaks. The late strong-event snapshots fail
the Pe<0.1 / radius-change<1% / gap-change<10% working screens, even for imposed
capture. They are retained **formal frozen-motion-omitted diagnostics**. A zero
instantaneous velocity at tc does not justify a frozen response. The 1.3 bar
incompressible particle law is additionally a formal extrapolation of liquid flow.

## 4. A/B/C and the limited gas calculation

The gas uses the read-only MR4 saturated air/H2O inventory, ideal EOS and NASA9
heat capacities. Vapor is counted once. Windows initialize at tc-3 µs from saved
MR4 temperature and request tc+2 µs on the unchanged MR2 mechanics. They are
observed-window diagnostics, never new five-cycle global temperatures.

A common nonlinear baseline is integrated once. To avoid subtractive noise from
tiny one-particle effects, the first-order sensitivity obeys
`delta_Tdot=f_T(T_base)*delta_T+q_laser/C(T_base)` with analytic NASA C'. Separate
incremental work/loss/source ledgers check its first law. A is a **transient early-heating diagnostic**, with an early diffusion pulse
(-2 µs); it does not construct a new equilibrium radius, vapor pressure or gas
composition. B uses near-collapse diffusion (-20 ns) and an after-collapse control (+200 ns).
C is a finite-normalized Gaussian centered -5 ns, sigma=2 ns on [-20,+10] ns;
its energy is either computed B transfer over 5 µs or explicitly ideal 100% local
absorption. No energy or Gaussian delay is fitted to a desired temperature.

At 1.15 bar, increments at the baseline thermal peak are approximately:
A 1.887e-5 K; B 4.680e-7 K; equal-B-energy C 5.602e-5 K;
ideal C 0.003506 K. The after-collapse pulse contributes exactly zero at that
earlier peak. The largest B increment later in the window is about 6.150e-6 K;
these two metrics are not mixed. The incremental pressure changes are small, but
the baseline mixture/polytropic pressure mismatches reach approximately 16%,119%,193%
at 1.0,1.15,1.3 bar. Thus even the baseline fixed mechanics is not self-consistent.

At 1.3 bar, the no-laser diagnostic actually stops near **18.3408990706 µs** at
6000 K; laser cases retain their own actual stop times. No coefficient-domain
continuation is supplied, and 6000 K is not claimed as a physical ceiling.
The baseline peak may lie beyond a laser-censored interval; its corresponding
increment is unavailable rather than extrapolated. No extra coupled check is
represented as validated in this stage.

## 5. Verification and important quantitative limits

- **31 automatic tests pass** in production and the fresh locked environment: geometry/contact, zero intensity/count/exchange, h=0,
  strong-contact monotonicity, energy and serialized-field consistency, diffusive
  SI scaling, passive invariant versus particle ODE, finite Gaussian normalization,
  signs, optical energy, vapor normalization, NASA bounds, causality, finite-design
  selection, source serialization/domains and upstream fingerprints.
- Fifteen mesh controls; production versus finest gas-energy difference <=0.196%,
  gold peak <=1.206%, facing liquid peak <=1.192% in selected cases. These empirical
  discrepancies are not rigorous global error bounds.
- Finest isolated Au gold-peak regression versus MR5 <=0.0893%.
- Independent Radau versus tight BDF: about 1.41e-9 gold-response difference,
  1.91e-12 relative gas-energy difference. Energy balance <=2.66e-14 of absorption.
- Independently integrated sampled wall/bath fluxes are retained alongside ODE
  ledgers. Gas increment 2 ns versus 0.5 ns cap changes selected maxima by <=8.54e-7
  relative. Gaussian energy normalization error <=2.22e-16 relative.
- Infinite-water finite-h point-sphere analytic limit is derived and checked.
  At a=20 nm,d=500 nm the finite-Au control differs by about 0.052% at 2 µs.
  Near d=1 nm,a=100 nm the center-point model badly misses the early finite-sphere
  channel. Exponentially small early far-source tails are explicitly unresolved.
- Domain/boundary tests: 1->2 µm padding changes reference gas delivery about 0.2%;
  0.5 µm padding is insufficient for that comparison. Bath-versus-water storage
  partition is domain dependent even where gas delivery has stabilized.
- The original same-machine full recalculation is retained as
  `clean_reproduction_original.json`. It did not establish cross-platform remesher
  determinism. The independent audit found 1405 stored versus 1404 regenerated
  nodes in event_1.00. The fixed test now uses the saved nodes/triangles/edges as
  numerical input and checks energy, flux and temperature invariants. No equality
  between freshly generated meshes is required. Three saved-mesh replays use a
  2e-5 max-norm observable tolerance (100 times the 2e-7 time-solver tolerance),
  separately from the existing spatial convergence study.
- `handoff_numeric_audit.csv` checks 59 rounded claims/bounds against original
  CSV/NPZ/API values. `correction_input_manifest.json` checks 55 original numeric
  files byte-for-byte and preserves all numerical arrays in 24 renamed sources.
- Fresh-environment correction checks are in `clean_reproduction.json` and the
  validation journal. They rerun tests, bounded saved-mesh controls, saved-data
  audits and presentation/PDF construction; no large parameter sweep is repeated.

Please review the actual `validation_summary.json`, convergence CSVs, gas increment
residuals, optical contract and event statuses rather than relying only on rounded
values here. The dispatch provides final exact file/checksum counts and remote
identity after publication.

## 6. Interface readiness and questions requiring independent judgment

`src/coupling.py` offers a stationary-domain solver and precomputed
`ThermalHistory`. The v2 source contract exports `q_wall_reference_baseline_W`,
`q_wall_reference_sum_W` and `q_wall_laser_increment_W`; it does **not** export
`q_wall_total_W`. The reference sum is the stored MR4 closure at explicitly
archived reference R(t),T(t), plus the frozen laser increment. In the retained
legacy diagnostics, reference T is the paired exported Tb+delta_T, not a future
supplied gas state. The reference sum is never called the full physical heat flux
for changed Tg. Use only `history.laser_increment_W(t)` in an ODE that already
evaluates its own current-state MR4 loss. A regression test changes gas temperature
and verifies that the stored baseline is not counted a second time. Legacy source
fields load as reference diagnostics; new exports use only explicit v2 names.
Time/geometry guards and formal applicability metadata remain in force.

A bounded analytic estimate now evaluates omitted liquid heating by prescribed
saved delta_Tg, without updating gas or mechanics. For an unobstructed stationary
sphere in infinite water, the surface step response is

`S(t) = h/(h+k/R) * [1-erfcx((1/R+h/k)*sqrt(alpha*t))]`.

Convolution with the saved gas increment estimates the extra positive return
`H*theta_surface`, H=4*pi*R²*h. The gas ODE already includes its own negative `-4*pi*k_g*R(t)*delta_Tg`
loss derivative; no second direct loss is added in this estimate. For B, estimated integrated return/source energy is about 0.2323% at
1.0 bar, 0.09284% at 1.15 bar and 0.000225% on the censored 1.3 bar interval.
The first-order response is left unchanged. A maximum-principle bound based on
`H*integral(max_so_far(delta_Tg))` is much looser: about 74.3%, 12.1%, 0.127%
respectively. Thus the small analytic estimate is not a rigorous error bound for
the original event geometry or a proof that feedback is universally negligible.
Finite-Au obstruction, domain truncation, changing R/advection and post-censoring
tails are excluded from this ideal control. `gas_feedback_estimate.csv` retains
all nine A/B/ideal-C controls, the already counted gas-loss term and refinement.

Please assess:

1. Conservative FEM weights, source/flux signs, field-storage consistency,
   geometry/refinement convergence and the limits of point/planar controls.
2. Adequacy of h and controlled baseline matching; the absolute unsteady liquid
   baseline is unresolved and the static mismatch exceeds the tiny gas effects.
3. The zero-gas-perturbation tagged channel, first-order gas linearization and
   omission of incremental gas-to-liquid feedback; conditions for safe MR7 use.
4. Pe/displacement and contact gates, passive invariant, unmodeled forces and the
   failure to claim a retained passive near-wall particle through strong collapse.
5. Actual local optical energy, near-boundary homogeneous-water Mie uncertainty,
   multiparticle independence and exclusion of the macroscopic layer budget.
6. Finite objective distinctions, wavelength plateaus, endpoint lag censoring,
   numerical-ranking uncertainty and all formal extreme failure labels.
7. NASA actual censoring times, fixed-mixture/ideal-EOS/transport limits, pressure
   mismatch, causal timing, separate incremental and baseline energy residuals.

The existing **Mini-Research 6 (Web)** was identified and read by its current
title. The verified MR5 handoff mechanism is the app's `send_message_to_thread`,
followed by reading the destination to confirm delivery. The final publication
dispatch will use that mechanism and will claim delivery only after acknowledgment.
Please return concrete corrections before any official closure. No RoadMap update
or subsequent research stage is authorized by this handoff.
