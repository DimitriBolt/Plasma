# Saved thermodynamic results (SI)

`bar` in amplitude columns is the only non-SI input display unit: 1 bar = 100000 Pa.
CSV floats retain double precision. Readers needing narrow-event timing should use
`pandas.read_csv(..., float_precision='round_trip')`. NaN/null means unavailable,
never a zero. NPZ arrays are compressed double precision. Their time grids are
event-aware and nonuniform; no PDF/plot sampling is used to define extrema.
For eta=0 the relaxation time is infinite (CSV `inf`, JSON `null` because JSON has
no infinity literal); this is identifiable from eta and does not mean zero time.

- `input_inventory.json`: paths relative to the repository root, SHA256, byte size,
  and verification against the accepted MR2 snapshot. Includes all tracked MR1/MR2
  inputs, Library files and the plan. These are external inputs, not copied here.
- `input_admission.csv`: the upstream admission record plus completion, observed
  final time, solver status, stopping reason, extrema scope and accepted-node count.
  `diagnostic_partial_interval_only` is never promoted to a full five-cycle result.
- `adiabatic_summary.csv`: both mechanics at all 56 published amplitudes. `observed_Tmax_K`
  applies to the available interval; `global_Tmax_K` exists only if `completed=True`.
  `Tmax_from_refined_Rmin_K` uses the authoritative upstream refined minimum directly.
  Its difference from the maximum of a new Hermite interpolant is explicitly retained.
  Upstream wall Mach/pressure proxy and admission/numerical categories are preserved.
- `A_{RP,KM}_X.XXXbar.npz`: selected air-like Model A time, T, pg, volume and gas
  density. Radius/velocity are read from MR2 and are not duplicated in these exports.
- `collapse_temperatures.csv`: every selected radius minimum/maximum, temperature,
  acoustic cycle, and separate flags for first strong and global deepest. First strong
  means the first local minimum R<R0/2. No such event exists at 1.0 bar.
- `RP_KM_cycle_temperature_comparison.csv`: inherited independently selected cycle
  extrema plus their Model A temperatures. Sharing a cycle is not event matching.
- `adiabatic_sensitivity.csv`: KM trajectories at 1.0, 1.15 and 1.3 bar for gamma
  1.2, 1.4, 1.55, 5/3 and T0 273.15, 293.15, 313.15 K. Water constants and mechanics
  are fixed; n is recomputed by equilibrium. The changed reconstructed pressure does
  not equal the old mechanical closure when gamma differs from 1.4.
- `thermal_summary.csv`: 30 Model B runs = 3 amplitudes x 2 coupling choices x
  (4 air-like eta values + 1 argon-like case). The source is always 0 W. Gas name,
  gamma, T0, conductivity and eta identify each row. `B_one_way_post_processing`
  is intentionally named; `B_fully_coupled` evolves (R,U,T). The extrema scope,
  completion and solver message must be read with all maxima and durations.
- `{one_way,coupled}_KM_...json`: full parameters/settings, source commit, summary
  and every interval above 1000, 10000, 100000 and 1000000 K. Each interval is [a,b]
  in seconds on the observed domain. Counts and durations cover all cycles.
- Corresponding `.npz`: t_s, R_m, U_m_s, T_K, pg_Pa, V_m3, density_kg_m3, Qloss_W,
  work_J, signed_heat_loss_J, source_energy_J. R/U for Model B define the reconstructed
  or newly coupled path used by that calculation. `*_accepted.npz` contains the
  thermal solver's accepted times and states, not dense polynomials. For one-way,
  states are (theta,W/E0,Qloss/E0,Qsource/E0); coupled states are
  (R,U,T,W/E0,Qloss/E0,Qsource/E0). E0=nCvT0. No later interpolation is claimed exact.
- `adiabatic_hot_durations.csv`, `thermal_hot_durations.csv`: cumulative hot duration
  in seconds and hot volume-time exposure in m^3 s, not an emission pulse width.
  Entire bubble volume is hot under the uniform-temperature assumption. Zero duration
  means threshold not attained, not a fabricated censored tail.
- Gas pressure, volume and mass density at the temperature peak need not coincide
  with their own global extrema. `max_pg_Pa` is sampled on the event-aware grid.
  Model B wall Mach and pressure/bulk-modulus proxy are also sampled diagnostics.
- `t_Tmax_s` is found from the actual temperature function and derivative roots.
  `t_Rmin_s` is the independently selected global radius minimum. Their difference
  is not automatically a local lag or a comparison of the same collapse.
- `max_abs_pressure_relative_mismatch` is max abs(pg_new/pg_old-1) for one-way.
  `pressure_mismatch_at_Tmax` is the signed ratio minus one at the temperature peak.
  Old pg is the unchanged gamma=1.4 mechanical pressure, including for argon-like
  post-processing. Large mismatch flags inconsistency, not a mechanical correction.
- `signed_heat_loss_J` is integral Qloss dt; it can be negative. Positive is net
  loss by the gas. `work_by_gas_J` is integral pg dV; positive is work by the gas.
  `energy_residual_scaled` divides the largest absolute gas-energy residual by
  max(E0, max|DeltaUg|, max|W|, max|Q|). Zero source contributes exactly zero.
- `solver_convergence.csv`: coarse, production, tight DOP853 and independent Radau
  comparisons. All variants get production/tight checks; the two eta=0/1 air-like
  regimes at all three amplitudes also get Radau. Differences are observations,
  not rigorous error bounds or physical uncertainties.
- `adiabatic_recovery.csv`: tight coupled eta=0 compared with saved mechanics and
  with the analytic temperature invariant. No original mechanical sweep is rerun.
- `input_grid_convergence.csv`: full export, every second/fourth export sample and
  accepted nodes. Authoritative radius extrema and final endpoints are added to all
  reconstructions; these data contain information beyond the coarsened nodes.
- `interpolation_checks.csv`: maximum withheld-export radius and velocity differences.
  Hermite velocity is the radius derivative, not an unrelated interpolation of U.
- `gauss_energy_*`: independent 16-point Gauss quadrature per accepted step, separate
  from integrated ledgers. Inward/outward heat are nonnegative magnitudes used only
  to describe bidirectional exchange; the physical flux is never clipped.
- `parameter_provenance.json`, `source_audit.md`: assumptions, constants, exact
  conductivity reference temperature and links, composition and source interpretation.
- Report tables/macros and all JPG figures are generated from these saved outputs.

Formal RP results, effective-gamma frozen-mechanics sensitivities and all extreme
states remain conditional equation outputs. No measured temperature or experimental
validation is present. No arbitrary T/Mach/radius filter suppresses an extremum.
