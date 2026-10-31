# Data schema and definitions

All calculations use SI units. CSV and JSON fields end in explicit units where
applicable. Model identifiers and paired reference run IDs define the comparison;
the model hierarchy is documented in README and `source_audit.md`.

## Trajectories

Each `trajectories/<run_id>.npz` contains `t_s`, `R_m`, `U_m_s`, `T_K`, `p_Pa`
(pressure of all inventory in the energy balance), `density_kg_m3`, `delta_U_J`,
`work_by_gas_J`, `signed_heat_loss_J`, `gross_compression_J`,
`energy_residual_J` and `temperature_stationary_times_s`. A legacy run's actual
mechanical interior pressure is `p_Pa + pv`; a mixture has no additional pv.
`pressure_relative_mismatch` is (new mechanical interior pressure / old
mechanical interior pressure - 1) for one-way runs; coupled runs have NaNs here.
The old denominator includes both its polytropic dry pressure and constant pv.

Thermal extrema are evaluated on solver dense output by derivative root finding,
not inferred from saved plotting points. The final series includes those refined
times. `<run_id>_accepted.npz` retains accepted time and integration state (T/T0,
scaled ledgers, and radius/velocity for coupled runs). Dense polynomials are
reconstructed by rerunning, not claimed to be recoverable exactly from the NPZ.

## Summary and censoring

`summary.csv` and per-run JSON agree. `observed_Tmax_K` is the refined maximum
on `[start_s,end_s]`. `global_Tmax_K` is present only for a complete requested
five-cycle interval. A coefficient or technical stop records its actual
`solver_message`, time and `boundary_temperature_K`. Missing global maxima are
NaN in CSV / null in JSON, not zero. `peak_cycle` indexes the acoustic interval
containing the temperature event, and is independent of the radius-minimum cycle.
An unbounded analytic upper-temperature domain is infinity in CSV and null in
JSON; the run loader restores that domain value on cached reproduction.
`t_Rmin_s` is the separately refined global radius minimum on the observed interval.

`p_at_Tmax_Pa` is the pressure used in the gas energy balance. For legacy runs,
`p_mechanical_at_Tmax_Pa` adds the separate constant pv. Formal ideal-gas density,
number density, maximum wall Mach, pressure/bulk-modulus proxy, thermal peak width,
and ideal acoustic crossing time are diagnostics, not admissibility filters.
The half-height width is relative to T0 around the selected thermal peak; it is
undefined if one half-height crossing is missing. No physical shock is resolved.

Near-peak candidate counts and first/last times use a declared relative band of
1e-8 below the maximum. This is a conditioning diagnostic, not an error bound or
a modification of the actual argmax. Weak periodic peaks can nearly tie and
numerical refinement can select another cycle. `peak_time_conditioning.csv`
keeps those changes explicit. CSV numeric reads use round-trip precision so that
decimal parsing does not perturb narrow-event times.

## Energy and contrasts

`delta_U_*_J` is sensible U(T)-U(T0) for the fixed composition. Work is by the
gas, heat loss is signed, and the source and chemical changes are zero in the ODE.
The residual is delta-U + work-by-gas + signed-heat-loss. Scaling uses the largest
absolute ledger/sensible-energy scale or n-total*Ru*T0 over the observed interval.
`gross_compression_*_J` integrates positive compression input, not signed work.
`collapse_compression_J` covers the preceding radius maximum to the thermal peak.
The independently integrated quadrature endpoints are the refined peak and the
observed end; 8/16 point rules are per accepted integrator subinterval.

`sensitivity.csv` defines global delta-T only when both paired trajectories finish
five cycles. `observed_max_ratio_diagnostic_only` compares observed maxima on
potentially different censored intervals; it is not a global physical effect.
`main_comparison.csv` uses the published MR3 table, rather than rounded prompt
values, as the reference. Regression and convergence retain signed relative
temperature differences, peak-time and observed-end-time differences. A small
thermal-boundary peak difference alone is not evidence for an uncomputed tail.

## Chemical assessment

Reaction enthalpies/internal energies are standard 298.15 K reference scales,
while atomic ionization thresholds are isolated ground-state energies.
`energy_for_complete_conversion_J` is a hypothetical inventory demand, compared
separately with each frozen run's energy/work. Full and first-step dissociation
channels overlap and must not be summed. Ionization of H/O/N presupposes their
hypothetical production by dissociation; no actual population is assigned.
The exponential energy diagnostic is not a fraction or rate. No chemical
temperature, equilibrium abundance or kinetics solution is stored or implied.
