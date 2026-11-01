# Saved data, units, signs and finite domains

All numeric API values are SI. CSV display amplitude is bar (100000 Pa);
figure axes use stated SI prefixes. CSV values use 17 significant digits;
round-trip parsing is recommended. NaN means unavailable, never zero.
All temperatures in the liquid model are **laser-induced rises**. Gas `T_K`
is absolute temperature. `baseline_T_K` is the common no-laser gas state;
`delta_T_K` is the directly integrated first-order increment.

- `reference.npz`: local 1 J/m², 532 nm, 20 nm Au radius, 100 nm surface gap,
  6 µm bubble, 10 ns top hat, h=4400 W/(m² K), G=1e8 W/(m² K). `time_s`
  starts at laser onset and ends at 2 µs. `gold_rise_K`,
  `nearest_water_rise_K` (liquid at the Au pole facing the bubble),
  `mean_Au_surface_water_rise_K`, `water_rise_K` [water node,time],
  `xy_m` [water node, cylindrical radius/height], `triangles`,
  `boundary_edges`, `boundary_markers` (1 bubble,2 bath,3 axis,4 Au),
  `wall_node_area_m2` and h. Temperature at the bubble wall differs from the
  named facing-Au water temperature. `q_wall_laser_increment_W` and
  `q_bath_W` are signed outgoing sink powers. `absorbed_so_far_J`,
  `gold_energy_J`, `water_energy_J`, `gas_energy_J`, `bath_energy_J`
  and `balance_residual_J` keep all channels distinct. `minimum_water_rise_K`
  records signed finite-element far-tail undershoots, without clipping.
- `gap_scan.csv`: d=1,2,5,10,20,50,100,200,500,1000 nm, fixed a=20 nm,
  R=6 µm, 10 ns, and unit 1e-15 J absorption. Fluxes/temperatures scale linearly
  to another local absorbed energy while coefficients and geometry remain fixed.
  Peak times are sampled, not claimed exact dense-output extrema.
- `fields/*.npz`: complete temperature/flux/energy time series and meridional
  fields for the near-contact finite sphere (a=100 nm,d=1 nm,h=1e6), far particle
  (a=20 nm,d=500 nm), and the 1.0/1.15 bar event geometries frozen at their tc
  radii. All use local F=1 J/m²,532 nm,10 ns. Keys match `reference.npz`, with
  JSON geometry/status metadata. Event `absolute_time_s` aligns laser center at
  tc-20 ns while geometry remains frozen at R(tc), as in the gas diagnostic;
  it is distinct from `event_snapshots.csv`, which freezes at each pulse center.
  These fields make Au/water/power curves and nonuniform wall flux auditable for
  several geometries, without claiming moving-boundary evolution.
- `finite_design.csv`: all 216 a x d/a x duration x wavelength combinations,
  with fluence continuously chosen under the declared 30 K/intensity/fluence
  allowances. h,G,k,R are fixed. Energies and peaks are scaled from unit
  responses. `gas_fraction` refers to 2 µs after pulse start.
  `incident_local_beam_J` uses a 25 µm equivalent beam radius at the local plane.
  `moderate_static_candidate` requires d/a>=5 and the specified equilibration
  screen. It does not admit changing mechanics or validate boundary optics.
  `design_extrema.json` selects separate objectives; within 1e-12 relative
  objective ties the largest absorption cross section is preferred. Thermal-cap
  wavelength plateaus are not evidence of a unique spectral optimum.
- Shared response summary fields: `energy_absorbed_J` (unit or physical scale
  as stated), `gas_energy_J` by 2 µs, `max_gold_rise_K`, `max_water_rise_K`,
  maxima times from pulse start, `max_q_W`, independent balance, nodes and
  geometric screens. `q_peak_at_window_end=True` means the gas-power maximum
  is not internally bracketed; the lag may exceed the response window.
  `q_peak_bracket_width_s` gives neighboring saved-sample spacing.
- `transport_sensitivity.csv`: one-at-a-time h=0 through1e10 W/(m² K),
  k=0.3/0.6/1.2 W/(m K), G=0 through1e12 W/(m² K), with all other parameters
  fixed to the reference. Unit absorbed energy is 1 fJ. These are separate
  scenario parameters; no global multidimensional optimum follows.
- `formal_regimes.csv`: near-wall a=20 nm,d=1 nm, 526 nm; 0.1 ps/1 ns/100 ns;
  F=0,0.1,1,10,100,1000 J/m². Flags enumerate the selected time, continuum,
  optical, capture, moderate-rise and material-reference failures. Formal rises
  are not actual gold/water temperatures. Nothing is filtered by those flags.
- `multiparticle.csv`: N is a **local count**, not m^-3. Centers lie on a ring
  normal to the illuminating beam, with no projected shadowing. It records
  source energy, single-pass extinction/beam allowance, center separation,
  a diffusion-tail overlap proxy and explicit separation screens. The envelope
  number density is only N/(4*pi*r_np²*2a), a placement descriptor, not a uniform
  random dispersion. N=0 supplies exactly zero absorption and gas energy.
- `timing_scan.csv`: pulse-center offset delta_t=t_laser-tc; windows are
  before tc and [tc-10,tc+10] ns. Static reference shape is translated exactly,
  with zero source before onset and no extrapolation after the saved response.
  Lags are relative to pulse center. `gap_timing_map.csv` repeats actual heat
  solves for each declared gap; no fitted Gaussian lag is used.
- `mechanical_events.csv`: full five-cycle completion/admission checked against
  MR2. tc selection is first strong R<R0/2 where present, else cycle-1 deepest.
  R<2Rc width and saved spacing are separate from the explicitly named MR4
  coupled thermal width/global thermal peak; these are not matched trajectories.
- `kinematics.csv`: sampled event-aware particle paths, surface gap and invariant
  r_np³-R³. Negative gaps are contact/intersection statuses, not heat-model states.
  The passive invariant is not expected for fixed-center or capture scenarios.
- `event_snapshots.csv`: fixed geometry at each pulse center for three motion
  assumptions. Contact-stopped rows have unavailable thermal fields. Pe and
  displacement use event-aware input nodes. Changes over diffusion time and over
  2 µs are separate. NaN diffusion-time change means that time exceeds the saved
  mechanical interval. All motion-failed results remain formal equation outputs.
- `gas_diagnostic.csv` and `histories/gas_<amplitude>_<mechanism>.npz`: bounded
  windows initialized from the named published MR4 baseline. Same saved R(t)
  for all mechanisms. Tmax is the **observed window** maximum, never a five-cycle
  global temperature. NASA boundary records `end_s`, actual censoring K and reason;
  `completed_window=False` means the requested tail is absent. Pressure mismatch
  is total mixture pressure / (old polytropic dry pressure + old constant pv)-1.
  The laser pressure increment relative to the common baseline is separate.
  Work is by the gas; MR4 heat loss is signed outward from gas. Source heat is
  positive inward. `signed_J`, `positive_J`, `negative_J` integrate total wall
  heat over the available window; negative part stays signed. Source energies
  cannot be compared as full-window totals after a censored stop.
- Gas arrays additionally retain `baseline_T_K`, `delta_T_K`, common+increment
  work/loss, laser energy, total residual and first-order incremental residual.
  The incremental first law uses C(T_base)*delta_T; the unmodeled quadratic
  energy is not silently treated as first-order transfer. `max_delta_T_K` and
  `delta_T_at_baseline_peak_K` are different metrics. The latter is unavailable
  if a laser-induced censoring event precedes that baseline peak.
- `histories/*_source.npz`: absolute time, signed baseline, laser increment,
  total W, JSON metadata with tc, finite interval, geometry, optical fluence,
  pulse center, model identity and validity status. These are first-order gas
  diagnostics using frozen liquid geometry. `ThermalHistory` never extrapolates.
  `refined_` histories and summary repeat the 1.15 bar increment with 0.5 ns
  maximum step; `gas_increment_convergence.csv` compares them with the 2 ns cap.
- Validation CSVs: `mesh_convergence`, `time_method_convergence`,
  `domain_convergence`, `isolated_MR5`, `analytic_point_controls`,
  `planar_spherical_comparison`, `mechanics_interpolation`,
  `gaussian_normalization`. All tested points remain, including early point
  tails that have meaningless large relative errors at negligible energy.
  The analytic point model does not include the finite Au sphere or its delay.
- `optical_contract.json`: actual MR5 optical API value and published moderate
  row, preserving its different G convention. `parameters.json` defines the
  production model/design. `environment.json` records actual versions.
  `input_inventory.json` fingerprints every tracked accepted upstream file.
  Validation and clean-reproduction records state empirical numerical scope.

All eight JPG figures and report numbers read these files. FILE_LIST includes
itself and SHA256SUMS; the checksum manifest hashes every listed file except
itself. LaTeX intermediate files, environments and caches are excluded.
