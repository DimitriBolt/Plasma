# Sources, scope and read-only interface audit

The entire accepted tree (1209 files) matches initial 493dbdc213f3789b7195d93760bd5c0c867cb98c.
GitHub connector fetch_commit(main) independently returned that same SHA; direct
git networking initially failed DNS. Source hashes and initial untracked status
are preserved in input_inventory.json. The prompt itself is outside the edited scope.

All specifically requested files exist and were read: plan Part G and H; MR1/2
mechanics, MR3 signed heat/gas and saved-mechanics definitions; MR4 model.py,
thermodynamics.py, tests, species_coefficients.json, parameters, summary and schema;
MR5 plasmonics.py, model_details, parameters and schema; MR6 coupling.py, tests,
model_details, schema, parameters, reference and four saved field archives. The
named MR4 and MR6 PDFs were text-extracted, including their equations and limits.

Read-only MR4 reuse: Mixture, NASA coefficient archive, capacity, sensible-energy
integral and Parameters. The complex directional extension uses the identical
piecewise NASA9 heat-capacity polynomial, selected by the real baseline T; it
is used for derivatives only. Accepted gas states are checked inside the original
species domain. Read-only MR5 reuse: optical_cross_sections, absorbed_power,
pulse_energy, ThermalProperties and thermal_scales. No static thermal efficiency
enters the moving RHS. MR6 source contract frozen-reference-baseline-plus-laser-v2
is respected: reference baselines/sums are never physical flux inputs. Level 1
remains the already published diagnostic. MR7 solves one current-state signed flux.

## Verified primary literature

- McBride, Zehe and Gordon, NASA/TP-2002-211556 (2002),
  https://ntrs.nasa.gov/citations/20020085330 and its full report.
  Coefficient format and thermodynamic basis; actual component domains are taken
  from MR4's archived data, not the database's overall 200--20000 K range.
- Baffou and Rigneault, Phys. Rev. B 84, 035415 (2011),
  https://doi.org/10.1103/PhysRevB.84.035415 ; complete author-hosted article:
  https://guillaume.baffou.com/publications/015-Baffou-PRB84_035415.pdf .
  Constant-property diffusion, one-temperature Au, interface and electron-phonon
  times are moderate-heating assumptions. This is theory, not experimental validation.
- Johnson and Christy, Phys. Rev. B 6, 4370 (1972),
  https://doi.org/10.1103/PhysRevB.6.4370 . The already accepted measured film table
  in MR5 is the numerical source, unchanged. No hot-Au or near-wall optical fit.
- Keller and Miksis, JASA 68, 628--633 (1980),
  https://doi.org/10.1121/1.384720 . The actual equation convention comes from
  the accepted MR2/MR4 derivation, including viscosity in the denominator.
- Kitz et al., Biomedical Optics Express 2, 291--304 (2011),
  https://doi.org/10.1364/BOE.2.000291 ; full primary experimental article:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3038445/ . It studies 90 nm diameter
  gold particles and 5 ns illumination, including substrates and suspensions.
  The reported 532 nm bubble-formation threshold is about 60 +/- 5 mJ/cm^2
  (600 +/- 50 J/m^2). This is a context point with different particle size,
  environment and observable, not a threshold or calibration for MR7. Their
  wavelength dependence also warns against absorption-only nucleation claims.

No parameter is calibrated against a desired Tmax. There is no matched acoustic,
laser, particle-placement and composition dataset for experimental validation.
A measured laser-induced nanobubble threshold does not validate the peak uniform
fixed-composition temperature of an acoustically driven six-micrometre bubble.

## Explicit next physical level

A defensible local hot-spot/material prediction needs moving axisymmetric liquid
transport with a resolved finite Au surface and a consistent flow, contact and
phase-change treatment. If the gas crosses its coefficient range, continuation
needs a separately formulated reacting EOS with formation/ionization energy and
species balances. NASA CEA data availability alone is not an energy-consistent
reaction model; equilibrium versus finite-rate chemistry must be decided against
nanosecond compression and acoustic crossing times. Neither extension is hidden
inside the present calculation. High-Mach flow additionally needs a compressible
liquid/shock model and nonsphericity tests; first-order KM is not certified by a
numerically converged high-temperature result.

Additional phase-reference source: NIST-JANAF liquid H2O at 298.15 K,
https://janaf.nist.gov/tables/H-063.html , Hf=-285.830 kJ/mol; gas H2O from the
accepted MR4/NIST value -241.826 kJ/mol. Their difference is a standard-state
enthalpy scale, not a high-pressure nucleation threshold. The reference internal
energy estimate subtracts Ru*298.15 for ideal vapor and negligible liquid pV.
No phase fraction, reaction extent or modified gas temperature is inferred.
