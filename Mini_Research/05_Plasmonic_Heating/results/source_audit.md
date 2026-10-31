# Source and interpretation audit

All four mandatory articles were obtained as complete readable PDFs and checked
against their titles, authors and DOI. SHA-256 fingerprints and full-text access
URLs are in `../data/provenance.json`. Third-party article PDFs remain outside
the publication package; the public-domain numerical optical table is included.

| Primary source | Direct use and boundary |
|---|---|
| [Johnson & Christy, 1972](https://doi.org/10.1103/PhysRevB.6.4370); [readable full PDF](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.6.4370) | Au column of Table I, 49 n/k rows; room-temperature evaporated films; relative epsilon computed as (n+i*k)^2. Film measurements do not establish the optical response of a heated, ligand-coated or very small sphere. |
| [Baffou & Quidant, 2013](https://doi.org/10.1002/lpor.201200003); [author-hosted full PDF](https://guillaume.baffou.com/publications/020-Baffou-LPR.pdf) | Physical heat-release framework, CW spherical control, interface jump, representative interface conductances, relaxation mechanisms and phase-change cautions. Volume is 2013; early online publication/copyright is 2012. |
| [Jain et al., 2006](https://doi.org/10.1021/jp057170o); [full PDF](https://pubs.acs.org/doi/pdf/10.1021/jp057170o) | Sphere-size/scattering trends and quantitative **calculated** comparison values, 40 nm diameter absorption at 528 nm and 80 nm diameter scattering at 560 nm. Their size-corrected epsilon/cubic interpolation differs from the uncorrected bulk/PCHIP baseline. |
| [Baffou & Rigneault, 2011](https://doi.org/10.1103/PhysRevB.84.035415); [author-hosted full PDF](https://guillaume.baffou.com/publications/015-Baffou-PRB84_035415.pdf) | Table I thermophysical reference constants, diffusion/interface scales, ~1.7 ps electron–phonon time, moderate-rise conduction context and short-pulse limitations. The implemented water finite-volume/matrix method is independently derived; it does not copy a printed optical formula without checking extinction/scattering definitions. |

The required raw table is copied byte-for-byte from the database at commit
`c5c2f188e848453def5970e347399d653df2ffc2`, file
`database/data/main/Au/nk/Johnson.yml`. The original CC0 header is retained.
The downloaded main-branch bytes were checked against this immutable version.
Its hash, date and direct immutable URL are archived. Raw wavelengths are μm;
the code converts to m. No rows are added or extrapolated. The numerical source
is sufficient for offline calculation; retrieving the articles is not a runtime
requirement after this audit.

Scientifically necessary supplemental sources:

- [miepython 3.0.2 official documentation](https://miepython.readthedocs.io/en/3.0.2/01_basics.html):
  independent implementation, host index and negative-imaginary sign convention.
  The comparison passes a conjugated index explicitly. The local implementation
  uses positive-imaginary loss and outgoing first-kind Hankel functions.
- [NIST ITS-90](https://www.nist.gov/si-redefinition/kelvin/kelvin-its-90):
  bulk gold point 1337.33 K, used only as a diagnostic, not a nanoparticle threshold.
- [NIST Au vapor-pressure data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440575&Mask=4):
  Stull Antoine A=5.46951, B=17292.476, C=-70.978 for pressure in bar and
  temperature in K, valid 2142–3239 K. Solving for one bar gives about 3233 K
  within the tabulated interval. It is not an actual destruction fluence or a
  model of evaporation. The fit is never extrapolated to room temperature.

Parameter separation: the densities/capacities/conductivities come from the stated
thermal table; real host indices, beam geometry, concentration lattice and laser
upper bounds are explicitly selected design parameters. The 30 K screen and
independence/overlap screens are conservative working assumptions, not published
universal material thresholds. All formula dimensions and loss conventions are
derived explicitly in the report.

Limited literature comparison does **not** constitute experimental validation.
The chosen four-source core does not provide a matched raw experimental thermal
trace/spectrum with all nanoparticle/coating/fluence conditions needed for such
validation. Only theoretical/estimated values are compared quantitatively; this
limitation remains an independent-review question.

Context from accepted studies is read-only: the MR1–MR4 reports/README and plan
were inspected, and their full tracked input inventory was verified against the
starting published commit. The 30 kHz acoustic frequency and the saved MR4
`M2_sat_air_1.15bar_coupled_eta1` temperature-half-height width (17.85146 ns)
provide only a timing illustration. No gas-temperature or bubble ODE is rerun.
