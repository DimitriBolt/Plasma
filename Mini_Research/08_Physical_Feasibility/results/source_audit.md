# Primary-source and numerical provenance audit

All seven exact named supporting PDFs are present, readable and extracted in
full. Their hashes are in `source_inventory.json`; the complete preexisting
tracked tree is separately protected. Section names in the Final Report's
footnotes and guide were checked against the actual final TeX/PDF headings.
`Plasma-plan.tex` was read and preserved, including its original title/hypothesis;
the new reports use Mini-Research 1–8 numbering and distinguish objective from
demonstrated result.

## Principal original sources

- **Keller & Miksis (1980), DOI 10.1121/1.384720.** Complete six-page journal
  reprint retrieved from the [DTIC public archive](https://archive.org/download/DTIC_ADA095754/DTIC_ADA095754.pdf),
  including the original incident-field treatment and incompressible limit.
  The adopted local-in-time full-pressure-jump convention remains explicitly
  identified; it is not claimed identical to every original forcing convention.
- **Sadighi-Bonabi & Gheshlaghi (2016), DOI 10.1063/1.4945343.** Complete exact
  local file `Library/035218_1_online.pdf` read. Equations (1)–(15), laser
  description, Tables II–IV and the conclusion were checked. Table III's cold
  acetone hydrochemical and Table IV's laser conditions have different radii,
  frequencies/amplitudes. Their 1.6e6 and 1.9e6 K entries are literature
  simulations, not an isolated AuNP enhancement under matched conditions.
  The hotter-core discussion is not a directly measured core temperature.
- **Flannigan & Suslick (2005), DOI 10.1038/nature03361.** Complete exact
  `Library/nature03361.pdf` read, identifying the named article within pp. 52–55
  despite adjacent unrelated material in the file. Acid-host Ar/SO/ionic
  emission, effective excitation temperatures and the authors' distinction
  between emission and a proposed hotter core are preserved. No measured
  thermonuclear ion temperature is attributed to this work.
- **Xu & Suslick (2010), DOI 10.1103/PhysRevLett.104.244301.** Exact local
  `Library/PhysRevLett.104.244301.pdf`; phosphoric-acid/molecular thermometry is
  contextual, not a quantitative calibration of the water model.
- **Denner (2021), DOI 10.1016/j.ultsonch.2020.105307.** Complete author manuscript
  [arXiv:2008.04332v2](https://arxiv.org/abs/2008.04332v2) read, especially
  equations (10)–(25), the implementation appendix and water Tait reference.
  Only its Gilmore/barotropic mechanical formulation is implemented here.
  The published Tait exponent 7.15 is retained; B is derived to match the
  existing rho0/c0, rather than mixing two ambient parameter sets. No
  Gilmore–NASG temperature or validated Tait liquid temperature is claimed.
- **NASA/TP-2002-211556.** Complete
  [NASA PDF](https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf)
  read together with the archived MR4 CEA species coefficients and metadata.
  General database coverage to 20000 K is not confused with the narrower
  6000 K domain of the selected H2O entry. No coefficients are extrapolated.
- **Baffou & Rigneault (2011), DOI 10.1103/PhysRevB.84.035415.** Complete
  [author-hosted paper](https://guillaume.baffou.com/publications/015-Baffou-PRB84_035415.pdf)
  read: Table I, electron/phonon relaxation, conduction and single-temperature
  limits. The inherited Au/water constants and 1.7 ps scale are references;
  finite-domain discretizations and MR8 audits are independently specified.
- **Kitz et al. (2011), DOI 10.1364/BOE.2.000291.** Complete
  [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3038445/) retrieved and
  read: 90 nm diameter, 5 ns excitation, supported single-particle threshold
  around 600 J/m² at 532 nm. No threshold is transferred to a 20 nm-radius
  particle in a moving acoustic geometry. Earlier static thresholds remain
  diagnostic context, not material admissibility certificates.
- **NIST:** [Au ITS-90 point](https://www.nist.gov/si-redefinition/kelvin/kelvin-its-90),
  [Au vapor-pressure data](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440575&Mask=4),
  [Ar ground-state ionization](https://physics.nist.gov/PhysRefData/Handbook/Tables/argontable1.htm)
  and archived MR4 formation-energy inputs. The bulk Au reference is 1337.33 K;
  the vapor-pressure fit is valid only over 2142–3239 K. Reference costs are
  not converted into equilibrium fractions or high-density temperature caps.
- **Nagy & Hegedus (2025), DOI 10.1016/j.ultsonch.2025.107651.**
  [Full article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12641221/) checked for
  convention and limited DNS-comparison context. Its benchmark ranges are
  not universal acceptance thresholds for this system.
- **Toegel et al. (2000), DOI 10.1103/PhysRevLett.85.3165.** The complete exact
  published four-page PDF was found in the preserved local primary-source cache
  `/tmp/mr4-audit/toegel.pdf` and read in full. Its vapor trapping, exchange and
  heat-capacity discussion supports the qualitative distinction from a fixed
  vapor inventory. Fresh publisher/university downloads returned 403, so the
  readable exact cached article was used; no excerpt or different paper was
  substituted. No numerical result from it is imported into the present runs.
- **D–D reactivity:** [UKAEA PROCESS official scientific documentation](https://ukaea.github.io/PROCESS/source/reference/process/models/physics/fusion_reactions/)
  checked for branch-specific Bosch–Hale parametrization context and ion
  temperature dependence. No reaction coefficient is evaluated or extrapolated,
  since isotope density, ion state, volume and lifetime are not supplied.

Third-party papers downloaded for inspection remain outside the publication
package. Their access URLs and fingerprints are recorded in the external-source
inventory; the published coefficients and all numerical continuation inputs
needed by this package are already in the repository. The inaccessible original
Gilmore 1952 PDF is not represented as read: the new implemented equation is
derived from the explicitly identified complete Denner primary manuscript.

## Numerical source map

The MR3 comparison and formal-extremum rows come from its `thermal_summary.csv`
and `adiabatic_summary.csv`; MR4 values from its final report/tables and coefficient
archive; MR5 optical/thermal values from its final report and unchanged source
table; MR6 resolved fractions/locality from `results/report_values.tex`; MR7
increments from `stable_pairs.csv`, diagnostics from `physical_diagnostics.csv`,
timing from `transport_timing.csv`, and energies directly from the selected
public reference NPZs. New results are labeled MR8. Rounded report values are
not used as new solver inputs.
