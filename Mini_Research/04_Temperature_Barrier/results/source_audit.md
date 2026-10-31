# Sources and upstream audit

The starting local and remote `main` were both
`1a12fac295dc904569ab2a57ce4f68a8ea75f0c4`. The accepted MR3 scientific snapshot
`6dd2a7b034a1361a82aa44507dd73688b8b69523` is its parent; the latest change only
simplifies an introductory paragraph. All 611 tracked MR1--MR3, Library and plan
inputs match the starting commit byte for byte (`input_inventory.json`). The
preexisting untracked user files are inventoried and are outside this package.
No accepted code, data, figure, report or plan is edited. Publication completeness
is established by matching the remote commit and the tracked package inventories.

## Local required documents

All three named local PDFs are present, readable and extracted with Poppler.
`Plasma-plan.tex`, the Temperature-Barrier Physics section, motivates vapor,
excitation and ionization. Its suggestions of a reduced temperature and a compact
Saha equation are not assumptions imposed on the answer. The explicit current
specification takes priority: frozen-mixture dynamics plus quantitative chemical
energy assessments, with no mass exchange or equilibrium chemistry.

- Sadighi-Bonabi and Gheshlaghi (2016), `Library/035218_1_online.pdf`,
  [DOI](https://doi.org/10.1063/1.4945343), Sec. II and Tables I--IV: their
  excluded-volume, transport, liquid constants and pressure convention differ.
  Their hydrochemical water result is not a calibration target for this ideal
  mixture. Their acetone and laser scenarios are outside the implemented scope.
- Flannigan and Suslick (2005), `Library/nature03361.pdf`,
  [DOI](https://doi.org/10.1038/nature03361), journal pp. 52--55 (the PDF contains
  neighboring articles): emission measurements in concentrated sulfuric acid
  concern a moving bubble and species-dependent excitation. They do not measure
  this model's uniform temperature in water.
- Xu and Suslick (2010), `Library/PhysRevLett.104.244301.pdf`,
  [DOI](https://doi.org/10.1103/PhysRevLett.104.244301), four pages: OH emission
  in phosphoric acid motivates water-vapor chemistry, while its liquid, drive
  and spectroscopic temperature differ from the present conditions.

## Required primary research

Full PDFs were retrieved and their model formulations examined, rather than
substituting search excerpts. Retrieval hashes are in `retrieved_sources.json`.
Publisher full-text download links were accessible even when DOI/harvest requests
failed. Third-party summaries are not coefficient sources.

- Toegel et al. (2000), [full text](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.85.3165),
  pp. 3166--3167, Eqs. (2)--(4) and the internal-energy/heat-capacity expressions:
  vapor trapping and molecular vibration affect compression. Their mass-variable
  open-system enthalpy term is essential to their model; it is not copied into a
  closed fixed-inventory model. No latent-heat term is invented here.
- Yuan (2005), [author full text](https://lsec.cc.ac.cn/~lyuan/pre2005.pdf),
  Secs. II and V: compressible multicomponent gas dynamics, interfacial transfer,
  reference formation energies, vibrational excitation and finite chemical rates
  distinguish its results from a frozen uniform gas. Its chemistry can lower the
  temperature, but this is not a fitted constraint on MR4.
- Hammer and Frommhold (2002),
  [full text](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.65.046309),
  Secs. II--IV, Eqs. (8), (15)--(20): vapor exchange, chemical energy, varying
  particle number, ionization and reverse processes must be treated consistently.
  Their rate equations and density effects are reasons not to interpret an
  isolated energy subtraction as equilibrium or kinetic ionization.

## Thermodynamic coefficients and ranges

McBride, Zehe and Gordon (2002), NASA/TP-2002-211556,
[NASA record](https://ntrs.nasa.gov/citations/20020085330),
[full report](https://ntrs.nasa.gov/api/citations/20020085330/downloads/20020085330.pdf),
pp. 2--3 and Appendices B--D, defines NASA9 cp/R and integrated h/(RT), s/R.
The N2, O2 and H2O records were checked against the Appendix D listing.
The coefficient numbers are also distributed by NASA's own
[CEA repository](https://raw.githubusercontent.com/nasa/cea/main/data/thermo.inp).
The retrieved source includes later annotations; the three heat-capacity records
used for the dynamics retain the published numbers. Exact raw species records,
their units, intervals and the hash of the full downloaded file are saved.

N2 and O2: 200--1000, 1000--6000, 6000--20000 K. H2O: 200--1000 and 1000--6000 K.
No accepted variable-capacity trajectory crosses these domains. Internal trial
stages may evaluate the nearest polynomial solely to bracket a terminal event;
they are not saved as an extension of the physical model. Analytic integration
of Cv, with interval offsets, makes sensible energy exactly continuous despite
small cp fitting discontinuities. Formation-energy constants cancel for fixed
composition and are not counted as extra sensible energy. Dry air is an explicitly
idealized 79/21 mole-percent N2/O2 mixture, with no chemical dissociation.

Argon uses the ground-state ideal monatomic identity cp=5Ru/2, cv=3Ru/2, not
an artificial polynomial dependence. Its analytic continuation to high T is a
formal frozen-electronic-state model, not a claim of physical validity. NASA's
high-T atomic excitation fit is not used or extrapolated. Ionization and electronic
excitation remain absent from the ODE. The constant-capacity controls evaluate
NASA cp at 293.15 K and then define a separate calorically perfect model at all
positive T; they are not high-temperature extrapolations of the NASA polynomial.

MR3 conductivities are retained identically in all same-gas comparisons:
air 0.0264 and argon 0.0177 W/(m K), from the 300 K entries of
[Huber's NIST table](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=926438).
They are controlled parameters, including in wet mixtures, not validated
temperature-dependent mixture transport coefficients. Water properties and pv
are inherited from MR3. [IAPWS-95](https://iapws.org/technical-guidance/release/IAPWS-95)
is the water reference, not a new liquid EOS calculation at collapse conditions.

## Chemical energy inputs

[NIST Chemistry WebBook](https://webbook.nist.gov/) gas formation enthalpies at
298.15 K: H=217.998, O=249.18, N=472.68, OH=38.99, H2O=-241.826 kJ/mol;
the atomic values and water use the indicated CODATA reviews,
OH uses Chase (1998). N2/O2 are standard reference elements. Individual URLs
are saved in `chemical_parameter_sources.json`. Reaction enthalpies and
delta-U=delta-H-delta-nu*Ru*T are computed, not mislabeled zero-K bond energies.
These are fixed reference energy scales, not activation barriers or complete
high-temperature reaction enthalpies. The NASA tables for H/O/N/OH are archived
only as provenance for possible verification, not used to evolve chemistry.

NIST atomic handbook first ground-state ionization energies: Ar=15.759610,
H=13.598433, O=13.61805, N=14.5341 eV. These are isolated-atom thresholds,
not density-corrected in-bubble potentials. Exact links are saved with the inputs.
Complete-conversion energies are compared separately to available frozen-model
energy and work. Mutually dependent dissociation channels must not be added.
No reaction extent, LTE fraction, kinetic rate or chemically corrected Tmax is
claimed. The exp(-E/kT) diagnostics do not contain partition functions or density
factors and are explicitly neither equilibrium populations nor rates.
