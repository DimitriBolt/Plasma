# Source and convention audit

All three specifically named local publications, the general plan (especially Part C),
and the accepted MR1/MR2 packages were present and readable before new model code was
written. The full local PDFs were extracted with Poppler; relevant equation and
spectroscopy pages were also rendered and visually inspected. No missing named
publication was replaced by its abstract. `input_inventory.json` hashes the unchanged
upstream package and verifies its bytes against MR2 commit
`f1f3e812d77fd7fa058ce111ff484b1c8658f937`. MR1's historical numerical baseline is
`d59bc0d0119b514152f763e2537a0477c7c63663`; authorized documentation changes already
present in the accepted MR2 snapshot are preserved.

## Mechanics and thermodynamics

- `Plasma-plan.tex`, Part C, specifies the compression/energy/heat-transfer progression.
  Its schematic pressure normalization `p0` is replaced here by the required equilibrium
  noncondensable-gas pressure `pg0=p0+2 sigma/R0-pv`, as in accepted mechanics. The plan
  itself is unchanged. No later part is implemented.
- MR2 `src/model.py` fixes the adopted KM pressure prefactor `(1+U/c)` and the derivative
  of the complete wall-pressure jump. Its viscous derivative contains both
  `-4 mu A/R` and `+4 mu U^2/R^2`. MR3 retains this convention and substitutes the
  energy-consistent derivative of gas pressure; it does not regenerate old mechanics.
- Sadighi-Bonabi & Gheshlaghi (2016), DOI
  [10.1063/1.4945343](https://doi.org/10.1063/1.4945343), local
  `Library/035218_1_online.pdf`: Eq. (1) prints a different pressure prefactor and
  Eq. (4) uses a different forcing phase. Their excluded-volume quasi-adiabatic model
  and hydrochemical model are not the closures implemented here. Table II's water
  case uses R0=6 micrometres, f=30 kHz, Pa=1.3 bar, T0=20 C and reports 1.0e4 K;
  Table III's water hydrochemical case reports 2.7e4 K, with argon and host vapor
  participating in the model. Their liquid constants differ from this study.
  The deuterated-acetone hydrochemical result 1.6e6 K is associated with -28.5 C,
  R0=8 micrometres, 38 kHz, and Table III's 14.4 bar; nearby text uses atmospheres,
  an internal units inconsistency that must not be silently reconciled. Table IV's
  laser-modified model uses other radii, forcing and laser energy. The proposed
  1e7 K core extrapolation is not a measured uniform gas temperature and is not
  a calibration target or evidence of fusion.
- Zhou & Prosperetti (2020), DOI
  [10.1017/jfm.2020.645](https://doi.org/10.1017/jfm.2020.645), readable
  [author institutional full text](https://ris.utwente.nl/ws/files/240998140/Zhou2020modelling.pdf).
  Sections 2--3 distinguish nearly uniform pressure from a spatial temperature
  distribution. Their quadratic/biquadratic closures retain spatial information.
  Sections 5--6 show that thermal exchange changes damping and radius histories.
  MR3's single-temperature resistance is a simpler phenomenological approximation;
  its eta is not fitted to their model and their validation does not transfer to it.

## Conductivity and composition

Marcia L. Huber, *Thermal Conductivity of Gases*,
[NIST-hosted table, page 1](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=926438),
gives at 300 K: air 26.4 mW/(m K), argon 17.7 mW/(m K). The air entry is at 100 kPa;
the argon entry is the dilute-gas limit. We use these **300 K constants**, without
pretending they are exact properties at 293.15 K or at extreme collapse states.
Gamma=1.4 and 5/3 represent idealized diatomic and monatomic heat capacities. Effective
gamma=1.2 and 1.55 are mathematical sensitivity parameters, not chemical identifications.
The adopted molar masses 0.02897 and 0.039948 kg/mol are bookkeeping assumptions for
mass density. The gas amount is derived from the initial pressure, volume and temperature,
not fitted; composition at the same initial state changes Cv, conductivity and mass,
but not the number of moles. Changing T0 while holding water and mechanics fixed is
explicitly a conditional sensitivity calculation.

## Experiments: temperature definitions matter

- Flannigan & Suslick (2005), DOI
  [10.1038/nature03361](https://doi.org/10.1038/nature03361), local
  `Library/nature03361.pdf`: the named article occupies journal pp. 52--55;
  the PDF also contains neighboring articles. In 85 wt% aqueous sulfuric acid,
  atomic Ar emission gives approximately 8000, 10700 and 15200 K at reported
  2.3, 2.5 and 2.8 bar. Because the bubble moves away from the antinode, those
  hydrophone pressure estimates are upper limits to the pressure it experiences.
  Ar/Ne mixtures at 3 bar exhibit strong conductivity dependence. Emission
  temperatures need not describe the opaque core or a uniform bubble average;
  the authors discuss nonthermal excitation and continuum-fit limitations.
- Xu & Suslick (2010), DOI
  [10.1103/PhysRevLett.104.244301](https://doi.org/10.1103/PhysRevLett.104.244301),
  local `Library/PhysRevLett.104.244301.pdf`: in 65% phosphoric acid regassed with
  50 torr Ar, OH rovibronic emission gives 6200--9500 K as acoustic pressure
  rises from 1.9 to 3.1 bar. At 2.4 bar, the reported fit is 8300 +/- 200 K.
  The bubble moves in the resonator. Noble-gas substitution changes emission
  temperature, conductivity and spectral intensity; water vapor decomposition
  is a limiting process. These are emission fits, not direct measurements of
  the uniform noncondensable-gas temperature used here.

No water/acid/acetone numerical match is claimed without matching host, gas, forcing,
inventory and closure. High T alone establishes neither plasma nor nuclear reactions.
