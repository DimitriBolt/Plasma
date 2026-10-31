# Final verification record

The completed checks establish reproducibility and internal consistency of the
specified models. They do not establish independent scientific acceptance or
experimental validation. No MR6 implementation was added.

- The starting remote main was independently verified as
  `803079ef68c0d3b02514dc88e211892079beeaff`.
- All 1008 tracked inputs from MR1–MR4, Library and the project plan were compared
  byte-for-byte with that commit and SHA-256 inventoried. No accepted model,
  data or result was modified. Pre-existing unrelated user files were preserved.
- All four required primary papers were fully readable, examined and
  fingerprinted. The original 49-row Johnson.yml table is included and checked
  against its immutable database version. Its omission from an initial packaging
  whitelist was found by clean export testing and fixed before publication; the
  packager now requires this data file explicitly.
- All 173088 optical grid points have positive absorption and satisfy
  extinction = absorption + scattering. No result is clipped. Lossless and
  index-matched controls, the dipole limit, full-grid series refinement,
  independent miepython comparisons and interpolation checks pass.
- Piecewise-linear/PCHIP sensitivity remains separate from numerical convergence.
  Competing spectral maxima and all formal high-fluence cases remain in the data.
- Depth integration, optically thin and opaque mathematical limits and the
  finite-beam energy partition pass. The independent opaque-layer quadrature
  explicitly resolves absorption lengths, avoiding an unresolved boundary zone.
- Thermal checks cover 90/180/360 against 720 water cells at 20/80/100 nm,
  0.1/1/10/100 ns and three interface conductances; a larger outer bath, independent
  BDF integration, analytic steady/short-pulse limits and independently integrated
  interface flux provide separate controls. Water and bath storage are distinguished.
- All 18 automatic tests pass in the production and fresh locked environments.
  `tests.txt`, `clean_environment_tests.txt` and `validation_summary.json`
  preserve the actual outcomes.
- A new virtual environment was installed from requirements-lock.txt. In a
  separate exported copy without .git or the development environment, the full
  optical, thermal, finite-design and formal scans plus numerical validation were
  rerun. All 19 CSV tables, all arrays in the two NPZ files, numeric JSON results,
  all eight JPG figures and the extracted report text agree exactly. The PDF
  recompiles successfully. See `clean_reproduction.json`; same-machine exact
  reproduction does not imply cross-platform byte identity.
- Every page of the nine-page final PDF was rendered to JPEG and visually
  inspected: pages 1–7 scientific text/figures, page 8 verification appendix,
  page 9 references. Equations, authorship, captions, colorbars, tables and page
  boundaries are readable; no clipping or overlapping text was found. The energy
  schematic's border margins were corrected and page 1 reinspected. LaTeX has
  no overfull box or unresolved references; a harmless underfull bibliography
  line reflects URL wrapping. PDF metadata names Dimitri Bolt.
- The selected maxima are explicit design-grid results, checked against saved
  candidate tables. The one-particle and finite-layer objectives are distinct.
  The 30 K screen, finite intensity/fluence, geometry, spacing, scattering and
  overlap assumptions are listed. Constant coefficients, cold bulk optics and
  lack of experimental calibration remain material limitations.
- The scientific report has no Git SHA, Codex instructions, agreement-status
  labels or internal correspondence. It does not use Part E as its title or
  numbering. Publication identity, manifests and dispatch evidence are kept
  outside the scientific text.

The final packaging step regenerates FILE_LIST.txt and SHA256SUMS; SHA256SUMS
hashes every listed file except itself. The staged Git archive is additionally
checked for file-list equality, required data/PDF/code/figure presence, manifest
integrity and test success before commit and publication. Actual publication
SHA and delivery acknowledgment are reported by the dispatch, not guessed here.
