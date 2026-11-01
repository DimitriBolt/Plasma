# Source and model-convention audit

Accepted starting local/remote main: `0204ce7a8429406b7ec235b02cec8581274fc8c6`.
The remote was read with `git ls-remote`, not assumed current. All 1083 tracked
starting files are byte-compared to their accepted Git blobs and separately
SHA-256 inventoried. The original unrelated untracked files are not inputs or
publication candidates. No prior model, result, report or plan is rewritten.

The requested source files are present and readable: `Plasma-plan.tex` coupling
section; root and all stage READMEs; MR1-MR5 scientific report sources and validation
journals; MR2 `data_schema.md`, `admission.csv`, `cycle_comparison.csv` and selected
complete trajectories; MR3 `heat_transfer.py`, `saved_mechanics.py` and
`bubble_thermo.py`; MR4 `thermodynamics.py`, `model.py`, data schema, summary,
NASA coefficient archive and selected saved histories; MR5 `plasmonics.py`, model
details, data schema, parameters, review handoff, optical tables and raw Johnson data.
These are specific scientific inputs, not similarly named substitutes.

Important inherited distinctions verified in the sources:

1. MR2 C1 Hermite reconstruction is new interpolation; accepted nodes/exports do
   not recover the original dense integration polynomials. Full five-cycle
   completion is checked before event use. 1.3 bar is formal, not validated flow.
2. MR3/MR4 signed heat loss can heat the gas during expansion. MR4 saturated
   vapor is in the common NASA mixture and cannot be added again as constant pv.
3. MR5 optical signatures use `radius_m`, vacuum `wavelength_m` and n_medium=1.33
   for water. Its radial liquid bath has no bubble and its net-liquid energy is
   not gas energy. The 235.051 µJ layer scenario is not a local bubble source.
4. The 3.182 pJ MR5 moderate control used G=50 MW/(m² K); the default MR6
   G=100 MW/(m² K) finite-design numbers need not equal that value. Isolated
   regression explicitly matches the coefficients instead of forcing a result.

Primary external sources checked for the model context:

- Johnson and Christy, Phys. Rev. B 6,4370 (1972),
  https://doi.org/10.1103/PhysRevB.6.4370. The raw optical table already archived
  in MR5 is reused byte-for-byte; no replacement data or hot-metal fit is introduced.
- Baffou and Rigneault, Phys. Rev. B 84,035415 (2011),
  https://doi.org/10.1103/PhysRevB.84.035415; readable author's full text
  https://guillaume.baffou.com/publications/015-Baffou-PRB84_035415.pdf.
  Its diffusion, uniform-Au and interface context supports the inherited moderate
  parameter model, not an attained extreme temperature or bubble efficiency.
- Shewchuk's primary Triangle documentation,
  https://www.cs.cmu.edu/~quake/triangle.html and
  https://www.cs.cmu.edu/~quake/triangle.help.html. `p` fits the supplied straight
  boundary graph; `q20` refines triangle angles. No geometric exactness beyond
  chordal input boundaries or universal mesh-error theorem is claimed.
- NASA/TP-2002-211556, https://ntrs.nasa.gov/citations/20020085330.
  The actual species coefficients/domains are those fingerprinted in MR4.
  H2O limits the fixed mixture to 6000 K; no extrapolated accepted gas state follows.

The point-sphere Robin impulse formula is independently derived by transforming
the radial diffusion equation with w=r*T and inverting its elementary Laplace
transform. It is a separate analytic benchmark, not an approximation asserted
valid at d<<a. The local finite-particle model is axisymmetric, not angularly uniform.
No proprietary paper PDF is duplicated into MR6.
