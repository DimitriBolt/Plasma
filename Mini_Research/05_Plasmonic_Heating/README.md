# Mini-Research 5: Plasmonic Absorption and Photothermal Heating

**Dimitri Bolt**. Academic supervisor: **Prof. Dr. Ildar R. Gabitov**.

[English report (PDF)](report/05_Plasmonic_Heating.pdf) ·
[LaTeX source](report/05_Plasmonic_Heating.tex) ·
[numbered results notebook](05_Plasmonic_Heating.ipynb) ·
[model and interface details](results/model_details.md).

For an isolated Au sphere, `P_abs = I_local sigma_abs`, `E_abs = F_local sigma_abs`.
For a finite layer,

\[
E_{abs}=FA_b\frac{\sigma_{abs}}{\sigma_{ext}}
\left(1-e^{-N\sigma_{ext}L}\right)\leq FA_b.
\]

The passive dipole limit is compared with Mie theory using archived experimental
Johnson–Christy constants. A transient radial conduction model separates gold
storage, liquid storage, interface transfer and energy reaching the outer bath.
The standalone [MR6 interface](src/plasmonics.py) contains no bubble dynamics or
gas thermodynamics. MR1–MR4 remain byte-identical to the recorded starting state.

Within the declared 30 K gold-rise screen, a 100 nm radius sphere absorbs
**3.182 pJ** during a 100 ns pulse and transfers **3.137 pJ** to liquid by another
100 ns. A separately optimized dilute layer absorbs **235.051 μJ** from a
**785.398 μJ** pulse. Its absorbed fraction is 29.928%; the distinct fraction
maximum is 90.201%. These are conditional model results, not experimentally
validated operating points, a universal maximum or a gas-temperature prediction.
The formal extreme cases and their failures remain available in the data.

## Reproduction

Use Python 3.12 (verified: 3.12.14). From this folder:

```bash
sha256sum -c SHA256SUMS
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip check
export MPLCONFIGDIR=/tmp/plasma-mr5-matplotlib
python -m src.experiments
python -m src.validate
python -m src.presentation
python -m unittest discover -s tests -v
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/05_Plasmonic_Heating.tex
python -m src.package
```

The original raw table is already included; numerical rebuilding needs no network
after dependencies are installed. The sibling MR1–MR4 directories, `Library/`
and `Plasma-plan.tex` are read and hashed only. A Git clone verifies their base
objects; an exported archive verifies their saved inventory without requiring
`.git`. No old ODE is integrated. LaTeX uses the standard packages named in the
source; Poppler supplies PDF inspection. Eight figures are JPG only, 300 dpi.

`experiments` regenerates the optical cube, all finite-design tables and formal
scans. `validate` independently checks Mie, interpolation, spatial resolution,
thermal balances and attenuation quadrature. `presentation` reads saved results
to regenerate figures and numerical LaTeX inserts. The notebook reads saved data
and does not rerun the experiments. Notebook use optionally needs a normal
Jupyter environment, not required for the command-line reproduction.

Manifest checks apply to the published bytes **before** rebuilding; PDF creation
metadata and timing text may change on rebuild. Numerical tables are compared
separately. [Fresh-environment reproduction](results/clean_reproduction.json)
records the tested clean archive and numerical agreement.

## Evidence

- [Parameters with SI units](results/parameters.csv), [full design](results/parameters.json),
  [source audit](results/source_audit.md), [raw-data provenance](data/provenance.json).
- [Optical cross-section samples](results/cross_sections.csv),
  [spectral maxima](results/resonances.csv), [thermal responses](results/thermal_response.csv).
- [Single-particle design](results/bounded_single_particle.csv),
  [ensemble optima](results/bounded_ensemble.csv), [formal regimes](results/formal_regimes.csv),
  [ideal energy-transfer upper-bound scenarios](results/ideal_upper_bound.csv).
- [Numerical validation](results/validation_summary.json),
  [test log](results/tests.txt), [data schema](results/data_schema.md),
  [unchanged upstream inventory](results/input_inventory.json).
- [Independent-review handoff](results/review_handoff.md),
  [file list](FILE_LIST.txt), [checksums](SHA256SUMS).

Publication and successful numerical checks do not constitute independent
scientific acceptance. MR6 is not started here.
