# Mini-Research 7: Integrated Laser–AuNP–Bubble Model

**Dimitri Bolt** · Academic supervisor: **Prof. Dr. Ildar R. Gabitov**.

[Scientific report](report/07_Integrated_Model.pdf) ·
[Numbered results notebook](07_Integrated_Model.ipynb) ·
[Mathematical decisions](results/model_decisions.md) ·
[SI data schema](results/data_schema.md).

This reduced model evolves Keller–Miksis mechanics, the accepted NASA9 fixed gas
inventory, positive-capacity material water volumes and dynamic Au storage.
The coordinate is xi=(r³−R³)/3, with r=(R³+3xi)^(1/3). One signed current-state
power Qlg=4πR²hg(Twall−Tgas) enters gas and water with opposite signs; it replaces
the old clamped-water loss. No frozen transfer fraction is an input.

The laser comparison uses the same moving-water model for both members. The
five-cycle peak difference is undefined if either run is censored. Finite nonlinear
pairs, unit-fluence derivatives, direct-difference numerical noise, mechanical
feedback and local particle motion are reported separately. The 1.3 bar cases
remain formal; 6000 K is the wet-gas coefficient boundary, not a physical ceiling.
Spherical Au-source averaging does not resolve a local water hotspot or guarantee
actual nanoparticle material validity. No fusion conclusion follows.

The package preserves all accepted MR1–MR6 files and the general plan. The work
is prepared for independent scientific assessment; official stage closure is a
separate decision. [Numerical validation](results/validation.md) and
[review handoff](results/review_handoff.md) contain the detailed evidence.

## Compact publication and continuation inputs

This repository contains the selected public package, rather than the roughly
1 GB local numerical archive. The exact [selection](results/publication_scope.json)
and [per-array inventory](results/public_data_inventory.json) state what is shipped.
All model source, parameters, aggregate tables, numerical checks, run metadata,
the report and figures remain available. No scientific value is refitted or rounded.

The [references](results/references/) retain complete accepted scaled states for
the no-laser baselines at 1.0, 1.15 and 1.3 bar (including the separately used
1.3 bar paired reference). These support future model comparisons through
`src.stable_pairs.archived`. Selected scalar traces retain the principal finite
nonlinear laser responses, formal dry extremes and all inputs needed to regenerate
the figures. Five liquid profiles are exact original snapshots. All 115 dynamic
runs retain their configurations and diagnostic summaries in [metadata](results/metadata/).

Other high-dimensional sweep outputs remain local. The full archive is recoverable
from the locally retained research snapshot and can be regenerated with the explicit
expensive script below; it is not an ancestor of the compact public commit.
Historical full-archive checks and current compact checks have separate scopes.

## Reproduce

Use Python 3.12 with `requirements-lock.txt`; the tests use standard `unittest`.
The sibling accepted stages supply read-only gas and optical inputs. LaTeX and
Poppler are system dependencies. Set `PYTHONDONTWRITEBYTECODE=1` and a writable
`MPLCONFIGDIR`, and run from this directory.
Keep the complete repository checkout: MR4 imports MR3 and MR2 transitively,
static controls read MR6, and the preservation test checks all accepted files,
including root files and Library. An MR7-only export is not standalone.

```bash
python -m unittest discover -s tests -v
python -m src.reference_verify
python -m src.presentation
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/07_Integrated_Model.tex
sha256sum -c SHA256SUMS
```

Figure/table regeneration reads saved data only. The expensive numerical replay
is explicitly separate:

```bash
bash reproduce_all.sh
```

Do not run duplicate experiment groups against the same output directory.
[`src`](src/) separates gas/mechanics, material transport, optics/particles,
experiments, cancellation-controlled paired differences and presentation.
All new scientific figures are print-quality JPG. The archived coefficients,
source inventory, environment and raw fields support offline reproduction.

For a bounded fresh numerical check, use `python -m src.replay_check`. It
reintegrates one no-laser cycle and compares a common prefix with the saved
five-cycle reference. It does not repeat the full parameter design. The full
replay script deliberately separates expensive integration from ordinary
figure/PDF regeneration and should run in a separate checkout.

The selected public scalar histories retain every original float64 sample.
`field_indices` maps the five shipped liquid profiles to scalar time. No-laser
references retain full accepted states for paired reconstruction. The earlier
[storage record](results/storage_compaction.json) describes the larger local
archive, not the smaller public subset. Set `MR7_DATA_MODE=full` only when using
or regenerating `results/raw`; that directory is deliberately ignored by Git.

For the nominal pulse, the stable Level 3 changes are about **1.41 µK at 1.0 bar**
and **2.18 µK at 1.15 bar**, conditional on spherical source averaging and an
unverified prescribed particle position. The late 1.15 bar tail is unresolved
under refinement. Static integrated transfers do not validate local Au heating.
