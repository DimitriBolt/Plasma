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

## Numerical completion and physical applicability

The reported finite nonlinear increments retain their original numerical values:

| Amplitude | Five-cycle change | Numerical completion | Physical applicability |
|---|---:|---|---|
| 1.00 bar | +1.4088163 µK | Complete, five cycles | Conditional reduced-model result; prescribed AuNP retention at a 100 nm gap is not physically validated |
| 1.15 bar | +2.1848791 µK | Complete, five cycles | Formal fixed-composition scenario; not a validated physical gas state |
| 1.30 bar, wet | undefined | Censored at 6000 K gas-coefficient boundary | Formal, incomplete trajectory; no global five-cycle difference |

For 1.15 bar, prescribed far-field absolute pressure becomes negative during part
of each period, gas pressure reaches approximately **1.30×10⁸ Pa**, and fixed
composition omits dissociation and ionization. Phase changes and the applicability
of first-order mechanics and incompressible liquid transport remain unresolved.
Negative liquid pressure alone is not a mathematical proof of impossibility:
metastable tension, cavitation and its stability are separate physical questions.
Completing numerical integration does not resolve them.

In the matched static comparisons, MR7 predicts peak Au temperature rises of only
**21–30% of resolved finite-Au MR6 values**. Integrated gas-energy delivery differs
by about 0.04–1.68%; that agreement does not validate the local Au temperature
field. Melting, boiling and nearby phase-transition thresholds in MR7 are
diagnostic references and cannot establish physical admissibility. The static
discrepancy is not a universal error bound for the moving model.

The late event-centered 1.15 bar response changes by **89.6%** from 48 to 96 cells
and remains **unresolved**, distinct from the stable early response. Wet 1.30 bar
censoring is preserved; 6000 K is not a physical ceiling and no missing tail is
filled. Dry formal baselines remain 6263.579 K (air) and 22808.716 K (constant-Cv Ar),
not temperatures of the wet bubble or experimentally achieved states.

[Final acceptance editorial checks](results/final_acceptance_review.md) describe
the current correction separately from historical full-archive verification.
