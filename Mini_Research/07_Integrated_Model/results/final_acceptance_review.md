# Final acceptance editorial review

Baseline: `658908ddad0a2b12a922f6a51ce7a58259a9cccd`.
Scientific input snapshot: `8547e28fd21042fd92d4907e2a1fc25e38988954`.

The three requested presentation corrections are satisfied. This is the local
editorial self-check; a repeated independent Web verdict is requested separately.
It does not close MR7 on RoadMap's behalf or authorize MR8.

## Corrections and locations

| Correction | Corrected locations | Preserved interpretation |
|---|---|---|
| Numerical completion versus physical applicability | PDF abstract, Table 2 and adjacent text on page 7, conclusions on page 12; generated table and narratives; README status table; notebook sections 2 and 6 | 1.00 bar: completed but conditional with unverified AuNP confinement. 1.15 bar: completed but formal fixed-composition. Wet 1.30 bar: coefficient-censored with undefined global difference. |
| Local Au heating and material thresholds | PDF page 11 immediately beside melting, boiling and phase-transition discussion; README; notebook sections 5 and 6; model-decision clarification | Peak Au rises are only 21–30% of resolved static MR6 values. Integrated gas-energy differences of 0.04–1.68% do not validate local material admissibility. No universal moving-model error bound is inferred. |
| Unresolved and extreme states | PDF pages 9 and 11; README; notebook sections 4 and 6; unchanged numerical tables and arrays | 89.6% late 1.15 bar mesh sensitivity remains unresolved. Wet 6000 K stop is a coefficient boundary, not a physical ceiling. Dry 6263.579/22808.716 K baselines remain formal. No tail is filled. |

The 1.15 bar explanation explicitly includes prescribed negative far-field absolute
pressure during part of a cycle, gas pressure near 1.30e8 Pa, fixed-composition and
phase-change limitations, first-order mechanics and incompressible transport.
Negative liquid pressure alone is not called a mathematical impossibility proof;
metastable liquid tension, cavitation and stability remain separate physical questions.

## Machine-verifiable preservation

Run from a Git checkout containing the baseline:

```bash
python results/final_acceptance_check.py --visual-reviewed --verify-manifest
```

The visual flag records the completed manual inspection; it is not an automated
substitute for viewing pages. `final_acceptance_preservation.json` lists each
protected file's baseline/current Git blob identity and SHA256.

- 1468 protected preexisting files are byte-identical, including all 1209 upstream
  files, 28 non-presentation source modules, all 45 CSV/NPZ files, parameters,
  metadata, compact selection and historical verification records.
- All 14 displayed equation/align blocks and all 10 original JPG figures are unchanged.
- Finite increments remain +1.4088163 and +2.1848791 microkelvin; the wet 1.30 bar
  five-cycle difference remains undefined.
- Existing full-archive validation records are neither edited nor relabeled as
  fresh compact checks. The earlier 13-page PDF records describe the old snapshot;
  the new check above applies to this corrected 14-page PDF.

## Current checks

The corrected PDF has 12 main-text pages and two appendix pages. All 14 pages
were rendered to JPG and inspected; only pages 7 and 12 changed during the final
wording pass and were rechecked. The final LaTeX log has no undefined references
or overfull boxes. Fonts are embedded and tables are not clipped.

All 34 existing quick tests pass (`final_acceptance_tests.txt`). A bounded fresh
one-cycle replay passes (`final_acceptance_replay.json`), with maximum gas-temperature
difference 1.5960e-5 K and radius difference 3.0183e-13 m on the declared common
prefix. Full parameter sweeps were not repeated. The first preservation-check
attempt required NUL-separated Git paths for Unicode filenames; no scientific
file was changed to resolve that path-encoding issue.

The generated table/notebook and prose remain reproducible with
`python -m src.presentation`. The correction changes only the ten editorial files
listed by the preservation record, inventories, and new `results/final_acceptance*`
documentation/check records. `FILE_LIST.txt` supplies the complete publication list.
The compact selection remains unchanged; no full local archive or its commits are added.
