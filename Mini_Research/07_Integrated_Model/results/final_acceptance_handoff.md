delivered_and_verified

# MR7 final acceptance corrections — independent re-audit requested

The correcting scientific/documentation commit is published and independently verified. The complete submitted handoff was read back from the verified existing Mini-Research 7 (Web), with an exact match after stripping outer whitespace (2026-11-01T05:21:59.470416+00:00; 8,526 received characters). The submitted-payload hash is recorded in final_acceptance_delivery.json. The delivery status was updated after this readback. A repeated independent binary verdict is requested; no RoadMap closure is asserted and MR8 has not been started.

Delivery marker: MR7-FINAL-ACCEPTANCE-b2c37a2467e863bbc444278c06b69b060a70762d

Correcting SHA and remote main at correction publication: `b2c37a2467e863bbc444278c06b69b060a70762d`.
Parent/baseline: `658908ddad0a2b12a922f6a51ce7a58259a9cccd`. Scientific input remains
`8547e28fd21042fd92d4907e2a1fc25e38988954`.

## Immutable publication and exact diff

- [Correcting commit](https://github.com/DimitriBolt/Plasma/commit/b2c37a2467e863bbc444278c06b69b060a70762d)
- [Corrected PDF](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/report/07_Integrated_Model.pdf)
- [Exact complete comparison against the baseline](https://github.com/DimitriBolt/Plasma/compare/658908ddad0a2b12a922f6a51ce7a58259a9cccd...b2c37a2467e863bbc444278c06b69b060a70762d)
- [Generated main table](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/table_main.tex)
- [README](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/README.md) and [numbered notebook](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/07_Integrated_Model.ipynb)
- [Executable preservation check](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/final_acceptance_check.py)
- [Per-file preservation proof](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/final_acceptance_preservation.json)
- [Editorial review and locations](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/final_acceptance_review.md)
- [Quick tests](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/final_acceptance_tests.txt), [bounded replay](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/final_acceptance_replay.json), [LaTeX log](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/final_acceptance_latex.txt)
- [File inventory](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/FILE_LIST.txt) and [SHA256 manifest](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/SHA256SUMS)
- [Unchanged compact selection](https://github.com/DimitriBolt/Plasma/blob/b2c37a2467e863bbc444278c06b69b060a70762d/Mini_Research/07_Integrated_Model/results/publication_scope.json)

## Three corrections

1. PDF Table 2 and adjacent text on page 7, the abstract, conclusions (page 12),
   README and notebook now separate **Numerical completion** from **Physical
   applicability**. The 1.00 bar five-cycle result is conditional on the reduced
   model and unvalidated prescribed AuNP confinement at a 100 nm gap. The 1.15 bar
   pair completes five numerical cycles but is a **formal fixed-composition
   scenario**, not a validated physical gas state. The accompanying text gives
   prescribed negative far-field absolute pressure during part of a period,
   gas pressure near 1.30e8 Pa, omitted chemistry/phase changes and limitations
   of mechanics and incompressible liquid transport. Negative liquid pressure
   alone is explicitly not a mathematical impossibility proof; metastable tension,
   cavitation and stability remain separate unresolved physical questions.
2. The warning immediately beside material thresholds on PDF page 11 states that
   static MR7 **peak Au temperature rises** are only 21–30% of resolved finite-Au
   MR6 values. Integrated gas-energy differences of 0.04–1.68% do not validate local
   Au temperatures, melting, boiling or nearby phase-transition admissibility.
   Material thresholds are diagnostic references. The static discrepancy is not
   a universal error bound for the moving model. README/notebook are consistent.
3. Page 9 explicitly retains **89.6%** late 1.15 bar mesh sensitivity and
   **unresolved** status. Wet 1.30 bar remains coefficient-censored at 6000 K,
   with undefined five-cycle global difference. Dry 6263.579 K air and 22808.716 K
   constant-Cv Ar baselines remain formal. No tail is filled, no wet trajectory is
   substituted, and 6000 K is not a physical temperature ceiling.

## Numerical preservation and current verification

The finite nonlinear differences remain +1.4088163 and +2.1848791 microkelvin;
the wet 1.30 bar global difference remains undefined. No equation, source/flow
model, solver, parameter, CSV, NPZ or accepted trajectory has changed.
The executable comparison against `658908ddad0a2b12a922f6a51ce7a58259a9cccd` confirms:

- 1468 protected preexisting files byte-identical, including all 1209 upstream
  files, 28 non-presentation source modules and all 45 CSV/NPZ files;
- all 14 equation/align blocks and all 10 JPG figures unchanged;
- compact selection, metadata and historical full-archive checks unchanged.

All 34 existing quick tests pass. The requested bounded one-cycle replay passes,
with maximum gas-temperature difference 1.5960e-5 K and radius difference
3.0183e-13 m on its common prefix. Full sweeps were not repeated. All 277 manifest
entries pass. GitHub independently confirms all 278 MR7 file blob identities at
the correcting SHA. The package is 45,833,286 bytes; the full local archive and
its unpublished commits are not included.

The rebuilt PDF has 14 pages (12 main, two appendix), all visually inspected.
Final LaTeX has no overfull boxes or undefined references; fonts are embedded.
Historical 13-page checks remain frozen as records of the previous snapshot.
The new final_acceptance_check.py applies to the corrected PDF and does not
relabel historical full-archive verification as current compact verification.

## Changed files in the correcting commit

- `Mini_Research/07_Integrated_Model/07_Integrated_Model.ipynb`
- `Mini_Research/07_Integrated_Model/FILE_LIST.txt`
- `Mini_Research/07_Integrated_Model/README.md`
- `Mini_Research/07_Integrated_Model/SHA256SUMS`
- `Mini_Research/07_Integrated_Model/report/07_Integrated_Model.pdf`
- `Mini_Research/07_Integrated_Model/report/07_Integrated_Model.tex`
- `Mini_Research/07_Integrated_Model/results/data_schema.md`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_check.py`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_latex.txt`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_preservation.json`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_prompt.md`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_replay.json`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_review.md`
- `Mini_Research/07_Integrated_Model/results/final_acceptance_tests.txt`
- `Mini_Research/07_Integrated_Model/results/model_decisions.md`
- `Mini_Research/07_Integrated_Model/results/narrative_conclusion.tex`
- `Mini_Research/07_Integrated_Model/results/narrative_primary.tex`
- `Mini_Research/07_Integrated_Model/results/table_main.tex`
- `Mini_Research/07_Integrated_Model/src/presentation.py`

Only ten existing editorial files, the two inventories and new final_acceptance
verification/documentation records changed. All public reference arrays are the
same selected arrays as before; other high-dimensional sweep outputs remain local.

## Requested independent verdict

Please independently check the three corrections, numerical/physical separation,
preservation proof, exact diff and corrected PDF. Return one of:

- `MR7 accepted within declared scope`;
- `MR7 requires further corrections`.

This request is a repeated independent review, not a declaration of RoadMap closure.
Do not begin MR8. The local Codex editorial task ends after verified delivery;
RoadMap retains the official stage decision.

The delivery-receipt commit changes only new handoff/verification metadata and inventories. Scientific links above pin the correcting commit b2c37a2467e863bbc444278c06b69b060a70762d. Historical receipts, numerical results and computational code remain unchanged.
