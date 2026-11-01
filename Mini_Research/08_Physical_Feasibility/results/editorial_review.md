# Editorial revision of the published MR8 package

Scientific reference: `083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2`.
Publication/delivery reference before this revision:
`de146b8e4ee18588d2ebe05b82e20bfc8737e02d`.

This is document maintenance, with no numerical model execution. Scientific
source modules, inputs, saved trajectories, results, figures, tables, notebooks
and historical verification/delivery records remain unchanged. The original
scientific commit remains the reference for full numerical reproduction and
its frozen-tree tests. Those historical guards have not been relaxed to accept
the deliberately revised MR2 README.

## Changes

- Added `report/02_Compressible_Bubble_Dynamics.pdf` in MR2 as an exact copy of
  `report/main.pdf`, retaining the original. Both have SHA-256
  `61486641c6155d8b90c3a247119d8587983225b2e40ce2043e9f6b2bba73debd`.
  The MR2 README, both `support_links.tex` files and the Final Report guide use
  the explicit filename. MR2's current file list and manifest include the copy.
- Condensed Final Report sections 3–7 from 4252 to 3365 English tokens
  (20.86%) with the same regular-expression estimate on the LaTeX section span.
  Equation/caption tokens are included; this is an editorial size estimate.
  The model hierarchy, principal results, qualifications, supporting sections,
  comparison values and new Gilmore–Tait evidence remain available standalone.
- Made completed versus censored/formal extrema explicit in the conclusion.
  Added section 11.3, “Comparison with Experimental Evidence and Validation
  Status,” using the existing primary references. Spectroscopic emission
  temperatures, other authors' acetone calculations, component agreement,
  numerical convergence and missing full-model experimental validation are
  distinguished. No physical maximum, DD yield or impossibility is asserted.
- Translated MR8's README into Russian, preserving its original Bash command
  block, relative-link targets, filenames, numerical values, units and models.
  Both scientific reports remain English. Updated current manifests rather
  than overwriting historical checks.

## Document and preservation checks

`editorial_verify.py` is an artifact-only audit. It reads saved files and Git
objects; it imports no scientific models and performs no integrations. Run
from the repository root, using Python with `pypdf`, plus LaTeX and Poppler:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd Final_Report/Consolidated_Final_Scientific_Report.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd Mini_Research/08_Physical_Feasibility/report/08_Physical_Feasibility.tex
python3 Mini_Research/08_Physical_Feasibility/results/editorial_verify.py
```

The six checks pass. They compare 1562 protected original files byte-for-byte
with the scientific reference, verify the unchanged MR2 original and new copy,
README preservation, scientific statements, section-size reduction, PDF links,
English PDF text, embedded fonts, page bounds and compilation diagnostics.
The saved current report hashes and results are in `editorial_checks.json`.

The Final Report has 24 pages and the Technical Report 9 pages. All 33 rebuilt
pages were rendered to JPG and visually inspected, including the contents,
footnotes, tables, guide and the new experimental/conclusion pages. The Final
Report's 52 internal links and 138 named destinations and the Technical
Report's 17 internal links and 64 named destinations resolve. Neither build
has overfull boxes or undefined references. MR2's copy requires no typesetting
change because it is byte-identical to the retained 13-page original.

Current checksum verification is performed separately in MR2, MR8 and
`Final_Report` using `sha256sum -c SHA256SUMS`. A clean staged-tree export is
also checked and both documents are rebuilt there, without rerunning numerical
experiments. Regenerated PDF bytes can vary with build metadata; model/data
preservation is established against Git objects, not by claiming platform-wide
bitwise PDF regeneration.

Publication is a separate delivery step. Remote main and raw PDF/README bytes
must be read back after pushing, including the seven supporting PDF targets
and both MR2 filenames. This record does not replace the preserved scientific
review or prior delivery records. Scientific/editorial acceptance remains with
Mini-Research 8 (Web); official program closure remains with RoadMap. MR9 is
not initiated.
