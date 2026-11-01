# Mathematical, numerical and artifact verification

## New scientific calculations

Eighteen one-cycle integrations cover three amplitudes, two liquid models and
three numerical controls. All complete their requested interval. Strict
DOP853 versus production DOP853/independent Radau gives maximum relative
minimum-radius discrepancy 4.2664e-10 and peak-speed discrepancy 3.6938e-9;
maximum selected-event difference is 1.4668e-16 s. This is numerical evidence,
not physical accuracy, a universal solver bound or a temperature comparison.

The Tait enthalpy integral/derivative, reference density/sound speed, positive
branch, full viscous implicit acceleration residual, zero-drive equilibrium
and RP incompressible limit are checked independently. The common KM event is
compared with the preserved MR2 first-cycle reference. No old thermal curve
is assigned to a new mechanical trajectory.

The saved gas, liquid, Au and summed thermal identities are reassembled from
their independently saved fields/ledgers. Maximum total residual divided by
the thermal scale is 1.3481e-6. The independently integrated nonlinear gas
caloric increment has residual below 4.6073e-8 of absorption. Signed positive
and negative heat are preserved. Original paired full-water-state residuals
are identified as inherited evidence, not a new reconstruction of missing states.

The ideal caloric operation uses the existing variable-capacity integral and
inverts it within its domain. It includes no unmodeled composition, high-T
coefficient extension, ionization fraction or reaction count. It is explicitly
not a dynamic/universal bound. Particle slip/Reynolds and transport times are
defined reference screens, not physical trapping or nucleation models.

Twelve automated tests pass; the complete output is `tests.txt`. Development
assembly correction is recorded separately in `development_note.md`.

## PDFs and standalone reading

Both LaTeX documents compile with no overfull boxes or undefined references.
Final Report has 24 pages (20 pages of introduction through conclusions),
Technical Report 9 pages. Every page was rendered to JPG and visually inspected;
the final one-page contents and corrected scenario cell were re-rendered.
All fonts are embedded. Structural checks confirm no words outside the page,
clickable PDF annotations for MR1–MR7, and absence of commit hashes, Codex,
RoadMap or delivery/acceptance procedures in the scientific PDFs.

`pdf_qa.json` records final hashes, per-page text sizes and link annotations.
`standalone_reader_test.md` maps the required scientific questions to complete
answers inside the Final Report. This is a completeness test, not independent
scientific approval. All seven exact supporting GitHub PDF links return 200
and their raw targets have PDF headers (`github_pdf_links.json`).

## Preservation, scope and reproducibility

All 1491 preexisting tracked files are checked byte-for-byte against the
preflight inventory. No accepted model, result, root plan or historical artifact
is edited. Only MR8 and Final_Report are candidates for publication. Initial
untracked user files remain outside staging. The approximately 1 GB full MR7
archive and its unpublished history are neither required nor uploaded.

Source, pinned environment, numerical inputs, compact data, tables, JPG figures,
LaTeX, numbered results notebook and explicit reproduction commands are provided.
The whole repository is required for read-only transitive model imports.
Public MR7 compact traces suffice. Clean archive verification and manifests
are recorded separately before publication; bitwise regeneration across
platforms is not claimed.

Independent scientific/editorial review remains necessary. No physical maximum,
plasma state, thermonuclear yield or official program closure is certified by
these numerical and artifact checks.
