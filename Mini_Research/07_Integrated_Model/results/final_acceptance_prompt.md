# Final scientific acceptance corrections — local implementation prompt

Work only on the compact MR7 package in `Mini_Research/07_Integrated_Model/`,
starting from published `658908ddad0a2b12a922f6a51ce7a58259a9cccd` after verifying
current remote main and preserving any unrelated edits. The scientific snapshot
is `8547e28fd21042fd92d4907e2a1fc25e38988954`; accepted MR6 is
`493dbdc213f3789b7195d93760bd5c0c867cb98c`.

1. In the generated main PDF table, README and numbered notebook, present
   **Numerical completion** separately from **Physical applicability**.
   The 1.00 bar pair completes five cycles but is a conditional reduced-model
   result with unvalidated AuNP confinement at a prescribed 100 nm gap.
   The 1.15 bar pair completes five cycles but is a **formal fixed-composition
   scenario**, not a validated physical gas state. Explain the prescribed negative
   far-field absolute pressure during part of each period, gas pressure near
   1.30e8 Pa, omitted chemistry/phase change, mechanical-model limits and
   incompressible liquid-transport limits. Negative liquid pressure alone is not
   a mathematical impossibility proof: metastable tension, cavitation and stability
   are separate unresolved physical questions. Wet 1.30 bar remains formal and
   coefficient-censored, with an undefined five-cycle global temperature difference.
2. Immediately beside melting/boiling/phase-transition discussion, state that
   static MR7 peak Au temperature rises are only approximately 21–30% of resolved
   finite-Au MR6 values. Integrated gas-energy delivery differences of about
   0.04–1.68% do not validate local Au temperatures or material admissibility.
   Material thresholds are diagnostic references, and the static discrepancy is
   not a universal error bound for the moving model.
3. Retain the unresolved late 1.15 bar response (89.6% mesh sensitivity), wet
   coefficient stop at 6000 K, and formal dry baselines 6263.579 K (air) and
   22808.716 K (constant-Cv Ar). Do not fill missing tails or call 6000 K a physical
   ceiling. Preserve finite increments +1.4088163 and +2.1848791 microkelvin and
   the undefined wet 1.30 bar global difference.

Edit only report/prose/table/notebook presentation, `src/presentation.py`, directly
related descriptive evidence, inventories and checksums. Preserve every equation,
other source module, parameter, CSV/NPZ, accepted trajectory, MR1–MR6 file,
compact selection and historical verification record byte-for-byte. Provide an
executable comparison against 658908d. Regenerate presentation from saved data,
verify unchanged figure hashes, rebuild and inspect every PDF page, run the quick
tests and bounded replay, and validate the complete SHA256 manifest. Do not repeat
the expensive sweeps or publish the full local archive or its commits.

Create a separate correcting commit and publish by ordinary fast-forward push.
Independently verify GitHub SHA, corrected PDF/table/README/manifest and compact
scope. Send a self-contained handoff with immutable links and preservation proof
to the verified existing **Mini-Research 7 (Web)**; confirm full readback before
recording successful delivery. Request an independent verdict of
`MR7 accepted within declared scope` or `MR7 requires further corrections`.
Do not close MR7 on RoadMap's behalf or start MR8.
