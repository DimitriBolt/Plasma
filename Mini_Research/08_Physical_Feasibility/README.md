# Mini-Research 8: Physical Feasibility and Extreme Temperature Limits

**Dimitri Bolt** · Academic supervisor: **Prof. Dr. Ildar R. Gabitov**.

[Technical report](report/08_Physical_Feasibility.pdf) ·
[Numbered results](08_Physical_Feasibility.ipynb) ·
[Consolidated final scientific report](../../Final_Report/Consolidated_Final_Scientific_Report.pdf).

The new bounded work compares one-cycle KM and Gilmore–Tait mechanics with a
common polytropic gas, audits saved signed thermal energies, inverts the existing
gas caloric law for a defined ideal direct-heating operation, and screens particle
retention and transport times. It does not change previous models or move gas
temperatures between mechanical trajectories.

At 1.30 bar, Gilmore gives a 0.513361 µm minimum versus 0.550809 µm for KM.
These are formal mechanical outputs. For nominal 34.3373 fJ absorption by one
particle, ideal fixed-volume wet-gas increments are 0.03347 and 0.02499 K;
the inherited coupled global changes remain 1.408816 and 2.184879 µK.
The 1.15 bar result is physically formal despite numerical completion. Wet
1.30 bar has no global difference after its 6000 K coefficient boundary.

The gas identity is ΔU = −∫p dV + ∫Q(l→g) dt. The combined thermal identity is
Δ(Ug+El+EAu)+Wby+Qbath=Eabs. Mechanical feedback, signed heat and storage
are kept separate; net gas heat is not a photon-delivery fraction.

**Conclusion:** no physically reliable universal upper temperature or
thermonuclear yield is established. Spatial hot regions, ion/electron states,
reacting fuel density and lifetime remain unresolved. Formal extremes are retained.

[Data schema and assumptions](results/data_schema.md) ·
[Primary-source audit](results/source_audit.md) ·
[Numerical validation](results/validation.md) ·
[Standalone-reader assessment](results/standalone_reader_test.md).

## Reproduction

Keep the whole repository: read-only imports use MR2 and the transitive MR1/MR3/MR4
modules. The public compact MR7 traces are sufficient; its full local archive is
not required or published. Use Python 3.12 and `requirements-lock.txt` in an isolated
environment. LaTeX/latexmk and Poppler are additional system tools.

From this directory:

```bash
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=/tmp/plasma-mr8-mpl
python -m src.run
python -m unittest discover -s tests -v
python -m src.presentation
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/08_Physical_Feasibility.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd ../../Final_Report/Consolidated_Final_Scientific_Report.tex
sha256sum -c SHA256SUMS
```

`src.run` performs 18 bounded mechanical integrations, followed by saved-data
analysis. `src.presentation` only reads results and regenerates JPGs/tables.
Numerical replay is compared with the archived values within the recorded
convergence differences, not by demanding bitwise equality across platforms.
Run the manifest check before regeneration to verify the published snapshot;
regenerated numerical files and PDFs need not have identical hashes.

Independent scientific/editorial review and RoadMap closure remain separate.
The delivery record and review handoff are technical supporting files.
