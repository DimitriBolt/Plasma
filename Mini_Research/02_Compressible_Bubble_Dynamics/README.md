# Mini-Research 2: compressible bubble dynamics

Five-cycle Keller-Miksis (KM) mechanics with the exact gas law and fixed water constants of accepted Mini-Research 1 (RP). Author: Dimitri Bolt. Academic advisor: Prof. Dr. Ildar R. Gabitov. No temperature, reactions, optical heating or nanoparticles are calculated. The general project plan is unchanged.

The scientific report is [report/main.pdf](report/main.pdf); its English LaTeX source is [report/main.tex](report/main.tex). Numerical results and figures are generated from saved data, including formal extreme trajectories and explicitly censored incomplete RP solutions. See [results/validation.md](results/validation.md), [results/source_audit.md](results/source_audit.md), and [results/data_schema.md](results/data_schema.md).

## Reproduce from a clone

Python 3.12 with the locked dependencies is the clean verification environment. Production sweep metadata records the Python 3.14 environment used for its initial generation. A clean Python 3.12 rerun is compared independently. System LaTeX requires latexmk/pdflatex and lmodern, amsmath, amssymb, booktabs, graphicx, geometry, caption, float, microtype and hyperref. Poppler is used for visual QA only.

From the repository root:

```bash
cd Mini_Research/02_Compressible_Bubble_Dynamics
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip check
python -m src.experiments --force
python -m src.extended_validation
python -m src.presentation
python -m src.audit
python -m pytest -q
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/main.tex
```

`--force` repeats ALL mechanical runs. Without it, `src.experiments` reuses per-amplitude saved trajectories but reruns convergence and the c-limit studies. Stages are selectable with `--stage sweep`, `--stage convergence`, `--stage c-limit`, `--stage metadata`. To rebuild only figures/tables/report from saved data, run the last four commands (audit/tests/LaTeX do not rerun the sweep, except the deliberately bounded reproduction/equilibrium tests).

The unchanged sibling `../01_Bubble_Dynamics/src/bubble_dynamics.py` and its accepted CSVs are required. The three named local papers in `../../Library/` and `../../Plasma-plan.tex` are checked and hashed. Remote literature is documented with exact readable full-text links in the source audit, so a numerical rebuild itself needs no network after installation.

## Conventions and interpretation

The main equation uses `(1+U/c)*Delta_p/rho + R/(rho*c)*dDelta_p/dt`, including the full surface and viscous derivatives. Its explicit denominator is `D=R*(1-U/c)+4*mu/(rho*c)`. A zero is possible on the outward branch, not simply when an inward wall reaches abs(U)=c. No pressure, speed or denominator is clipped and no Mach threshold removes a trajectory.

Global extrema include endpoints. The characteristic collapse remains the deepest local minimum and its preceding maximum, as in MR1. Separate per-cycle pairs and the first minimum below R0/2 avoid confusing different cycles. Incomplete solutions have no five-cycle global extrema; plotting inserts gaps and comparison differences are suppressed. Adaptive event-aware exports preserve the narrow formal RP core as well as the wider KM collapse.

Numerical convergence does not establish physical validity. Mach and pressure/bulk-modulus scales are continuous diagnostics. Shape stability, gas uniformity, fixed kappa, mass transfer, excluded volume and continuum/EOS validity remain explicitly unresolved. Negative absolute far-field water pressure can describe metastable tension and is not automatically impossible. A separate, validated Gilmore comparison is specified for 1.0, 1.15 and 1.3 bar, with no unfinished Gilmore code included.

## Package and source integrity

`src/`: model/driver, experiments, saved-data presentation, extended numerical checks, audit and manifest tools. `tests/`: equation, domain, events, numerical and saved-data tests. `results/`: all input parameters, trajectories, extrema, summaries, convergence, source audit and admission records. `figures/`: 13 English figures, each in PDF and PNG. `report/`: English source and compiled PDF. `FILE_LIST.txt` and `SHA256SUMS` list this package, excluding caches and LaTeX auxiliaries.

The only changes in Mini-Research 1 are its README status/advisor and the advisor name in validation.md. Its original SHA256SUMS remains a historical manifest of commit d59bc0d; those two authorized documentation edits naturally change their hashes. The accepted model and numerical artifacts are byte-for-byte preserved. User files outside this scope are not added or changed.
