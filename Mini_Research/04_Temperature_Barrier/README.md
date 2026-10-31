# Mini-Research 4: Temperature-Barrier Physics

**Dimitri Bolt**. Academic supervisor: **Prof. Dr. Ildar R. Gabitov**.

The [English report](report/04_Temperature_Barrier.pdf) and
[LaTeX source](report/04_Temperature_Barrier.tex) compare variable heat capacity,
fixed water-vapor inventory and mechanical feedback with Mini-Research 3.
The central result is [the MR3/MR4 temperature comparison](results/main_comparison.csv).
These are conditional model estimates prepared for independent scientific review.

For a closed frozen mixture,

\[
C_V^{tot}(T)\dot T=-p_g\dot V-4\pi\eta k_gR(T-T_\ell),\qquad
p_g=n_{tot}R_uT/V.
\]

The pressure derivative follows this same energy equation. NASA9 molecular
coefficients and their ranges are archived; Ar uses the constant ground-state
monatomic capacity. A coefficient boundary censors a trajectory, not the physical
temperature. One-way pressure mismatch and coupled KM validity diagnostics remain
available, including every 1.3 bar scenario. Chemistry is a quantitative **energy
assessment**, with no computed reacting composition or chemically corrected Tmax.
There is no mass exchange, latent heat, laser heating, nanoparticle physics or MR5.

## Controlled comparisons

| Saved model | Inventory / pressure | Heat capacity |
|---|---|---|
| M0 | MR3 dry inventory plus separate constant pv | Exact MR3 constants |
| M1_const / M1_variable | Same dry components and legacy pv | Frozen at T0 / NASA variable |
| M2_dry | Dry closed mixture; no separate pv; equilibrium total pressure | Variable |
| M2_sat_const / M2_sat | Saturated initial H2O inventory; no separate pv | Frozen at T0 / variable |
| M2_v5 / M2_v10 | Conditional trapped H2O, 5 or 10 times saturation fraction | Variable |

For vapor sensitivity, total initial pressure and total moles are fixed; dry gas
is displaced by H2O. Enhanced fractions are supersaturated hypothetical trapped
states, not phase equilibria at 293.15 K. All heat-transfer comparisons retain the
same gas-specific MR3 conductivity. The heat-transfer sensitivity uses eta=0.5,1,2,
including air at 1.15 bar to provide complete-interval comparisons alongside
censored extreme scenarios.
One-way paths are the same saved MR2 mechanics reconstructed in MR3, not newly
coupled MR4 paths. Fully coupled results solve radius and thermodynamics together.

## Reproduce

The read-only sibling MR1--MR3 folders, `../../Plasma-plan.tex`, and all three
named PDFs in `../../Library/` must be present. The source inventory checks their
bytes. The published start point is recorded in `results/input_inventory.json`.
Numerical rebuilding needs no network once Python dependencies are installed.
No old mechanics or accepted scientific results are regenerated or modified.

From this folder, using Python 3.12:

```bash
sha256sum -c SHA256SUMS
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip check
export MPLCONFIGDIR=/tmp/plasma-mr4-matplotlib
python -m src.run_experiments --force
python -m src.validate --workers 4 --force
python -m src.chemistry
python -m src.presentation
python -m unittest discover -s tests -v
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/04_Temperature_Barrier.tex
python -m src.package
```

Without `--force`, saved runs and validation cases are reused. After modifying
equations or settings, use `--force`. The validation worker count controls only
independent numerical calculations. LaTeX needs standard amsmath, amssymb,
booktabs, graphicx, geometry, lmodern, caption, fancyhdr, float, microtype and
hyperref; Poppler supports visual verification. Figures are JPG only, 300 dpi,
with English labels and white backgrounds. Generated metadata/PDF binary bytes
can change on rebuilding; numerical reproducibility is assessed with the saved
convergence tolerances, separately from checksums of the published snapshot.

## Evidence and data

- [Model / source audit](results/source_audit.md), [input hashes](results/input_inventory.json),
  [coefficient records and intervals](results/species_coefficients.json),
  [parameters](results/parameters.json), [schema](results/data_schema.md).
- [All runs](results/summary.csv), [MR3 regression](results/mr3_regression.csv),
  [paired sensitivities](results/sensitivity.csv), [boundary events](results/boundary_events.csv).
- [Energy budgets](results/energy_budgets.csv), [independent quadrature](results/independent_quadrature.csv),
  [solver and input convergence](results/convergence.csv), [validation](results/validation.md).
- [Peak-time conditioning](results/peak_time_conditioning.csv) distinguishes
  numerical time differences from selection among nearly tied periodic peaks.
- [Clean-copy verification](results/clean_snapshot_verification.json) and
  [fresh numerical controls](results/clean_reproduction.csv) document a separate
  source/data/PDF copy with the same locked dependencies. `python -m src.clean_verify`
  repeats the four bounded numerical controls from the current checkout.
- [Reaction energy scales](results/reaction_energy_scales.csv),
  [inventory energy demands](results/chemistry_inventory_budgets.csv),
  [chemical inputs and primary links](results/chemical_parameter_sources.json).
- [Independent-review handoff](results/review_handoff.md),
  [publication file list](FILE_LIST.txt), [SHA256 manifest](SHA256SUMS).

`results/trajectories/` contains each summary JSON, a saved physical series NPZ
and accepted integrator states. `results/validation_cases/` retains each stricter,
Radau or input-resolution check. The report tables/figures use these saved results.
Numerical convergence does not establish ideal-gas, KM, thermal-equilibrium or
experimental validity; neither table limits nor the finite design supply a
universal temperature ceiling.
