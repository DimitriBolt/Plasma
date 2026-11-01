# Mini-Research 6: Frozen-Geometry AuNP-Bubble Thermal Coupling

**Dimitri Bolt**. Academic supervisor: **Prof. Dr. Ildar R. Gabitov**.

[Scientific report](report/06_AuNP_Bubble_Coupling.pdf) ·
[numbered saved-results notebook](06_AuNP_Bubble_Coupling.ipynb) ·
[model and API](results/model_details.md) · [data schema](results/data_schema.md).

The finite-particle axisymmetric heat equation resolves local gas-wall delivery:

\[
C_{Au}\dot\theta_{Au}=P_{abs}-\int G(\theta_{Au}-\theta_l)dA,\qquad
\delta\dot Q_{lg}=\int_{S_b}h\theta_l\,dA.
\]

The computed geometry is frozen. The signed gas closure is **MR4 baseline plus
only the laser increment**. Gold/water storage, gas delivery and bath escape
are separate energy channels. The [reference](results/reference.npz) uses a
20 nm particle, 100 nm surface gap, 6 micrometre bubble, 532 nm, local 1 J/m²
and a 10 ns pulse. About **0.326%** of its 3.434 fJ absorption reaches the gas
within 2 microseconds in this static model. This is conditional on the effective
gas conductance and frozen liquid. Event-specific flow/displacement checks fail
near strong collapse, and passive particles move away from the wall.

The [finite design](results/finite_design.csv) separates maximum absorption,
gas energy at 2 microseconds and efficiency. Its moderate static gas-energy
maximum is about **3.95 fJ** from about **3.48 pJ** local absorption.
Wavelength plateaus, finite grid boundaries and numerical uncertainty remain
explicit. [Formal regimes](results/formal_regimes.csv) preserve near-contact,
short-pulse and large-temperature failures. None is an achieved temperature.

## Reproduction

Python **3.12.14**, full [dependency lock](requirements-lock.txt). The accepted
sibling stages and project plan are read only. From this directory:

```bash
sha256sum -c SHA256SUMS
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m pip check
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/plasma-mr6-mpl
python -m src.experiments
python -m src.validate
python -m src.presentation
python -m pytest -q
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/06_AuNP_Bubble_Coupling.tex
```

The final scope is stationary AuNP-to-liquid-to-bubble coupling. The 0.3258%
reference efficiency belongs only to the specified geometry and 2 microsecond
window. Moving boundaries, advection and self-consistent mechanics are reserved
for MR7. The retained event/A/B/C results are diagnostic context.

For the **bounded final-correction verification**, use the locked environment and
saved data, without repeating the original sweeps:

```bash
python -m src.final_corrections --replay
python -m src.presentation
python -m pytest -q
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/06_AuNP_Bubble_Coupling.tex
```

This checks 59 numerical handoff claims, independently integrates energy/flux
from the five archived spatial fields, evaluates the saved-trace analytic gas
feedback estimate, and replays only three old thermal cases on saved meshes.
A platform's remesher may choose different nodes; exact regenerated-mesh equality
is not a test requirement. The 2e-5 replay tolerance is separate from the stored
spatial-convergence evidence. `src.clean_verify` remains an optional historical
same-machine exact full-recalculation comparator, not a cross-platform criterion.

The original full design, trajectories and convergence data remain preserved.
No new optimization or moving-boundary calculation is part of the correction.
The corrected figure and scientific PDF are rebuilt from saved results. Eight
figures are JPG only. `clean_reproduction_original.json` retains historical full
recalculation evidence; the [current clean record](results/clean_reproduction.json)
and [verification journal](results/validation.md) describe the bounded fresh run.
Manifest checks apply to published bytes before regeneration.

## Source interface and scope

[coupling.py](src/coupling.py) supplies a stateful `FrozenHeatModel` and a
`ThermalHistory` with `thermal_source(t_s, history, bubble_state=None)`.
The history explicitly labels the stored reference-state baseline/sum and the
laser increment; no `q_wall_total_W` is exported. A gas ODE already computing its
current-state heat loss adds **only** `history.laser_increment_W(t_s)`. The
reference sum is not a physical total for another Tg. Mechanism A is a transient
early-heating diagnostic without a new equilibrium radius, vapor pressure or
composition. The precomputed response stores its thermal memory and finite validity interval;
it rejects extrapolation and inconsistent supplied radius. Histories in
`results/histories/` separate signed baseline and laser power and record absolute
time, collapse time, units, geometry and applicability. They are **formal frozen
event diagnostics**, not sources already validated for changed MR7 mechanics.

The original work begins at accepted MR5 snapshot
`0204ce7a8429406b7ec235b02cec8581274fc8c6` and fingerprints all 1083 tracked
upstream inputs. No MR1-MR5 model/data/report, general plan or unrelated user
file is changed. This separate correction follows the original MR6 snapshot
`faa2828f554517a0e2e43d9ae42cff37d682fbd0`. It is submitted for independent review; publication
does not close the stage, update RoadMap or begin MR7.

[Independent review handoff](results/review_handoff.md) contains equations,
numerical results, limitations, checks and unresolved questions. Delivery status
is stated there and in the final publication message, without assuming a chat
exists or reusing a preceding stage's conversation identity.
