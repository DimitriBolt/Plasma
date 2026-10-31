# Mini-Research 1: Single-bubble dynamics under acoustic excitation

Author: Dimitri Bolt. Status: **prepared for independent review; stage not closed**.

This study integrates the spherical, incompressible Rayleigh-Plesset equation in water at 20 °C with a constant polytropic gas exponent. It contains original five-cycle simulations at 0, 0.2, 0.5, 0.8, 1.0, and 1.3 bar. The last regime is a **formal numerical stress test outside the physical applicability of the model**. No bubble-temperature, laser, nanoparticle, chemical, or nuclear-reaction model is implemented.

## Files

- `src/bubble_dynamics.py`: dimensional SI model, positivity-aware DOP853, event localization, and independent analytic linear response.
- `src/run_experiments.py`: amplitude sweep, convergence including Radau, equilibrium and linear checks, IAPWS reference-property audit, source hashes and metadata.
- `src/plot_results.py`: nine reproducible illustrations (vector PDF and 300 dpi PNG) and LaTeX tables, using saved results only.
- `tests/test_bubble_dynamics.py`: analytic checks, invalid inputs, floor termination, linear limit, and second-method agreement.
- `results/trajectories.csv`: time, radius, velocity, far-field pressure and wall Mach number, all in SI units.
- `results/amplitude_summary.csv`: extrema, local-event counts, deepest-collapse time and duration, completion status and physical warnings.
- `results/convergence.csv`: four DOP853 accuracy settings and an independent Radau calculation at 1.3 bar.
- `results/extrema.csv`: refined radius and velocity extrema; excludes the initial state.
- `results/linear_comparison.csv`: independent linear-response comparison at 100 and 500 Pa.
- `results/run_metadata.json`: solver settings, installed versions, equilibrium and linear errors, IAPWS audit, and SHA256 provenance.
- `results/report_*.tex`: tables and numeric macros generated from saved results.
- `report/mini_research_01.tex` and `.pdf`: English scientific report and supplementary figures.
- `requirements.txt`: checked direct dependencies; `requirements-lock.txt`: complete tested environment.

## Reproduce

From the repository root, with Python 3.12 (tested: 3.12.14), a network connection for installation, and a local LaTeX installation with `latexmk`, `pdflatex`, `amsmath`, `booktabs`, `graphicx`, `geometry`, `hyperref`, `caption`, and `float`:

```bash
cd Mini_Research/01_Bubble_Dynamics
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock.txt
python -m src.run_experiments
python -m src.plot_results
python -m pytest -q
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/mini_research_01.tex
```

`requirements.txt` can be used instead of the full lock for the direct dependency versions. Scripts resolve paths relative to their own location, independently of the current working directory. By default, the simulation checks the existing repository inputs `Library/035218_1_online.pdf`, `Library/nature03361.pdf`, and `Plasma-plan.tex`. For another arrangement pass `--source-root PATH_TO_REPOSITORY` to `src.run_experiments`. The source PDFs and plan are read and hashed, never modified. To remove only build intermediates use `latexmk -c -cd report/mini_research_01.tex`.

## Interpretation and safeguards

The equation uses dynamic viscosity in Pa s, not kinematic viscosity. Gas partial pressure excludes the separately added constant vapor pressure. The acoustic sign is `p_infty = p0 - Pa*sin(2*pi*f*t)`, so the first quarter-cycle lowers pressure. All sweep parameters other than amplitude remain fixed.

Production settings: five periods; DOP853 with `rtol=1e-10`, componentwise `atol=(1e-16 m, 1e-8 m/s)` and maximum step `T/1000`. Invalid internal trial radii are rejected and retried at one quarter of the step by a small DOP853 subclass; no radius clipping or change to the ODE is used. A terminal floor at `0.001 R0 = 6 nm` is technical, not a physical radius prediction. Reaching it censors the result. A failed or incomplete sweep raises an error.

Events detect zeros of velocity and acceleration. Dense-output roots are refined with `brentq`, absolute time tolerance `1e-24 s`, because SciPy's default absolute event-time tolerance is insufficient for the unphysical, sub-picosecond final rebound of the 1.3 bar stress test. This tolerance is a root-search setting, not a claim of absolute physical timing accuracy. A common time origin remains in seconds; floating-point and ODE errors are assessed by convergence. Exports combine a uniform grid, adaptive-step subdivisions, and refined event times; plots do not re-integrate or infer extrema from a coarse plotting grid.

Global extrema include initial and final states. Local counts omit `t=0`. The characteristic collapse is the deepest post-initial local minimum within five cycles; its duration is measured from the preceding local maximum. This definition is not the empty-cavity Rayleigh collapse time. Rebounds produce several local extrema per acoustic cycle. Five cycles do not establish asymptotic periodicity or stability.

Wall Mach numbers of order 0.1 already warrant caution; values near or above one invalidate incompressible collapse predictions. Negative far-field absolute pressure is flagged where `Pa > p0`, rather than clipped. In particular, the 1.3 bar extrema are numerical properties of this idealized ODE, not credible physical collapse predictions. Agreement between integrators verifies the calculation only. Heat/mass transfer, varying vapor pressure, gas chemistry, shape modes, finite gas volume, shocks and liquid compressibility remain outside the scope. These calculations do not replicate the 2016 thermodynamic simulations or quantitatively validate against the 2005 sulfuric-acid experiment.

## Sources

- Sadighi-Bonabi and Gheshlaghi (2016), Section II, equations (1)-(4), Tables I-IV: https://doi.org/10.1063/1.4945343
- Flannigan and Suslick (2005), experimental motivation in concentrated aqueous sulfuric acid: https://doi.org/10.1038/nature03361
- IAPWS-95: https://iapws.org/technical-guidance/release/IAPWS-95
- IAPWS surface tension: https://iapws.org/technical-guidance/release/Surf-H2O
- IAPWS viscosity: https://iapws.org/technical-guidance/release/viscosity
- Open-source reference implementation (independent of the bubble model): https://github.com/jjgomera/iapws
- SciPy integration documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

The rounded model properties and the independently evaluated IAPWS values are both recorded. `Plasma-plan.tex` supplies broader context; this deliverable covers only the prescribed radius-dynamics subset. Stage closure requires separate mathematical, reproducibility, and PDF review.
