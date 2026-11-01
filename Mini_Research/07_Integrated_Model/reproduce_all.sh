#!/usr/bin/env bash
# Explicit expensive numerical replay. Run in a separate checkout, sequentially.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
export MR7_DATA_MODE=full
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/mr7-mpl}"
mkdir -p results/raw
python -m src.benchmarks
python -m src.analytic_control
python -m src.static_domain_control
python -m src.experiments main
python -m src.experiments sensitivity
python -m src.local_controls
python -m src.experiments convergence
python -m src.extra_controls
python -m src.regression
python -m src.independent_radau
python -m src.high_sensitivity
python -m src.high_mesh
python -m src.stable_pairs --names base_1.00 base_1.15 extra_high_1.30_BDF_base
python -m src.scientific_audit
python -m src.physical_scales
python -m src.transport_view
python -m src.numerical_scope
python -m src.compact
python -m unittest discover -s tests -v
python -m src.presentation
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd report/07_Integrated_Model.tex
# Inspect the new PDF before deliberately refreshing its publication manifest.
