import numpy as np
import pytest
from src.bubble_dynamics import (Parameters, external_pressure, integrate,
                                 linear_solution, rhs, summarize)


@pytest.fixture
def p():
    return Parameters()


def test_equilibrium_balance_and_trajectory(p):
    assert p.pg0 > 0
    assert abs(p.pg0+p.pv-p.p0-2*p.sigma/p.R0) < 1e-9
    sol = integrate(0,p)
    assert sol.success and sol.status == 0
    y = sol.sol(np.linspace(0,5*p.T,10001))
    assert np.max(abs(y[0]/p.R0-1)) < 1e-8
    assert abs(y[1]).max() < 1e-6
    assert summarize(sol,0,p)["local_minima_count"] == 0


def test_forcing_phase_and_bar_conversion(p):
    Pa = 1e5
    np.testing.assert_allclose(external_pressure(np.array([0,p.T/4,3*p.T/4]),Pa,p),
                               [p.p0,p.p0-Pa,p.p0+Pa],atol=1e-9)
    assert external_pressure(p.T/4,1.3e5,p) < 0


def test_acceleration_matches_pressure_force(p):
    # At rest and R0 all static terms cancel, leaving acoustic forcing / rho R0.
    np.testing.assert_allclose(rhs(p.T/4,[p.R0,0.],500.,p),
                               [0.,500./(p.rho*p.R0)],rtol=1e-12)


@pytest.mark.parametrize("radius", [0.,-1e-6])
def test_invalid_radius_rejected(p,radius):
    with pytest.raises(ValueError):
        rhs(0,[radius,0],1e5,p)


@pytest.mark.parametrize("amplitude", [-1.,np.nan])
def test_invalid_amplitude_rejected(p,amplitude):
    with pytest.raises(ValueError):
        integrate(amplitude,p)


def test_floor_terminates_without_false_complete():
    p = Parameters(floor_fraction=.9)
    sol = integrate(1e5,p)
    assert sol.status == 1 and len(sol.t_events[0]) == 1
    assert sol.y[0,-1] == pytest.approx(p.R_floor,rel=1e-8)
    summary = summarize(sol,1e5,p)
    assert not summary["completed"]
    assert "metrics_censored" in summary["warnings"]


def test_linear_limit_is_quadratic(p):
    t = np.linspace(0,5*p.T,10001)
    errors = []
    for Pa in (100.,500.):
        nonlinear = integrate(Pa,p).sol(t)[0]
        analytic = linear_solution(t,Pa,p)[0]
        errors.append(np.max(abs(nonlinear-analytic))/p.R0)
    assert errors[0] < 2e-6
    assert 20 < errors[1]/errors[0] < 30


@pytest.mark.parametrize("Pa", [1.0e5, 1.3e5])
def test_extrema_resolved_and_method_agreement(p, Pa):
    a = summarize(integrate(Pa,p,steps_per_period=500),Pa,p)
    b = summarize(integrate(Pa,p,method="Radau",steps_per_period=500),Pa,p)
    assert a["completed"] and b["completed"]
    assert a["local_minima_count"] > 5
    assert a["collapse_time_s"] > 0
    assert a["Rmin_m"] > p.R_floor
    for key in ("Rmax_m","Rmin_m","max_abs_U_m_s","collapse_time_s"):
        assert a[key] == pytest.approx(b[key],rel=2e-5)


def test_saved_sweep_matches_trajectory_extrema():
    """Catch stale summaries or a coarse export that loses an extremum."""
    import csv
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    data = np.genfromtxt(root / "results/trajectories.csv", delimiter=",", names=True)
    with (root / "results/amplitude_summary.csv").open() as f:
        summary = list(csv.DictReader(f))
    assert {float(r["Pa_Pa"]) for r in summary} == {0, 20000, 50000, 80000, 100000, 130000}
    for row in summary:
        d = data[data["Pa_Pa"] == float(row["Pa_Pa"])]
        assert row["completed"] == "True" and row["solver_status"] == "0"
        assert d["t_s"][0] == 0 and d["t_s"][-1] == pytest.approx(5/30000)
        assert np.all(d["R_m"] > 0)
        for observed, key in ((d["R_m"].min(), "Rmin_m"),
                              (d["R_m"].max(), "Rmax_m"),
                              (abs(d["U_m_s"]).max(), "max_abs_U_m_s")):
            assert observed == pytest.approx(float(row[key]), rel=1e-12, abs=1e-20)


def test_saved_convergence_for_both_amplitudes():
    """Verify independent amplitude references and recompute printed errors."""
    import csv
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "results/convergence.csv"
    with path.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 10
    for Pa in (1.0e5, 1.3e5):
        selection = [r for r in rows if float(r["Pa_Pa"]) == Pa]
        assert len(selection) == 5
        reference = next(r for r in selection if r["method"] == "DOP853" and float(r["rtol"]) == 1e-11)
        for row in selection:
            assert row["completed"] == "True"
            for key in ("Rmax_m", "Rmin_m", "collapse_time_s", "max_abs_U_m_s"):
                expected = abs(float(row[key])-float(reference[key]))/abs(float(reference[key]))
                assert float(row["relative_error_"+key]) == pytest.approx(expected, abs=1e-18)
                if float(row["rtol"]) <= 1e-10:
                    assert expected < 1e-6


def test_report_dependencies_and_english_sources():
    """Ensure the published report has its generated tables and all figures."""
    import re
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    tex = root / "report/mini_research_01.tex"
    text = tex.read_text()
    assert not re.search(r"[А-Яа-яЁё]", text)
    assert "stage not closed" not in text and "remains open" not in text
    dependencies = re.findall(r"\\input\{([^}]+)\}", text)
    assert len(dependencies) == 4
    for name in dependencies:
        dependency = tex.parent / name
        assert dependency.is_file(), name
        assert not re.search(r"[А-Яа-яЁё]", dependency.read_text())
    images = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text)
    assert len(images) == 9
    for name in images:
        assert (root / "figures" / name).is_file(), name
        assert (root / "figures" / Path(name).with_suffix(".png")).is_file(), name
