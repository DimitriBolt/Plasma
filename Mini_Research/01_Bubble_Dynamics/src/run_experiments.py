"""Reproduce the amplitude sweep, numerical checks, reference audit, and tables."""
import argparse
import csv
import hashlib
import json
import platform
import subprocess
from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
import numpy as np
from iapws import IAPWS95
from .bubble_dynamics import (Parameters, event_times, external_pressure,
                              integrate, linear_solution, sample_times, summarize)

ROOT = Path(__file__).resolve().parents[1]
AMPLITUDES = (0, .2, .5, .8, 1., 1.3)
SOURCES = {
    "IAPWS95": "https://iapws.org/technical-guidance/release/IAPWS-95",
    "surface_tension": "https://iapws.org/technical-guidance/release/Surf-H2O",
    "viscosity": "https://iapws.org/technical-guidance/release/viscosity",
    "AIP2016": "https://doi.org/10.1063/1.4945343",
    "Nature2005": "https://doi.org/10.1038/nature03361",
    "solve_ivp": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html",
}


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reference_audit():
    water = IAPWS95(T=293.15, P=.101325)
    saturation = IAPWS95(T=293.15, x=0)
    return dict(implementation="iapws 1.5.4; IAPWS95 and IAPWS transport correlations",
                temperature_K=293.15, pressure_MPa=.101325,
                rho_kg_m3=water.rho, mu_Pa_s=water.mu,
                c_m_s=water.w, sigma_N_m=saturation.sigma,
                saturation_pressure_Pa=saturation.P * 1e6)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path,
                        default=ROOT.parents[1], help="Repository containing Library/ and Plasma-plan.tex")
    args = parser.parse_args()
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    p = Parameters()
    summaries, events, solver_runs = [], [], []
    with (out / "trajectories.csv").open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["Pa_Pa", "t_s", "R_m", "U_m_s", "p_infty_Pa", "wall_Mach"])
        for bar in AMPLITUDES:
            Pa = bar * 1e5
            sol = integrate(Pa, p)
            metrics = summarize(sol, Pa, p)
            if not metrics["completed"]:
                raise RuntimeError(f"Incomplete sweep at {bar} bar: {metrics}")
            summaries.append(metrics)
            t = sample_times(sol, p)
            R, U = sol.sol(t)
            if np.any(R <= 0):
                raise RuntimeError("Nonpositive exported radius")
            writer.writerows(zip(np.full(len(t), Pa), t, R, U,
                                 external_pressure(t, Pa, p), abs(U)/p.c))
            if Pa:
                for index, kind in ((1,"radius_maximum"),(2,"radius_minimum"),(3,"velocity_extremum")):
                    for te in event_times(sol, index, p):
                        re, ue = sol.sol(te)
                        events.append(dict(Pa_Pa=Pa, kind=kind, t_s=te, R_m=re, U_m_s=ue))
            solver_runs.append(dict(Pa_Pa=Pa, nfev=sol.nfev, accepted_steps=len(sol.t)-1,
                                    min_step_s=float(np.diff(sol.t).min()),
                                    max_step_s=float(np.diff(sol.t).max()),
                                    export_rows=len(t), status=sol.status))
            print(f"{bar:.1f} bar: R/R0=[{metrics['Rmin_over_R0']:.6g}, "
                  f"{metrics['Rmax_over_R0']:.6g}], |U|max={metrics['max_abs_U_m_s']:.6g}, "
                  f"Mw={metrics['max_wall_Mach']:.6g}", flush=True)
    write_csv(out / "amplitude_summary.csv", summaries)
    write_csv(out / "extrema.csv", events)

    convergence = []
    settings = [("DOP853", 1e-7, (1e-13,1e-5),100),
                ("DOP853", 1e-9, (1e-15,1e-7),300),
                ("DOP853", 1e-10,(1e-16,1e-8),1000),
                ("DOP853", 1e-11,(1e-17,1e-9),2000),
                ("Radau", 1e-10,(1e-16,1e-8),1000)]
    for Pa in (1.0e5, 1.3e5):
        amplitude_rows = []
        for method, rtol, atol, spp in settings:
            sol = integrate(Pa, p, method=method, rtol=rtol, atol=atol, steps_per_period=spp)
            row = summarize(sol, Pa, p)
            if not row["completed"]:
                raise RuntimeError(f"Convergence integration failed at {Pa/1e5:g} bar")
            amplitude_rows.append(dict(method=method, rtol=rtol, atol_R_m=atol[0],
                                       atol_U_m_s=atol[1], max_step_s=p.T/spp, **row))
            print(f"Convergence: {Pa/1e5:g} bar, {method} rtol={rtol:g}, nfev={sol.nfev}", flush=True)
        reference = amplitude_rows[3]
        for row in amplitude_rows:
            for key in ("Rmax_m", "Rmin_m", "collapse_time_s", "max_abs_U_m_s"):
                row["relative_error_"+key] = abs(row[key]-reference[key])/abs(reference[key])
        convergence.extend(amplitude_rows)
    write_csv(out / "convergence.csv", convergence)

    eq = integrate(0, p)
    eq_y = eq.sol(sample_times(eq,p))
    eq_error = float(np.max(abs((eq_y[0]-p.R0)/p.R0)))
    linear_checks, linear_rows = [], []
    t = np.linspace(0,5*p.T,20001)
    for Pa in (100.,500.):
        nonlinear = integrate(Pa,p).sol(t)
        linear = linear_solution(t,Pa,p)
        error = float(np.max(abs(nonlinear[0]-linear[0]))/p.R0)
        signal = float(np.max(abs(linear[0]-p.R0))/p.R0)
        linear_checks.append(dict(Pa_Pa=Pa,max_radius_error_over_R0=error,
                                  max_linear_displacement_over_R0=signal,
                                  error_over_linear_displacement=error/signal))
        linear_rows.extend(dict(Pa_Pa=Pa,t_s=ti,R_nonlinear_m=rn,R_linear_m=rl)
                           for ti,rn,rl in zip(t,nonlinear[0],linear[0]))
    write_csv(out / "linear_comparison.csv", linear_rows)
    source_hashes = {}
    for name in ("Library/035218_1_online.pdf","Library/nature03361.pdf","Plasma-plan.tex"):
        path = args.source_root / name
        if not path.is_file():
            raise FileNotFoundError(f"Required source missing: {path}")
        source_hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    try:
        commit = subprocess.check_output(["git","-C",str(args.source_root),"rev-parse","HEAD"],
                                         text=True,stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    code_hashes = {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in sorted((ROOT / "src").glob("*.py"))}
    metadata = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                    scope="Mini-Research 1 only; radius dynamics; no temperature model",
                    python=platform.python_version(), platform=platform.platform(),
                    packages={name:version(name) for name in
                              ("numpy","scipy","matplotlib","pytest","iapws")},
                    parameters_SI=asdict(p),derived=dict(pg0_Pa=p.pg0,period_s=p.T,
                          natural_frequency_Hz=p.omega0/(2*np.pi),R_floor_m=p.R_floor),
                    production_solver=dict(method="DOP853",rtol=1e-10,atol=[1e-16,1e-8],
                          max_step_s=p.T/1000,cycles=5),
                    sampling="5001 uniform times plus four subdivisions per accepted adaptive interval plus event roots",
                    metric_definition="Global extrema over [0,5T] include endpoints; local counts exclude t=0. Collapse time is deepest post-initial local minimum; duration from preceding local maximum.",
                    convergence_reference="Per amplitude: DOP853 rtol=1e-11, atol=[1e-17,1e-9], max_step=T/2000",
                    convergence_amplitudes_Pa=[1.0e5, 1.3e5],
                    equilibrium=dict(relative_radius_error=eq_error,
                          max_abs_U_m_s=float(abs(eq_y[1]).max()),
                          pressure_residual_Pa=p.pg0+p.pv-p.p0-2*p.sigma/p.R0),
                    pressure_test=dict(Pa_Pa=1e5,t_s=[0,p.T/4,3*p.T/4],
                          values_Pa=external_pressure(np.array([0,p.T/4,3*p.T/4]),1e5,p).tolist()),
                    linear_checks=linear_checks,reference_properties=reference_audit(),
                    source_urls=SOURCES,source_sha256=source_hashes,
                    source_repository_commit=commit,code_sha256=code_hashes,solver_runs=solver_runs)
    (out / "run_metadata.json").write_text(json.dumps(metadata,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(equilibrium=metadata["equilibrium"],linear_checks=linear_checks,
                          reference_properties=metadata["reference_properties"]),indent=2))


if __name__ == "__main__":
    main()
