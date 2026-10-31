"""Independent algebra, event/domain behavior, and saved-data verification."""
import csv
import gzip
import hashlib
import json
from dataclasses import replace
import numpy as np
import pytest
from src.model import (ROOT, RP_PATH, Parameters, Settings, rp, km_rhs,
                       pressure_jump, pressure_derivative, denominator,
                       implicit_residual, gas_pressure, solve, summarize)


def test_reference_code_is_accepted_byte_for_byte():
    assert hashlib.sha256(RP_PATH.read_bytes()).hexdigest() == (
        'f8c1baaf5a4fc891383b7203b67c026908f4bf005d5418e273b640784e161188')


def test_balance_phase_and_initial_acceleration():
    p=Parameters(); pa=1.3e5
    assert abs(pressure_jump(0,p.R0,0,pa,p)) < 1e-9
    t=np.array([0,p.T/4,p.T/2,3*p.T/4])
    np.testing.assert_allclose(rp.external_pressure(t,pa,p),
                               [p.p0,p.p0-pa,p.p0,p.p0+pa],atol=1e-9)
    expected=p.R0*pa*2*np.pi*p.f/(p.rho*p.c*p.R0+4*p.mu)
    assert km_rhs(0,[p.R0,0],pa,p)[1] == pytest.approx(expected,rel=1e-13)
    assert abs(rp.rhs(0,[p.R0,0],pa,p)[1]) < 1e-6


def test_full_pressure_derivative_against_centered_difference():
    p=Parameters(); t=.37*p.T; R=2e-6; U=-75.; A=3e10; pa=1.2e5
    h=1e-13
    # Follow an arbitrary quadratic path; the test does not use the ODE.
    def jump(dt):
        return pressure_jump(t+dt,R+U*dt+.5*A*dt**2,U+A*dt,pa,p)
    numeric=(jump(h)-jump(-h))/(2*h)
    analytic=pressure_derivative(t,R,U,A,pa,p)
    assert numeric == pytest.approx(analytic,rel=2e-8)
    # Isolate the viscous derivative; this catches the missing acceleration term.
    def visc(dt): return -4*p.mu*(U+A*dt)/(R+U*dt+.5*A*dt**2)
    assert (visc(h)-visc(-h))/(2*h) == pytest.approx(
        -4*p.mu*A/R+4*p.mu*U**2/R**2,rel=2e-8)


@pytest.mark.parametrize('R,U,t',[(6e-6,0.,0.),(2e-6,-500.,8e-6),(8e-7,250.,16e-6)])
def test_explicit_acceleration_satisfies_original_implicit_equation(R,U,t):
    p=Parameters(); pa=1.3e5; a=km_rhs(t,[R,U],pa,p)[1]
    residual=implicit_residual(t,R,U,a,pa,p)
    scale=abs(R*a)+U**2+abs(pressure_jump(t,R,U,pa,p)/p.rho)+1
    assert abs(residual)/scale < 1e-14


def test_dimensional_scale_invariance():
    # Rescale base units numerically: L'=a L, T'=b T, M'=d M.
    # This is an independent dimensional test of EVERY term in acceleration.
    p=Parameters(); a,b,d=3.7,2.3,5.1
    q=replace(p,R0=p.R0*a,f=p.f/b,p0=p.p0*d/(a*b*b),
              rho=p.rho*d/a**3,mu=p.mu*d/(a*b),
              sigma=p.sigma*d/b**2,pv=p.pv*d/(a*b*b),c=p.c*a/b)
    t,R,U,pa=1e-5,3e-6,-72.,1.2e5
    x=km_rhs(t,[R,U],pa,p)[1]
    y=km_rhs(t*b,[R*a,U*a/b],pa*d/(a*b*b),q)[1]
    assert y == pytest.approx(x*a/b**2,rel=1e-13)
    assert denominator(R*a,U*a/b,q)==pytest.approx(a*denominator(R,U,p),rel=1e-14)


def test_algebraic_incompressible_limit():
    p=Parameters(c=1e30)
    for t,R,U in [(0,6e-6,0),(7e-6,2e-6,-200.),(1e-5,1e-6,500.)]:
        np.testing.assert_allclose(km_rhs(t,[R,U],1.3e5,p),
                                   rp.rhs(t,[R,U],1.3e5,p),rtol=1e-14,atol=1e-6)


def test_equilibrium_trajectory():
    p=Parameters(); s=solve(0); r=summarize(s)
    assert r['completed'] and r['local_minima_count']==0
    assert np.max(abs(s.y[0]/p.R0-1)) < 1e-12
    assert np.max(abs(s.y[1])) < 1e-8


def test_denominator_zero_is_outward_sonic_not_inward_mach_one():
    p=Parameters(); R=1e-6
    Uzero=p.c+4*p.mu/(p.rho*R)
    assert abs(denominator(R,Uzero,p))/R < 1e-15
    assert denominator(R,-p.c,p)>2*R
    assert denominator(R,p.c,p)>0
    with pytest.raises(ValueError,match='denominator'):
        solve(1e5,y0=[R,Uzero])


def test_incomplete_run_keeps_last_point_and_censors_metrics():
    # Deliberately high technical radius guard, no modified pressure or velocity.
    cfg=Settings(km_radius_guard_m=5.9e-6)
    sol=solve(1e5,settings=cfg); r=summarize(sol)
    assert not r['completed'] and sol.status==1
    assert sol.t[-1]<5*Parameters().T
    assert sol.y[0,-1]==pytest.approx(cfg.km_radius_guard_m,rel=1e-12)
    assert r['extrema_scope']=='observed_interval_censored'
    assert r['last_time_s']==sol.t[-1]


@pytest.mark.parametrize('R',[0.,-1e-6])
def test_no_radius_clipping(R):
    with pytest.raises(ValueError): km_rhs(0,[R,0],1e5,Parameters())


def test_reproduce_from_saved_parameters_and_independent_events():
    path=ROOT/'results/trajectories/KM_1.300bar.json'
    if not path.exists(): pytest.skip('Run experiments first')
    data=json.loads(path.read_text()); p=Parameters(**data['parameters'])
    cfg=Settings(**data['settings']); sol=solve(1.3e5,p=p,settings=cfg)
    row=summarize(sol)
    for k in ('Rmax_m','Rmin_m','max_abs_U_m_s','collapse_time_s'):
        assert row[k]==pytest.approx(data['summary'][k],rel=1e-10)
    for e in sol.events:
        t=e['time_s']; h=max(4*np.spacing(t),1e-16)
        if e['kind']=='minimum':
            assert sol.sol(t-h)[1]<0 and sol.sol(t+h)[1]>0
        elif e['kind']=='maximum':
            assert sol.sol(t-h)[1]>0 and sol.sol(t+h)[1]<0


def csv_rows(name):
    with (ROOT/'results'/name).open() as f: return list(csv.DictReader(f))


def test_every_saved_trajectory_contains_all_summary_extrema_and_positive_radius():
    path=ROOT/'results/amplitude_summary.csv'
    if not path.exists(): pytest.skip('Run sweep first')
    for row in csv_rows('amplitude_summary.csv'):
        key=f"{row['model']}_{float(row['Pa_bar']):.3f}bar"
        with gzip.open(ROOT/'results/trajectories'/f'{key}.csv.gz','rt') as f:
            x=np.genfromtxt(f,delimiter=',',names=True)
        assert np.all(x['R_m']>0) and np.all(np.diff(x['t_s'])>0)
        if row['completed']=='True':
            assert x['t_s'][-1]==pytest.approx(5/30000,rel=1e-14)
        else:
            assert row['extrema_scope']=='observed_interval_censored'
            assert x['t_s'][-1]<5/30000
            assert row['solver_message']!='Reached requested end time'
        for val,k in [(min(x['R_m']),'Rmin_m'),(max(x['R_m']),'Rmax_m'),
                      (max(abs(x['U_m_s'])),'max_abs_U_m_s')]:
            assert val==pytest.approx(float(row[k]),rel=2e-12,abs=1e-20)


def test_saved_RP_reproduction_and_convergence():
    if not (ROOT/'results/convergence.csv').exists(): pytest.skip('Run convergence first')
    for r in csv_rows('rp_reproduction.csv'):
        if float(r['Pa_bar'])>0:
            assert float(r['relative_error']) < 2e-7
    for r in csv_rows('convergence.csv'):
        assert r['completed']=='True'
        if float(r['rtol'])<=1e-10:
            for k in ('Rmin_m','Rmax_m','max_abs_U_m_s','collapse_time_s','collapse_duration_s'):
                assert float(r['relative_error_'+k]) < (1e-4 if r['model']=='RP' else 2e-5)


def test_saved_c_limit():
    if not (ROOT/'results/c_limit.csv').exists(): pytest.skip('Run c-limit first')
    for bar in (1.,1.3):
        group=[r for r in csv_rows('c_limit.csv') if float(r['Pa_bar'])==bar]
        assert len(group)==6
        assert float(group[-1]['relative_Rmin_error'])<float(group[-2]['relative_Rmin_error'])
        assert float(group[-1]['relative_Rmin_error'])<.01
        assert float(group[-1]['relative_max_speed_error'])<.02


def test_generated_comparisons_do_not_mix_cycles():
    if not (ROOT/'results/cycle_comparison.csv').exists(): pytest.skip('Run sweep first')
    p=Parameters()
    for row in csv_rows('cycle_comparison.csv'):
        for m in ('RP','KM'):
            t=float(row[m+'_time_s'])
            if np.isfinite(t): assert int(t/p.T)+1==int(row['cycle'])


def test_censored_run_differences_are_suppressed():
    if not (ROOT/'results/comparison.csv').exists(): pytest.skip('Run sweep first')
    for r in csv_rows('comparison.csv'):
        if r['RP_completed']=='False' or r['KM_completed']=='False':
            assert np.isnan(float(r['relative_difference_Rmin_m']))
            assert np.isnan(float(r['KM_minus_RP_collapse_time_s']))


def test_report_has_generated_tables_and_english_figures():
    import re
    text=(ROOT/'report/main.tex').read_text()
    assert not re.search('[А-Яа-яЁё]',text)
    for name in re.findall(r'\\input\{([^}]+)\}',text):
        path=ROOT/'report'/name
        assert path.is_file(),name
        assert not re.search('[А-Яа-яЁё]',path.read_text())
    figures=re.findall(r'\\fig\{([^}]+)\}',text)
    assert len(figures)==13
    for name in figures:
        for ext in ('pdf','png'): assert (ROOT/'figures'/f'{name}.{ext}').is_file()


def test_extended_limit_and_independent_clean_repetition():
    if not (ROOT/'results/clean_reproduction.csv').exists(): pytest.skip('Run extended checks')
    last=csv_rows('c_limit_extended.csv')[-1]
    assert float(last['relative_Rmin_error'])<1e-8
    assert float(last['max_R_difference_over_R0'])<.01
    assert float(last['max_U_difference_over_reference_max'])<.1
    for row in csv_rows('clean_reproduction.csv'):
        if float(row['Pa_bar'])>0:
            assert row['completed']=='True'
            assert float(row['relative_error'])<1e-6
