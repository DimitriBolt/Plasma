"""Physical invariants and limiting cases; no tests reproduce implementation algebra."""
from dataclasses import replace
from functools import lru_cache
import json
import numpy as np
import pandas as pd
import pytest
from scipy.integrate import solve_ivp,quad
from src.coupling import (Geometry,FrozenHeatModel,ThermalHistory,thermal_source,particle_trajectory,
                          gaussian_source,point_sphere_fraction,point_robin_fraction,integrated_parts)
from src.upstream import ROOT,optics,thermo

@lru_cache(None)
def heat(h=4400.,G=None):
    model=FrozenHeatModel(Geometry(gap_m=100e-9),h,refinement=.75,interface_W_m2K=G)
    t=np.unique(np.r_[0.,np.geomspace(1e-12,2e-6,500),10e-9])
    return model,model.pulse(1e-15,10e-9,t)

def test_geometry_contact_rejected():
    for d in [0.,-1e-9,np.nan]:
        with pytest.raises(ValueError): Geometry(gap_m=d)

def test_zero_intensity_exact():
    m,_=heat();o=m.pulse(0.,10e-9,np.array([0.,1e-8,1e-7]))
    assert np.array_equal(o['q_wall_laser_increment_W'],np.zeros(3))
    assert np.array_equal(o['gold_energy_J'],np.zeros(3))

def test_independent_energy_channels():
    _,o=heat();E=1e-15
    assert max(abs(o['balance_residual_J']))<E*2e-8
    assert o['gas_energy_J'][-1]+o['bath_energy_J'][-1]<E
    # Integrate the saved wall flux independently of the ODE's cumulative state.
    assert abs(np.trapezoid(o['q_wall_laser_increment_W'],o['time_s'])/o['gas_energy_J'][-1]-1)<.003

def test_SI_diffusive_scale_invariance():
    m,o=heat();g=m.geometry;s=2.
    big=FrozenHeatModel(Geometry(radius_m=g.radius_m*s,particle_radius_m=g.particle_radius_m*s,
                                 gap_m=g.gap_m*s,outer_padding_m=g.outer_padding_m*s),
                        h_W_m2K=m.h/s,refinement=.75,interface_W_m2K=m.G/s)
    b=big.pulse(1e-15*s**3,10e-9*s**2,o['time_s']*s**2)
    assert max(abs(b['gold_rise_K']-o['gold_rise_K']))<max(o['gold_rise_K'])*1e-5
    assert b['gas_energy_J'][-1]/s**3==pytest.approx(o['gas_energy_J'][-1],rel=1e-5)

def test_zero_gas_heat_transfer():
    _,o=heat(0.)
    assert np.array_equal(o['q_wall_laser_increment_W'],np.zeros(len(o['time_s'])))
    assert np.array_equal(o['gas_energy_J'],np.zeros(len(o['time_s'])))

def test_zero_Au_water_exchange():
    m,o=heat(4400.,0.)
    assert np.allclose(o['gold_energy_J'],o['absorbed_so_far_J'],rtol=1e-6,atol=1e-25)
    assert max(abs(o['water_energy_J']))<1e-28
    assert max(abs(o['gas_energy_J']))<1e-28

def test_stronger_contact_transmits_more():
    assert heat(1e8)[1]['gas_energy_J'][-1]>heat()[1]['gas_energy_J'][-1]>heat(1e3)[1]['gas_energy_J'][-1]

def test_passive_invariant_against_ODE():
    # Prescribed linear collapse is a nontrivial independent differential control.
    R0,R1=10e-6,1e-6;end=1e-6;U=(R1-R0)/end
    t=np.linspace(0,end,301);R=R0+U*t
    r,d,status=particle_trajectory(R,20e-9,100e-9,'passive')
    sol=solve_ivp(lambda s,y:[(R0+U*s)**2*U/y[0]**2],(0,end),[r[0]],rtol=2e-12,atol=1e-18,t_eval=t)
    assert max(abs(sol.y[0]-r))<5e-15
    assert np.ptp(r**3-R**3)<1e-28
    assert np.all(np.diff(d)>0) and np.all(status=='separated')

def test_fixed_contact_separate_status():
    _,d,status=particle_trajectory([6e-6,7e-6],20e-9,10e-9,'fixed')
    assert d[-1]<0 and status[-1]=='contact_or_intersection'

def test_capture_is_prescribed_gap():
    _,d,_=particle_trajectory([10e-6,1e-6],20e-9,10e-9,'capture')
    assert np.allclose(d,10e-9,rtol=1e-12,atol=1e-22)

def test_gaussian_finite_interval_normalization():
    interval=(-3e-9,1e-9);width=interval[1]-interval[0]
    integral=quad(lambda x:float(gaussian_source(interval[0]+width*x,1e-15,0.,2e-9,interval))*width,0,1,epsabs=1e-28)[0]
    assert abs(integral/1e-15-1)<1e-12
    assert gaussian_source(-4e-9,1e-15,0.,2e-9,interval)==0

def test_signed_flux_not_clipped():
    history=ThermalHistory([0.,1.,2.],[0.,1.,0.],[-2.,-2.,1.])
    q=thermal_source(1.,history)
    assert q['q_wall_total_W']==-1 and q['q_wall_negative_W']==-1
    assert q['q_wall_positive_W']==0 and q['q_wall_laser_increment_W']==1
    parts=integrated_parts([0.,1.,2.],[-2.,-1.,1.])
    assert parts['negative_J']<0 and parts['positive_J']>0
    assert parts['signed_J']==pytest.approx(parts['positive_J']+parts['negative_J'])

def test_history_memory_serialization_and_domain(tmp_path):
    h=ThermalHistory([0.,1.,2.],[0.,2.,0.],[-1.,-1.,-1.],{'frozen_radius_m':6e-6})
    path=tmp_path/'source.npz';h.save(path);loaded=ThermalHistory.load(path)
    assert loaded.thermal_source(.7)==h.thermal_source(.7)
    for t in [-1e-10,2.+1e-10]:
        with pytest.raises(ValueError): loaded.thermal_source(t)
    with pytest.raises(ValueError): loaded.thermal_source(1.,{'R_m':5e-6})

def test_point_analytic_limits():
    R=6e-6;r=7e-6;alpha=.6/(1000*4187)
    assert point_sphere_fraction(np.array([0.]),R,r,alpha)[0]==0
    assert point_sphere_fraction(np.array([1e8]),R,r,alpha)[0]==pytest.approx(R/r,rel=1e-6)
    assert point_robin_fraction(np.array([1e-6]),R,r,alpha,.6,0.)[0]==0
    a=point_robin_fraction(np.array([1e-6]),R,r,alpha,.6,1e12)
    b=point_sphere_fraction(np.array([1e-6]),R,r,alpha)
    assert a[0]==pytest.approx(b[0],rel=1e-5)

def test_optical_units_and_energy():
    c=optics.optical_cross_sections(20e-9,532e-9,1.33)
    assert c.absorption_m2==pytest.approx(3.4337324474063743e-15,rel=2e-13)
    assert optics.absorbed_power(1e8,c.absorption_m2)*10e-9==pytest.approx(optics.pulse_energy(1.,c.absorption_m2),rel=1e-14)
    assert c.absorption_m2+c.scattering_m2==pytest.approx(c.extinction_m2,rel=1e-12)

def test_mixture_equilibrium_no_double_vapor():
    m=thermo.Mixture();p=thermo.Parameters()
    assert m.p_constant==0
    assert m.pressure(p.R0,293.15)==pytest.approx(p.p0+2*p.sigma/p.R0)
    assert 'H2O' in m.fractions

def test_coefficients_stop_at_real_boundary():
    m=thermo.Mixture()
    assert m.Tmax==6000.
    with pytest.raises(ValueError): m.capacity(6001.)
    assert m.capacity(6000.)>0

def test_saved_gap_delay_and_photon_budget():
    gap=pd.read_csv(ROOT/'results/gap_scan.csv')
    small=gap.iloc[gap.d_m.argmin()];far=gap.iloc[gap.d_m.argmax()]
    assert far.t_q_peak_s>small.t_q_peak_s
    assert far.gas_fraction<small.gas_fraction
    multi=pd.read_csv(ROOT/'results/multiparticle.csv')
    assert (multi.local_absorbed_J<=multi.incident_beam_J).all()
    assert (multi.gas_energy_J<=multi.local_absorbed_J).all()
    assert multi[multi.N_particles==0].gas_energy_J.iloc[0]==0
    assert multi.nonoverlap_screen.all()

def test_late_pulse_causality():
    tab=pd.read_csv(ROOT/'results/timing_scan.csv')
    after=tab[tab.delta_t_s>10e-9]
    assert (after.gas_by_tc_J==0).all()
    assert (after.gas_late_window_J==0).all()

def test_saved_field_matches_energy_and_maximum_principle():
    # Check the serialized, rescaled physical field as well as the solver ledger.
    z=np.load(ROOT/'results/reference.npz');m=FrozenHeatModel(Geometry(gap_m=100e-9))
    assert np.array_equal(m.xy,z['xy_m'])
    energy=m.full_capacity[:-1]@z['water_rise_K']
    assert max(abs(energy-z['water_energy_J']))<z['absorbed_so_far_J'][-1]*1e-10
    assert max(z['nearest_water_rise_K'])<=max(z['gold_rise_K'])*(1+1e-5)
    assert max(z['water_rise_K'].ravel())<=max(z['gold_rise_K'])*(1+1e-5)
    assert np.allclose(z['gold_energy_J'],m.full_capacity[-1]*z['gold_rise_K'],rtol=1e-12,atol=1e-28)
    paths=sorted((ROOT/'results/fields').glob('*.npz'))
    assert len(paths)==4
    for path in paths:
        with np.load(path) as f:
            meta=json.loads(str(f['metadata_json']));g=Geometry(**meta['geometry'])
            model=FrozenHeatModel(g,h_W_m2K=meta['h_W_m2K'])
            assert np.array_equal(model.xy,f['xy_m'])
            Ew=model.full_capacity[:-1]@f['water_rise_K'];E=f['absorbed_so_far_J'][-1]
            assert max(abs(Ew-f['water_energy_J']))<E*1e-10
            assert max(f['water_rise_K'].ravel())<=max(f['gold_rise_K'])*(1+1e-5)
            assert max(abs(f['balance_residual_J']))<E*1e-8

def test_finite_design_objectives_and_admission():
    tab=pd.read_csv(ROOT/'results/finite_design.csv');data=json.loads((ROOT/'results/design_extrema.json').read_text())
    admitted=tab[tab.moderate_static_candidate]
    assert len(tab)==216 and len(admitted)>0
    for key,column in [('moderate_max_absorbed','energy_absorbed_J'),('moderate_max_gas_energy','gas_energy_J'),('moderate_max_fraction','gas_fraction')]:
        assert data[key][column]==pytest.approx(admitted[column].max(),rel=2e-14)
    assert (admitted.max_gold_rise_K<=30.*(1+1e-12)).all()
    assert (admitted.gap_over_a>=5.).all()

def test_saved_numerical_validation():
    data=json.loads((ROOT/'results/validation_summary.json').read_text())
    assert data['independent_Radau_gold_relative']<2e-5
    assert abs(data['independent_Radau_gas_relative'])<2e-5
    assert data['energy_balance_relative_max']<1e-8
    assert data['isolated_gold_peak_finest_relative_max']<.02

def test_gas_censoring_and_signed_loss():
    tab=pd.read_csv(ROOT/'results/gas_diagnostic.csv');c=tab[~tab.completed_window]
    assert len(c)>0 and (c.censoring_temperature_K-6000.).abs().max()<1e-5
    assert (c.end_s<c.requested_end_s).all()
    assert tab.max_pressure_mismatch.max()>.01
    # The net total wall flux includes a baseline with either sign.
    assert (tab.negative_J<0).all()
    assert tab.min_delta_T_K.min()>-1e-12
    complete=tab[(tab.Pa_bar==1.15)&tab.completed_window]
    assert complete[complete.mechanism=='B_after'].delta_T_at_baseline_peak_K.iloc[0]==0.
    refinement=pd.read_csv(ROOT/'results/gas_increment_convergence.csv')
    assert refinement.relative_difference.abs().max()<.002

def test_upstream_fingerprints():
    import hashlib
    data=json.loads((ROOT/'results/input_inventory.json').read_text())
    assert data['accepted_snapshot']=='0204ce7a8429406b7ec235b02cec8581274fc8c6'
    for row in data['files']:
        assert hashlib.sha256((ROOT.parents[1]/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
