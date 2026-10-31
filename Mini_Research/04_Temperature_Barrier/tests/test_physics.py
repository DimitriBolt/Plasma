import hashlib
import json
import math
import unittest
import numpy as np
import pandas as pd
from scipy.integrate import quad
from src.upstream import ROOT,REPO,Parameters,RU,old
from src.thermodynamics import Mixture,T0,volume,species_cp,DATA
from src.model import acceleration,solve,Settings
from src.diagnostics import evaluate
from src.run_experiments import CONFIGS,design,run_key

class ThermodynamicsTests(unittest.TestCase):
    def test_initial_mechanical_equilibrium(self):
        p=Parameters()
        for gas in ('air','argon'):
            for c in CONFIGS.values():
                m=Mixture(gas=gas,**c)
                self.assertAlmostEqual(m.pressure(p.R0,T0)+m.p_constant,p.p0+2*p.sigma/p.R0,places=8)
                self.assertAlmostEqual(acceleration(0,p.R0,0,T0,0,m,1),0,delta=1e-6)

    def test_partial_pressures_and_saturated_inventory(self):
        p=Parameters();m=Mixture();partial=m.partial_pressures(p.R0,T0)
        self.assertAlmostEqual(partial['H2O'],p.pv,places=9)
        self.assertAlmostEqual(partial['N2']+partial['O2'],p.pg0,places=9)
        self.assertAlmostEqual(sum(partial.values()),m.pressure(p.R0,T0),places=8)

    def test_eos_and_no_vapor_double_count(self):
        for convention,v in [('mixture',1),('legacy',0)]:
            m=Mixture(convention=convention,vapor_factor=v);R=2e-6;T=2000.
            self.assertAlmostEqual(m.pressure(R,T)*float(volume(R))/(m.nt*RU*T),1,places=14)
            self.assertEqual(m.p_constant>0,'H2O' not in m.fractions)
        with self.assertRaises(ValueError): Mixture(convention='legacy',vapor_factor=1)
        self.assertEqual(Mixture(vapor_factor=0).p_constant,0.)

    def test_cp_cv_identity_and_positivity(self):
        for species in ('N2','O2','H2O','Ar_ground'):
            upper=6000 if species=='H2O' else 20000
            for T in np.geomspace(200,upper,51):
                cp=species_cp(species,float(T));cv=cp-RU
                self.assertGreater(cv,0);self.assertAlmostEqual(cp-cv,RU,places=12)
        for gas in ('air','argon'):
            m=Mixture(gas=gas)
            for T in np.geomspace(200,6000,37): self.assertGreater(m.capacity(float(T)),0)

    def test_energy_derivative_integral_and_joins(self):
        for gas in ('air','argon'):
            for config in ('M1_variable','M2_sat','M2_sat_const'):
                m=Mixture(gas=gas,**CONFIGS[config])
                for T in (210.,293.15,700.,1500.,5500.):
                    h=T*1e-5;der=(m.du(T+h)-m.du(T-h))/(2*h)
                    self.assertAlmostEqual(der/m.capacity(T),1,delta=1e-8)
                    numerical=quad(lambda x:m.capacity(x),T0,T,points=[1000.] if T>1000 else None,epsabs=1e-22)[0]
                    self.assertAlmostEqual(m.du(T)/numerical if numerical else 1,1,delta=1e-9)
                for boundary in (1000.,6000.):
                    if boundary<m.Tmax:
                        h=1e-5
                        self.assertLess(abs(m.du(boundary+h)-m.du(boundary-h)),3*h*m.capacity(boundary))
                samples=np.array([210.,293.15,700.,1500.,5500.])
                np.testing.assert_allclose(m.cv_array(samples),[m.cv(float(t)) for t in samples],rtol=1e-13)
                np.testing.assert_allclose(m.du_array(samples),[m.du(float(t)) for t in samples],rtol=1e-12,atol=1e-22)

    def test_pressure_derivative(self):
        for gas in ('air','argon'):
            m=Mixture(gas=gas);R=2e-6;U=-120.;T=1800.;dT=m.temperature_rate(R,U,T,1)
            dt=1e-13
            fd=(m.pressure(R+U*dt,T+dT*dt)-m.pressure(R-U*dt,T-dT*dt))/(2*dt)
            self.assertAlmostEqual(fd/m.pressure_rate(R,U,T,dT),1,delta=1e-7)

    def test_heat_flux_signs(self):
        m=Mixture()
        self.assertLess(m.heat_loss(6e-6,200,1),0)
        self.assertGreater(m.heat_loss(6e-6,1000,1),0)
        self.assertEqual(m.heat_loss(6e-6,T0,1),0)
        self.assertEqual(m.heat_loss(6e-6,1000,0),0)

    def test_mr3_exact_equation_limit(self):
        bt=old('bubble_thermo');co=old('coupled_model')
        for gas,gamma,kg in [('air',1.4,.0264),('argon',5/3,.0177)]:
            m=Mixture(gas=gas,**CONFIGS['M0']);g=bt.Gas(gamma=gamma,kg=kg)
            for R,U,T in [(6e-6,0,T0),(2e-6,-300,4000),(12e-6,20,230)]:
                rate=m.temperature_rate(R,U,T,1)
                self.assertAlmostEqual(rate/bt.temperature_rate(R,U,T,-m.heat_loss(R,T,1),g) if rate else 1,1,delta=1e-14)
                a=acceleration(1.1e-5,R,U,T,1.3,m,1);b=co.acceleration(1.1e-5,R,U,T,1.3e5,g,1,Parameters())
                self.assertAlmostEqual(a/b,1,delta=2e-14)

    def test_argon_no_artificial_variable_heat_capacity(self):
        m=Mixture(gas='argon',cv_mode='variable',vapor_factor=0,convention='legacy')
        for T in (100,1000,10000,1e6): self.assertAlmostEqual(m.cv(T)/(1.5*RU),1,places=14)

    def test_domain_stops_and_preserves_observed_interval(self):
        from dataclasses import replace
        m=Mixture();object.__setattr__(m,'Tmax',294.)
        s=solve(1.,m,mode='coupled');summary,z=evaluate(s)
        self.assertFalse(summary['completed']);self.assertTrue(math.isnan(summary['global_Tmax_K']))
        self.assertEqual(s.message,'coefficient_upper_boundary');self.assertAlmostEqual(z['T_K'][-1],294,places=7)
        self.assertLess(s.t[-1],5/Parameters().f)

class SavedResultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        file=ROOT/'results/summary.csv'
        if not file.exists(): raise unittest.SkipTest('Run production before artifact tests')
        cls.df=pd.read_csv(file,float_precision='round_trip')

    def test_design_complete_and_consistent(self):
        self.assertEqual(set(self.df.run_id),{run_key(*x) for x in design()})
        for row in self.df.to_dict('records'):
            m=json.loads((ROOT/'results/trajectories'/f"{row['run_id']}.json").read_text())
            self.assertAlmostEqual(row['observed_Tmax_K'],m['observed_Tmax_K'],places=7)
            with np.load(ROOT/'results/trajectories'/f"{row['run_id']}.npz") as z:
                j=np.argmax(z['T_K']);self.assertAlmostEqual(z['t_s'][j],row['t_Tmax_s'],places=18)
                self.assertAlmostEqual(z['T_K'][j],row['observed_Tmax_K'],places=7)
                self.assertLessEqual(z['T_K'].max(),row['Tmax_domain_K']+1e-6)
                self.assertGreaterEqual(z['T_K'].min(),row['Tmin_domain_K']-1e-6)
                if row['completed']: self.assertAlmostEqual(z['t_s'][-1],5/Parameters().f,places=16)
                else: self.assertTrue(math.isnan(row['global_Tmax_K']))

    def test_energy_ledgers(self):
        self.assertLess(self.df.energy_residual_scaled.max(),3e-5)
        self.assertTrue((self.df.source_energy_J==0).all())
        self.assertTrue((self.df.chemical_energy_change_J==0).all())
        for row in self.df.to_dict('records'):
            residual=row['delta_U_peak_J']+row['work_by_gas_peak_J']+row['signed_heat_loss_peak_J']
            self.assertLessEqual(abs(residual),row['energy_residual_max_J']*1.001+1e-20)

    def test_regression_and_convergence(self):
        reg=pd.read_csv(ROOT/'results/mr3_regression.csv',float_precision='round_trip');conv=pd.read_csv(ROOT/'results/convergence.csv',float_precision='round_trip')
        self.assertEqual(len(conv),56)
        for row in conv.to_dict('records'):
            self.assertEqual(row['completed'],bool(self.df[self.df.run_id==row['run_id']].iloc[0].completed))
        self.assertEqual(len(reg),12);self.assertLess(abs(reg.relative_difference).max(),2e-5)
        self.assertLess(abs(conv[conv.check.isin(['strict','Radau'])].peak_relative_difference).max(),2e-5)
        self.assertLess(abs(conv[conv.check.str.startswith('input_')].peak_relative_difference).max(),2e-4)
        self.assertTrue((conv.energy_residual_scaled<3e-5).all())
        q=pd.read_csv(ROOT/'results/independent_quadrature.csv',float_precision='round_trip')
        self.assertLess(q[q.order==16].scaled_residual.max(),3e-5)
        for row in conv[abs(conv.peak_time_difference_s)>1e-9].to_dict('records'):
            base=self.df[self.df.run_id==row['run_id']].iloc[0]
            self.assertGreater(base.near_peak_candidate_count,1)
            self.assertTrue(base.near_peak_time_first_s-1e-9<=base.t_Tmax_s+row['peak_time_difference_s']<=base.near_peak_time_last_s+1e-9)

    def test_extrema_stationary_or_boundary(self):
        for r in self.df.to_dict('records'):
            with np.load(ROOT/'results/trajectories'/f"{r['run_id']}.npz") as z:
                candidates=z['temperature_stationary_times_s']
                at_boundary=min(abs(r['t_Tmax_s']),abs(r['t_Tmax_s']-r['end_s']))<1e-18
                self.assertTrue(at_boundary or np.min(abs(candidates-r['t_Tmax_s']))<1e-18)
                if not at_boundary:
                    j=np.searchsorted(z['t_s'],r['t_Tmax_s']);mix=Mixture(gas=r['gas'],**CONFIGS[r['config']])
                    R,U,T=[float(z[k][j]) for k in ('R_m','U_m_s','T_K')]
                    compression=-mix.pressure(R,T)*4*np.pi*R**2*U;heat=mix.heat_loss(R,T,r['eta'])
                    self.assertLess(abs(compression-heat)/max(abs(compression)+abs(heat),1e-30),5e-5)

    def test_report_tables_use_saved_data(self):
        table=pd.read_csv(ROOT/'results/main_comparison.csv',float_precision='round_trip')
        for r in table.to_dict('records'):
            m=self.df[self.df.run_id==r['MR4_run_id']].iloc[0]
            self.assertAlmostEqual(m.observed_Tmax_K,r['MR4_observed_Tmax_K'],places=7)
            self.assertEqual(m.completed,r['MR4_completed'])
        text=(ROOT/'results/table_main.tex').read_text()
        for v in table.MR4_observed_Tmax_K: self.assertIn(f'{v:.1f}',text)

    def test_upstream_unchanged(self):
        inventory=json.loads((ROOT/'results/input_inventory.json').read_text())
        for f in inventory['files']:
            self.assertEqual(hashlib.sha256((REPO/f['path']).read_bytes()).hexdigest(),f['sha256'])

if __name__=='__main__': unittest.main()
