"""Equation identities, transport signs, independent limits and saved-data integrity."""
import hashlib
import json
import unittest
import numpy as np
import pandas as pd
from src.bubble_thermo import (ROOT,REPO,Gas,Parameters,RU,mechanics,volume,
    adiabatic_temperature,adiabatic_pressure,gas_pressure,temperature_rate,pressure_rate)
from src.heat_transfer import heat_loss,relaxation_time
from src.heat_source import heat_source
from src.coupled_model import acceleration,rhs,solve_coupled
from src.one_way import solve_one_way,temperature
from src.numerics import Settings,InvalidState,hot_intervals,thermal_extrema
from src.saved_mechanics import SavedMechanics
from src.diagnostics import adiabatic_metrics


class EquationTests(unittest.TestCase):
    def setUp(self):self.p=Parameters();self.g=Gas()

    def test_01_static_pressure_balance(self):
        p=self.p
        self.assertAlmostEqual(p.pg0+p.pv,p.p0+2*p.sigma/p.R0,places=9)
        self.assertAlmostEqual(p.pg0,123232.66666666666,places=8)

    def test_02_ideal_gas_and_adiabatic(self):
        p,g=self.p,self.g;R=p.R0*np.geomspace(.01,8,100)
        T=adiabatic_temperature(R,g,p);pg=adiabatic_pressure(R,g,p)
        np.testing.assert_allclose(pg*volume(R),g.n(p)*RU*T,rtol=4e-15)
        np.testing.assert_allclose(pg,gas_pressure(R,T,g,p),rtol=4e-15)
        np.testing.assert_allclose(T*volume(R)**(g.gamma-1),g.T0*volume(p.R0)**(g.gamma-1),rtol=4e-15)

    def test_03_pressure_derivative_chain_rule(self):
        p,g=self.p,self.g;R=2.3e-6;U=-54.;T=2300.;q=-.002
        dT=temperature_rate(R,U,T,q,g,p);pg=gas_pressure(R,T,g,p)
        chain=pg*(dT/T-3*U/R)
        self.assertAlmostEqual(pressure_rate(R,U,T,q,g,p)/chain,1.,places=13)

    def test_04_zero_heat_polytropic_derivative(self):
        p,g=self.p,self.g;R=2e-6;U=-100.;T=float(adiabatic_temperature(R))
        old=-3*p.kappa*mechanics.gas_pressure(R,p)*U/R
        self.assertAlmostEqual(pressure_rate(R,U,T,0,g,p)/old,1.,places=13)

    def test_05_heat_flow_sign(self):
        p,g=self.p,self.g
        self.assertGreater(heat_loss(p.R0,g.Tl+10,g,1),0)
        self.assertLess(heat_loss(p.R0,g.Tl-10,g,1),0)
        self.assertEqual(heat_loss(p.R0,g.Tl,g,1),0)
        self.assertEqual(heat_loss(p.R0,1e9,g,0),0)
        with self.assertRaises(ValueError):heat_loss(p.R0,300,g,-1)

    def test_06_source_is_zero(self):
        self.assertEqual(heat_source(0,1e-6,1e7),0)

    def test_07_km_implicit_viscous_derivative(self):
        p,g=self.p,self.g;t=13e-6;R=2e-6;U=-120.;T=3100.;Pa=1.3e5;eta=1.
        A=acceleration(t,R,U,T,Pa,g,eta,p)
        pg=gas_pressure(R,T,g,p);dpg=pressure_rate(R,U,T,-heat_loss(R,T,g,eta),g,p)
        dp=pg+p.pv-(p.p0-Pa*np.sin(2*np.pi*p.f*t))-2*p.sigma/R-4*p.mu*U/R
        ddp=dpg+2*p.sigma*U/R**2+Pa*2*np.pi*p.f*np.cos(2*np.pi*p.f*t)-4*p.mu*A/R+4*p.mu*U**2/R**2
        left=(1-U/p.c)*R*A+1.5*(1-U/(3*p.c))*U**2
        right=(1+U/p.c)*dp/p.rho+R*ddp/(p.rho*p.c)
        self.assertLess(abs(left-right)/max(abs(left),abs(right)),1e-13)

    def test_08_km_rhs_matches_accepted_equation(self):
        p,g=self.p,self.g
        for t,R,U in [(0,p.R0,0),(18e-6,1e-6,-650),(23e-6,8e-6,20)]:
            T=float(adiabatic_temperature(R))
            a=acceleration(t,R,U,T,1.3e5,g,0,p)
            ref=mechanics.km_rhs(t,[R,U],1.3e5,p)[1]
            np.testing.assert_allclose(a,ref,rtol=2e-14,atol=1e-6)

    def test_09_unforced_equilibrium(self):
        sol=solve_coupled(0,eta=1,settings=Settings('DOP853',1e-9,100),cycles=1.)
        self.assertTrue(sol.completed)
        np.testing.assert_allclose(sol.y[0],self.p.R0,rtol=1e-12)
        np.testing.assert_allclose(sol.y[2],self.g.T0,rtol=1e-12)
        self.assertLess(np.max(abs(sol.y[1])),1e-10)

    def test_10_one_way_zero_heat_limit(self):
        p,g=self.p,self.g
        class Path:
            end=1/p.f;amplitude=0.
            def __call__(self,t):
                return p.R0*(1+.2*np.sin(2*np.pi*p.f*np.asarray(t))),p.R0*.2*2*np.pi*p.f*np.cos(2*np.pi*p.f*np.asarray(t))
        path=Path();sol=solve_one_way(path,eta=0,settings=Settings('DOP853',1e-10,100))
        t=np.linspace(0,path.end,200)
        np.testing.assert_allclose(temperature(sol,t),adiabatic_temperature(path(t)[0]),rtol=1e-14)

    def test_11_signed_heat_relaxation_analytic(self):
        p=self.p;g=Gas(T0=280,Tl=310)
        class Path:
            end=2e-6;amplitude=0.
            def __call__(self,t):return p.R0+np.zeros_like(t),np.zeros_like(t)
        sol=solve_one_way(Path(),g,eta=1,settings=Settings('DOP853',1e-10,200))
        t=np.linspace(0,Path.end,100)
        expected=g.Tl+(g.T0-g.Tl)*np.exp(-t/relaxation_time(p.R0,g,1,p))
        np.testing.assert_allclose(temperature(sol,t),expected,rtol=1e-9)

    def test_12_no_state_clipping(self):
        for state in [[-1e-6,0,300,0,0,0],[1e-6,0,-10,0,0,0]]:
            with self.assertRaises(InvalidState):rhs(0,state,1e5)

    def test_13_peak_and_threshold_off_plot_grid(self):
        T=lambda t:300+10000*np.exp(-((np.asarray(t)-.531)/.0002)**2)
        dT=lambda t:-2*(np.asarray(t)-.531)/.0002**2*(T(t)-300)
        grid=np.array([.530,.5309,.5311,.532])
        peak,time,_=thermal_extrema(T,dT,grid)
        self.assertAlmostEqual(peak,10300,places=6);self.assertAlmostEqual(time,.531,places=14)
        intervals=hot_intervals(T,np.r_[grid,time],1000)
        expected=2*.0002*np.sqrt(np.log(10000/700))
        self.assertAlmostEqual(sum(b-a for a,b in intervals),expected,places=13)

    def test_14_dimensions(self):
        # Dimensions are powers of (kg,m,s,K,mol); symbolic addition represents multiplication.
        def mul(*xs):return tuple(sum(a) for a in zip(*xs))
        def powd(x,n):return tuple(n*a for a in x)
        R=(0,1,0,0,0);U=(0,1,-1,0,0);T=(0,0,0,1,0)
        pg=(1,-1,-2,0,0);kg=(1,1,-3,-1,0);capacity=(1,2,-2,-1,0)
        power=(1,2,-3,0,0);rho=(1,-3,0,0,0);c=U
        self.assertEqual(mul(kg,R,T),power)
        self.assertEqual(mul(pg,powd(R,2),U),power)
        self.assertEqual(mul(power,powd(capacity,-1)),(0,0,-1,1,0))
        self.assertEqual(mul(U,powd(R,-1),T),(0,0,-1,1,0))
        self.assertEqual(mul(pg,U,powd(R,-1)),(1,-1,-3,0,0))
        self.assertEqual(mul(power,powd(R,-3)),(1,-1,-3,0,0))
        self.assertEqual(mul(pg,powd(rho,-1)),mul(U,U))
        self.assertEqual(mul(R,(1,-1,-3,0,0),powd(rho,-1),powd(c,-1)),mul(U,U))


class SavedResultTests(unittest.TestCase):
    def test_15_all_accepted_inputs_unchanged(self):
        inv=json.loads((ROOT/'results/input_inventory.json').read_text())
        self.assertTrue(inv['all_match'])
        for row in inv['files']:
            self.assertEqual(hashlib.sha256((REPO/row['path']).read_bytes()).hexdigest(),row['sha256'])

    def test_16_admission_and_censoring(self):
        a=pd.read_csv(ROOT/'results/input_admission.csv');s=pd.read_csv(ROOT/'results/adiabatic_summary.csv')
        self.assertEqual(len(a),len(s))
        self.assertTrue(s.loc[~s.completed,'global_Tmax_K'].isna().all())
        self.assertTrue((a.loc[~a.completed,'last_time_s']<5/Parameters().f).all())
        for amp,category in [(1.,'conditional_baseline'),(1.15,'conditional_sensitivity'),(1.3,'formal_extreme_scenario')]:
            self.assertEqual(a[(a.model=='KM')&(a.Pa_bar==amp)].iloc[0].next_stage_use,category)

    def test_17_extrema_and_temperature_agreement(self):
        a=pd.read_csv(ROOT/'results/adiabatic_summary.csv')
        self.assertLess(a.refined_vs_interpolant_relative_error.max(),2e-8)
        s=pd.read_csv(ROOT/'results/collapse_temperatures.csv')
        r=s[(s.Pa_bar==1.3)&s.global_deepest]
        self.assertEqual(dict(zip(r.mechanics,r.cycle)),{'KM':1,'RP':2})

    def test_18_saved_gas_identity(self):
        s=pd.read_csv(ROOT/'results/thermal_summary.csv')
        for row in s.to_dict('records'):
            g=Gas(gamma=row['gamma'],T0=row['T0_K'],kg=row['kg_W_m_K'])
            with np.load(ROOT/'results'/(row['run_id']+'.npz')) as z:
                np.testing.assert_allclose(z['pg_Pa']*z['V_m3'],g.n()*RU*z['T_K'],rtol=5e-15)
                self.assertTrue(np.all(z['T_K']>0))
                self.assertEqual(float(z['source_energy_J'][-1]),0.)

    def test_19_coupled_adiabatic_recovery(self):
        s=pd.read_csv(ROOT/'results/adiabatic_recovery.csv')
        self.assertTrue(s.completed.all())
        self.assertLess(s.max_temperature_invariant_relative_error.max(),1e-6)
        self.assertLess(s.relative_Tmax_vs_saved_mechanics.max(),1e-5)

    def test_20_solver_convergence(self):
        s=pd.read_csv(ROOT/'results/solver_convergence.csv')
        self.assertEqual(len(s),39)
        self.assertTrue(s.completed.all())
        # Coarse grid is deliberately diagnostic; production and independent methods must agree.
        checked=s[s.configuration.isin(['production','independent_Radau'])]
        self.assertLess(checked.relative_Tmax_error.max(),2e-4)
        self.assertLess(checked.relative_duration_1000K_error.max(),2e-3)

    def test_21_input_grid_convergence(self):
        s=pd.read_csv(ROOT/'results/input_grid_convergence.csv')
        self.assertEqual(set(s.Pa_bar),{1.,1.15,1.3})
        self.assertEqual(len(s),57)
        fine=s[(s.grid.isin(['export_stride2','full_export_production']))]
        self.assertLess(fine.relative_Tmax_error.max(),2e-4)
        self.assertLess(fine.relative_duration_1000K_error.max(),2e-3)

    def test_22_energy_residual(self):
        s=pd.read_csv(ROOT/'results/thermal_summary.csv')
        self.assertLess(s.energy_residual_scaled.max(),2e-4)
        for name in ['solver_convergence.csv','input_grid_convergence.csv']:
            t=pd.read_csv(ROOT/'results'/name)
            self.assertLess(t.gauss_energy_scaled_residual.max(),2e-4)

    def test_23_series_summary_and_hot_intervals(self):
        s=pd.read_csv(ROOT/'results/thermal_summary.csv',float_precision='round_trip')
        self.assertEqual(len(s),30)
        for row in s.to_dict('records'):
            with np.load(ROOT/'results'/(row['run_id']+'.npz')) as z:
                self.assertLess(abs(float(z['T_K'].max())/row['observed_Tmax_K']-1),2e-12)
                self.assertEqual(float(z['t_s'][-1]),row['last_time_s'])
            meta=json.loads((ROOT/'results'/(row['run_id']+'.json')).read_text())
            for item in meta['hot_intervals']:
                length=sum(b-a for a,b in item['intervals'])
                self.assertAlmostEqual(length,item['duration_s'],places=20)
                self.assertAlmostEqual(row[f"duration_above_{int(item['threshold_K'])}K_s"],length,places=20)

    def test_24_clean_environment_and_source(self):
        d=json.loads((ROOT/'results/clean_environment_verification.json').read_text())
        self.assertTrue(d['passed']);self.assertTrue(d['isolated_environment'])
        self.assertEqual(len(d['cases']),4)


if __name__=='__main__':unittest.main()
