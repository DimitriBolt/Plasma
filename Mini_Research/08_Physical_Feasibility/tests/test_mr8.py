import unittest,json
import numpy as np
import pandas as pd
from scipy.integrate import quad
from src.mechanics import Tait,gilmore_rhs,implicit_residual,mr2,solve
from src.upstream import ROOT,thermo
from src.run import inventory

class MathematicalChecks(unittest.TestCase):
    def test_reference_liquid(self):
        p=mr2.Parameters();e=Tait();rl,ri,c,h=e.quantities(p.p0,p.p0)
        self.assertAlmostEqual(rl/p.rho,1.,places=14)
        self.assertAlmostEqual(c/p.c,1.,places=14)
        self.assertEqual(h,0.)

    def test_enthalpy_integral_and_derivative(self):
        e=Tait()
        for pi,pl in [(-28675.,2e9),(101325.,101326.),(2e7,1e6)]:
            h=e.quantities(pl,pi)[-1]
            ref=quad(lambda q:1/e.density(q),pi,pl,epsabs=1e-10)[0]
            self.assertLess(abs(h-ref)/max(abs(ref),1e-10),2e-10)
            step=10.;dh=(e.quantities(pl+step,pi)[-1]-e.quantities(pl-step,pi)[-1])/(2*step)
            self.assertLess(abs(dh*e.density(pl)-1),1e-7)

    def test_outside_tait_branch_rejected(self):
        with self.assertRaises(ValueError):Tait().density(-Tait().B)

    def test_full_viscous_implicit_equation(self):
        for R,U,t in [(6e-6,0.,0.),(1e-6,-300.,18e-6),(3e-6,40.,20e-6),(30e-6,-1.,10e-6)]:
            residual,scale=implicit_residual(t,[R,U],1.3e5)
            self.assertLess(abs(residual)/scale,2e-14)

    def test_incompressible_limit(self):
        p=mr2.Parameters();y=[3e-6,-30.];t=17e-6;pa=1.15e5
        expected=mr2.rp.rhs(t,y,pa,p)[1]
        errors=[abs(gilmore_rhs(t,y,pa,Tait(sound_multiplier=m))[1]/expected-1) for m in [1e3,1e5,1e7]]
        self.assertLess(errors[-1],1e-7)
        self.assertTrue(errors[2]<errors[1]<errors[0])

    def test_unforced_stationary_state(self):
        p=mr2.Parameters()
        self.assertEqual(gilmore_rhs(0,[p.R0,0],0)[1],0.)
        row,sol,t=solve(0.,'Gilmore',rtol=1e-10,steps=200)
        self.assertLess(np.max(abs(sol.y[0]-p.R0))/p.R0,1e-12)

    def test_temperature_domain_and_caloric_identity(self):
        m=thermo.Mixture()
        with self.assertRaises(ValueError):m.du(6001.)
        for T in [550.,2798.,5900.]:
            h=.001
            derivative=(m.du(T+h)-m.du(T-h))/(2*h)
            self.assertLess(abs(derivative/m.capacity(T)-1),1e-8)

class SavedEvidenceChecks(unittest.TestCase):
    def test_mechanical_convergence(self):
        d=pd.read_csv(ROOT/'results/convergence.csv')
        self.assertEqual(len(d),18)
        self.assertLess(d.relative_Rmin.max(),2e-7)
        self.assertLess(d.relative_speed.max(),2e-6)

    def test_common_KM_event_matches_accepted_data(self):
        d=pd.read_csv(ROOT/'results/mechanics.csv')
        ref=d[(d.model=='KM')&(d.Pa_bar==1.3)&(d.control=='strict')].iloc[0]
        self.assertLess(abs(ref.Rmin_m/0.550808707142226e-6-1),1e-8)
        self.assertLess(abs(ref.t_Rmin_s-18.34106426519485e-6),1e-13)

    def test_energy_identities_and_signed_flux(self):
        d=pd.read_csv(ROOT/'results/energy_audit.csv')
        self.assertLess(d.total_scaled_residual.max(),1e-5)
        self.assertTrue((d.signed_gas_power_min_W<0).all())
        self.assertTrue((d.signed_gas_power_max_W>0).all())
        p=pd.read_csv(ROOT/'results/paired_energy_audit.csv')
        self.assertLess(p.independent_gas_increment_residual_over_absorption.max(),1e-6)

    def test_fixed_volume_energy_not_constant_capacity_bound(self):
        d=pd.read_csv(ROOT/'results/ideal_caloric_allowance.csv')
        for _,r in d.iterrows():
            m=thermo.Mixture(gas=r.gas,vapor_factor=r.vapor_factor)
            delta=m.du(r.baseline_T_K+r.ideal_all_absorbed_at_fixed_volume_delta_T_K)-m.du(r.baseline_T_K)
            self.assertLess(abs(delta/r.E_absorbed_J-1),1e-6)

    def test_preserve_all_preexisting_tracked_files(self):
        self.assertEqual(inventory()['changed'],[])

if __name__=='__main__':unittest.main()
