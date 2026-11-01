"""Independent physical invariants, limiting cases, terminal events and artifacts."""
import unittest
import hashlib
import json
from dataclasses import replace
import numpy as np
from scipy.integrate import solve_ivp
from src.upstream import ROOT,REPO,thermo,mr4,T0,Parameters
from src.model import Config,Model,TS
from src.particle import Particle
from src.transport import MaterialGrid,center_radius
from src.benchmarks import stationary

class PhysicsTests(unittest.TestCase):
    def test_quiescent_equilibrium_without_drive(self):
        m=Model(Config(amplitude_bar=0,cells=8,cycles=.1),Particle(count=0))
        s=m.solve();self.assertTrue(s.completed)
        np.testing.assert_allclose(s.sol(s.end),m.initial(),rtol=1e-11,atol=1e-12)

    def test_adiabatic_argon_invariant_with_zero_gas_heat(self):
        m=Model(Config(amplitude_bar=.2,gas='argon',vapor_factor=0,eta=0,cells=8,cycles=.2),Particle(count=0))
        s=m.solve();y=s.sol(s.t)
        np.testing.assert_allclose(y[2]*y[0]**2,1,rtol=2e-7,atol=1e-9)

    def test_insulated_liquid_wall_temperature_limit(self):
        from src.diagnostics import evaluate
        m=Model(Config(eta=0,cells=8,cycles=.03),Particle(count=0))
        s=m.solve();r,z=evaluate(s,do_quad=False)
        np.testing.assert_allclose(z['T_wall_K'],T0+z['T_l_rise_K'][0],rtol=1e-14)
        np.testing.assert_array_equal(z['Q_l_to_g_W'],0.)
    def test_material_invariant_against_particle_ode(self):
        R0=6e-6;U=-2.;d=1e-6;a=20e-9;t=np.linspace(0,1e-6,100)
        sol=solve_ivp(lambda t,y:[(R0+U*t)**2*U/y[0]**2],(0,t[-1]),[R0+a+d],t_eval=t,rtol=1e-12,atol=1e-18)
        exact=center_radius(R0+U*t,a,d,'passive')
        np.testing.assert_allclose(sol.y[0],exact,rtol=1e-10,atol=1e-17)
        self.assertLess(np.ptp(exact**3-(R0+U*t)**3),1e-29)

    def test_material_cell_capacities_do_not_move(self):
        g=MaterialGrid();c=g.capacity.copy()
        for R in (1e-6,6e-6,30e-6):
            _,r=g.radii(R);V=4*np.pi*np.diff(r**3)/3
            np.testing.assert_allclose(V*g.rho_cp,c,rtol=1e-10)
        self.assertTrue(np.all(c>0))

    def test_exact_spherical_steady_flux(self):
        g=MaterialGrid();R=4e-6;r,_=g.radii(R);_,K,_=g.conductances(R)
        T=17/r+23
        flux=K*(T[:-1]-T[1:])
        np.testing.assert_allclose(flux,4*np.pi*g.k*17,rtol=1e-12)

    def test_source_weights_exact_normalization_and_motion(self):
        g=MaterialGrid(cells=80)
        for R in (2e-6,6e-6,20e-6):
            for mode in ('passive','capture','fixed'):
                r=center_radius(R,20e-9,1e-6,mode)
                if r-20e-9<R:continue
                w=g.surface_weights(R,r,20e-9)
                self.assertAlmostEqual(w.sum(),1,places=12)
                self.assertGreaterEqual(w.min(),0)

    def test_liquid_internal_flux_cancellation(self):
        m=Model(Config(cells=24),Particle())
        R=4e-6;tl=np.linspace(30,-2,24);r=center_radius(R,20e-9,100e-9,'capture')
        w=m.grid.surface_weights(R,r,20e-9)
        rate,qb=m.grid.water_rate(R,tl,-2e-6,3e-6,w)
        self.assertAlmostEqual(m.grid.capacity@rate,3e-6+2e-6-qb,delta=1e-17)

    def test_no_geometric_heat_creation(self):
        g=MaterialGrid();T=np.zeros(g.cells)
        for R in (1e-6,6e-6,30e-6):
            rate,bath=g.water_rate(R,T,0.,0.,np.zeros(g.cells))
            np.testing.assert_array_equal(rate,0.);self.assertEqual(bath,0.)

    def test_flux_sign_and_interface_reconstruction(self):
        m=Model(Config(cells=16),Particle(count=0));y=m.initial()
        y[2]=2;self.assertLess(m.powers(0,y)[0],0)
        y[2]=.8;self.assertGreater(m.powers(0,y)[0],0)
        z=Model(replace(m.cfg,eta=0),m.particle)
        self.assertEqual(z.powers(0,y)[0],0)

    def test_clamped_MR4_equation_regression(self):
        m=Model(Config(clamp_water=True,cells=8),Particle(count=0));p=Parameters()
        for R,U,T in [(p.R0,0,T0),(2e-6,-100,1500),(12e-6,20,230)]:
            y=m.initial();y[:3]=[R/p.R0,U/p.c,T/T0];t=14e-6
            r=m.rhs(t/TS,y)/TS
            self.assertAlmostEqual(r[2]*T0,m.mix.temperature_rate(R,U,T,1),delta=1e-5)
            self.assertAlmostEqual(r[1]*p.c,mr4.acceleration(t,R,U,T,1,m.mix,1),delta=1e-2)

    def test_complex_direction_against_independent_finite_difference(self):
        m=Model(Config(cells=16),Particle());y=m.initial();y[1]=-.1;y[2]=2
        y[3:m.au]=np.linspace(.03,.001,16);y[m.au]=.02
        direction=np.random.default_rng(17).normal(size=m.dim)*.003
        cs=m.rhs(12,y.astype(complex)+1j*1e-20*direction).imag/1e-20
        eps=1e-5
        fd=(m.rhs(12,y+eps*direction)-m.rhs(12,y-eps*direction))/(2*eps)
        np.testing.assert_allclose(cs,fd,rtol=3e-5,atol=2e-7)

    def test_photon_budget_and_units(self):
        p=Particle(count=8);meta=p.metadata()
        self.assertAlmostEqual(meta['sigma_abs_m2']/3.4337324474063743e-15,1,places=12)
        self.assertLess(meta['N_sigma_ext_over_beam_area'],1)
        self.assertEqual(Model(particle=Particle(count=0)).abs_per_fluence,0)
        self.assertEqual(Model(particle=Particle(absorption_scale=0)).abs_per_fluence,0)

    def test_coefficients_are_not_extrapolated(self):
        for v,limit in [(1.,6000.),(0.,20000.)]:
            mix=thermo.Mixture(vapor_factor=v)
            self.assertEqual(mix.Tmax,limit)
            with self.assertRaises(ValueError):mix.capacity(limit+1)
        self.assertTrue(np.isinf(thermo.Mixture(gas='argon',vapor_factor=0).Tmax))

    def test_zero_G_retains_energy_in_Au(self):
        m=Model(Config(cells=12,cycles=.03),Particle(G_W_m2K=0))
        b=m.solve();s=m.solve(fluence=.1,center_s=.5e-6)
        np.testing.assert_allclose(s.sol(s.end)[:3],b.sol(b.end)[:3],rtol=2e-8,atol=1e-10)
        y=s.sol(s.end)
        self.assertAlmostEqual(y[m.au]*T0*m.particle.capacity/m.abs_per_fluence,.1,delta=1e-8)

    def test_zero_N_and_zero_absorption_exact_source(self):
        for particle in (Particle(count=0),Particle(absorption_scale=0)):
            m=Model(Config(cells=12,cycles=.03),particle)
            s=m.solve(fluence=1,center_s=.5e-6)
            self.assertEqual(s.sol(s.end)[m.ledger+4],0)

    def test_contact_terminates_before_intersection(self):
        m=Model(Config(cells=12,cycles=.2),Particle(motion='fixed',gap_m=1e-9))
        s=m.solve();self.assertEqual(s.status,'particle_contact');self.assertFalse(s.completed)
        R=m.state(s.end,s.sol(s.end))[0]
        self.assertAlmostEqual(center_radius(R,20e-9,1e-9,'fixed')-R-20e-9,0,delta=1e-17)
        with self.assertRaises(ValueError):s.sol(s.end+1e-8)

    def test_actual_gas_boundary_censors(self):
        m=Model(Config(cells=12,cycles=1),Particle(count=0))
        object.__setattr__(m.mix,'Tmax',294.)
        s=m.solve();self.assertFalse(s.completed);self.assertEqual(s.status,'gas_coefficient_upper_boundary')
        self.assertAlmostEqual(s.sol(s.end)[2]*T0,294.,places=7)
        from src.diagnostics import extrema
        self.assertTrue(np.isnan(extrema(s)['global_Tmax_K']))
        self.assertEqual(extrema(s)['selected_event'],'observed_interval_endpoint_reference')

    def test_static_conservation_independent_energy(self):
        v=stationary(cells=60)
        E=3.4337324474063743e-15
        self.assertLess(np.max(abs(v['balance_residual_J']))/E,3e-8)
        self.assertLess(abs(np.trapezoid(v['q_gas_W'],v['time_s'])/v['gas_energy_J'][-1]-1),.002)

    def test_full_source_inventory_preserved(self):
        inv=json.loads((ROOT/'results/input_inventory.json').read_text())
        for r in inv['files']:
            self.assertEqual(hashlib.sha256((REPO/r['path']).read_bytes()).hexdigest(),r['sha256'],r['path'])

if __name__=='__main__':unittest.main()
