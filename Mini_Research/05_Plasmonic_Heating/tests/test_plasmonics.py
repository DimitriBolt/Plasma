"""Physical limits and cross-implementation tests; no bubble integration."""
import hashlib
import json
import unittest
import numpy as np
from scipy.integrate import quad
from src.plasmonics import *


class OpticalTests(unittest.TestCase):
    def test_data_knots_and_bounds(self):
        d=GoldData()
        np.testing.assert_allclose(d.index(d.wavelength_m),d.table[:,1]+1j*d.table[:,2],rtol=1e-14)
        for w in [100e-9,2e-6,np.nan,-1]:
            with self.assertRaises(ValueError): d.index(w)
        with self.assertRaises(ValueError): d.index(500e-9,'cubic')

    def test_radiative_dipole_passivity(self):
        for eps in [-4+2j,-4+.01j,2+0j]:
            c=dipole_cross_sections(100e-9,532e-9,1.33,eps)
            self.assertGreaterEqual(c.absorption_m2,0)
            self.assertAlmostEqual(float((c.absorption_m2+c.scattering_m2)/c.extinction_m2),1,places=12)

    def test_mie_lossless(self):
        c=mie_cross_sections(50e-9,532e-9,1.33,1.5+0j)
        self.assertLess(abs(c.absorption_m2/c.extinction_m2),1e-12)

    def test_host_index_matching(self):
        c=mie_cross_sections(50e-9,532e-9,1.33,1.33+0j)
        self.assertLess(abs(c.extinction_m2),1e-40)

    def test_small_dipole_limit(self):
        c=optical_cross_sections(.1e-9,532e-9)
        d=optical_cross_sections(.1e-9,532e-9,model='dipole')
        self.assertLess(abs(d.absorption_m2/c.absorption_m2-1),2e-5)

    def test_passivity_and_balance_grid(self):
        g=np.load(ROOT/'results/optical_grid.npz')
        self.assertTrue(np.all(g['sigma_abs_m2']>0))
        np.testing.assert_allclose(g['sigma_abs_m2']+g['sigma_sca_m2'],g['sigma_ext_m2'],rtol=2e-14)

    def test_independent_library_and_series(self):
        v=json.loads((ROOT/'results/validation_summary.json').read_text())
        self.assertLess(v['max_miepython_abs_relative'],1e-7)
        self.assertLess(v['max_mie_series_relative'],1e-10)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): optical_cross_sections(-1e-9,532e-9)
        with self.assertRaises(ValueError): mie_cross_sections(1e-9,532e-9,1.33,1-1j)
        with self.assertRaises(ValueError): absorbed_power(-1,1e-15)
        with self.assertRaises(ValueError): ideal_energy_transfer_upper_bound(1,-1)


class EnergyTests(unittest.TestCase):
    def test_pulse_and_dimensions(self):
        self.assertAlmostEqual(float(absorbed_power(1e8,1e-15))/1e-7,1,places=14)
        self.assertEqual(pulse_energy(1,1e-15),1e-15)

    def test_beam_partition_and_thin_limit(self):
        c=optical_cross_sections(20e-9,532e-9)
        for n in [0,1e10,1e17,1e22]:
            e=ensemble_energy(n,20e-9,1e-3,1e-8,1,c)
            self.assertLessEqual(e['absorbed_J'],e['incident_J'])
            self.assertAlmostEqual(float((e['absorbed_J']+e['scattered_J']+e['transmitted_J'])/e['incident_J']),1,places=13)
        e=ensemble_energy(1e10,20e-9,1e-3,1e-8,1,c)
        self.assertLess(abs(e['absorption_fraction']/(1e10*c.absorption_m2*1e-3)-1),1e-7)

    def test_opaque_limit(self):
        c=optical_cross_sections(20e-9,532e-9)
        e=ensemble_energy(1e22,20e-9,1e-3,1e-8,1,c)
        self.assertAlmostEqual(float(e['absorption_fraction']),float(c.absorption_m2/c.extinction_m2),places=13)

    def test_zero_cross_section(self):
        e=ensemble_energy(0,20e-9,1e-3,1e-8,1,CrossSections(0.,0.,0.,.1))
        self.assertEqual(e['transmitted_J'],e['incident_J'])
        self.assertEqual(e['absorbed_J'],0)

    def test_heat_balance_and_positive_states(self):
        h=RadialHeatModel(20e-9,cells=90)
        time=np.r_[0,np.geomspace(1e-12,1e-6,200)]
        r=h.pulse(1,1e-15,10e-9,time)
        self.assertLess(np.max(abs(r['balance_residual_J']))/1e-15,1e-9)
        self.assertGreaterEqual(np.min(r['gold_rise_K']),-1e-12)
        self.assertGreaterEqual(np.min(r['water_rise_K']),-1e-10)
        self.assertTrue(np.all(r['liquid_net_energy_J']<=1.000000001e-15))
        self.assertTrue(np.all(r['liquid_net_energy_J']>=-1e-25))

    def test_adiabatic_short_pulse_limit(self):
        # Mathematical heat-equation limit only; not a physical fs Au model.
        h=RadialHeatModel(20e-9,cells=90)
        v=h.temperatures([1e-16],10.,1e-16)[0,0]
        self.assertLess(abs(v/(1e-15/h.capacity[0])-1),1e-5)

    def test_steady_resistances(self):
        h=RadialHeatModel(20e-9,cells=90)
        v=h.temperatures([1.],1e-7,2.)[0,0]
        expected=1e-7*((1/h.a-1/h.outer)/(4*np.pi*.6)+1/(4*np.pi*h.a**2*1e8))
        self.assertLess(abs(v/expected-1),1e-10)

    def test_heat_refinement_and_independent_balances(self):
        v=json.loads((ROOT/'results/validation_summary.json').read_text())
        self.assertLess(v['max_heat_gold_relative'],.003)
        self.assertLess(v['max_heat_water_surface_relative'],.01)
        self.assertLess(v['interface_quadrature_relative'],1e-8)
        self.assertLess(v['independent_BDF_gold_relative'],1e-7)
        self.assertLess(v['max_depth_quadrature_absolute'],1e-9)

    def test_source_hash_and_upstream_unchanged(self):
        source=json.loads((ROOT/'data/provenance.json').read_text())
        self.assertEqual(hashlib.sha256((ROOT/'data/Johnson.yml').read_bytes()).hexdigest(),source['table_sha256'])
        inventory=json.loads((ROOT/'results/input_inventory.json').read_text())
        for row in inventory['files']:
            self.assertEqual(hashlib.sha256((ROOT.parents[1]/row['path']).read_bytes()).hexdigest(),row['sha256'])

    def test_bounded_maxima_and_gates(self):
        import pandas as pd
        s=json.loads((ROOT/'results/summary.json').read_text())
        single=pd.read_csv(ROOT/'results/bounded_single_particle.csv')
        single=single[single.conservative_geometry]
        self.assertAlmostEqual(single.energy_absorbed_J.max()/s['bounded_single_particle']['energy_absorbed_J'],1,places=12)
        en=pd.read_csv(ROOT/'results/bounded_ensemble.csv')
        self.assertAlmostEqual(en.absorbed_J.max()/s['bounded_ensemble_energy']['absorbed_J'],1,places=12)
        self.assertAlmostEqual(en.absorption_fraction.max()/s['bounded_ensemble_fraction']['absorption_fraction'],1,places=12)
        self.assertTrue(np.all(en.tau_sca<=.1+1e-12))
        self.assertTrue(np.all(en.volume_fraction<=1e-3))
        self.assertTrue(np.all(en.absorbed_J<=en.incident_J))


if __name__=='__main__': unittest.main(verbosity=2)
