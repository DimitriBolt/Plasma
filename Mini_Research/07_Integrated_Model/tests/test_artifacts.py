"""Independent scientific acceptance of admitted saved data, not implementation replicas."""
import json
import unittest
import numpy as np
import pandas as pd
from src.upstream import ROOT
from src.data_access import data_path,metadata,available_dynamic_runs,compact
R=ROOT/'results'

@unittest.skipUnless((R/'all_runs.csv').exists(),'Run final scientific curation first')
class ArtifactTests(unittest.TestCase):
    def test_recorded_full_archive_power_replay_coverage(self):
        d=pd.read_csv(R/'constitutive_replay.csv')
        self.assertTrue(d.passed.all());self.assertGreater(len(d),60)
        self.assertEqual(len(d),len(pd.read_csv(R/'all_runs.csv')))

    def test_censoring_never_receives_a_global_temperature(self):
        d=pd.read_csv(R/'all_runs.csv')
        self.assertTrue(d.loc[~d.completed,'global_Tmax_K'].isna().all())
        p=pd.read_csv(R/'stable_pairs.csv')
        self.assertTrue(p.loc[~p.completed_pair,'global_delta_Tmax_K'].isna().all())
        self.assertTrue((d.loc[d.status=='gas_coefficient_upper_boundary','boundary_T_K']>=5999.99).all())

    def test_photon_allowance_and_signed_channel_accounting(self):
        d=pd.read_csv(R/'all_runs.csv')
        self.assertTrue((d.E_absorbed_J<=d.E_incident_J+1e-25).all())
        self.assertTrue((d.photon_allowance_ratio<=1).all())
        self.assertTrue((d.Q_l_to_g_J<0).any());self.assertTrue((d.Q_l_to_g_J>0).any())

    def test_shipped_optical_energies_follow_the_actual_pulse(self):
        checked=0
        for name in available_dynamic_runs():
            meta=metadata(name);s=meta['summary'];p=meta['particle']
            with np.load(data_path(name)) as z:
                if 'E_incident_J' not in z:continue
                checked+=1
                fraction=np.clip((z['t_s']-(s['pulse_center_s']-s['duration_s']/2))/s['duration_s'],0,1)
                full=s['fluence_J_m2']*p['count']*p['sigma_abs_m2']*p['absorption_scale']
                beam=s['fluence_J_m2']*np.pi*p['beam_radius_m']**2
                self.assertLess(np.max(abs(z['E_absorbed_J']-full*fraction)),max(full*3e-8,1e-25),name)
                self.assertLess(np.max(abs(z['E_incident_J']-beam*fraction)),max(beam*3e-8,1e-25),name)
        self.assertGreaterEqual(checked,8 if compact() else 60)

    def test_complete_native_weak_pair_agrees_with_stable_formulation(self):
        d=pd.read_csv(R/'independent_Radau.csv').iloc[0]
        s=pd.read_csv(R/'stable_pairs.csv')
        r=s[(s.Pa_bar==1)&(s.level==3)&(s.order==3)&(s.fluence_J_m2==10)].iloc[0]
        self.assertTrue(d.completed_pair);self.assertGreater(d.global_delta_Tmax_K,0)
        self.assertLess(abs(d.global_delta_Tmax_K/r.global_delta_Tmax_K-1),.005)

    def test_nonlinear_increment_and_quadrature_control(self):
        d=pd.read_csv(R/'stable_pairs.csv')
        g=d[(d.Pa_bar==1)&(d.level==3)&(d.fluence_J_m2==10)]
        self.assertLess(np.ptp(g.global_delta_Tmax_K),1e-10)
        a=d[(d.Pa_bar==1)&(d.level==3)&(d.order==3)]
        slope=a.global_delta_Tmax_K/a.fluence_J_m2
        self.assertLess(slope.max()/slope.min()-1,2e-4)

    def test_incremental_energy_is_resolved_on_the_absorbed_scale(self):
        d=pd.read_csv(R/'paired_energy.csv')
        self.assertLess(d.total_over_absorbed.max(),1e-6)
        self.assertLess(d.water_over_absorbed.max(),1e-8)
        self.assertLess(d.Au_over_absorbed.max(),1e-8)

    def test_weak_early_sign_and_unresolved_late_tail_are_distinct(self):
        d=pd.read_csv(R/'strong_collapse_resolution.csv')
        early=d[d.offset_s< -1e-8]
        self.assertTrue((early.cells96_slope>0).all());self.assertLess(early.relative_mesh_change.max(),.02)
        late=d[np.isclose(d.offset_s,0,atol=1e-15)].iloc[0]
        self.assertEqual(late.numerical_resolution,'unresolved_late_tail')
        self.assertGreater(late.relative_mesh_change,.1)

    def test_no_source_precedes_its_pulse(self):
        d=pd.read_csv(R/'timing_sensitivity.csv')
        late=d[(d.center_s-5e-9>d.first_cycle_t_Tmax_s)&d.first_cycle_complete]
        self.assertGreater(len(late),0)
        self.assertLess(abs(late.first_cycle_dT_dF_K_m2_J).max(),1e-15)
        self.assertTrue(d.loc[~d.completed,'global_dTmax_dF_K_m2_J'].isna().all())

    def test_clamped_MR4_and_transient_conduction_controls(self):
        d=pd.read_csv(R/'MR4_regression_full.csv')
        self.assertLess(abs(d.peak_relative_difference).max(),1e-6)
        self.assertTrue((d.MR7_completed==d.MR4_completed).all())
        a=pd.read_csv(R/'analytic_point_control.csv');a=a[(a.cells==320)&(a.t_s>0)]
        self.assertLess(abs(a.relative_error).max(),.005)
        self.assertLess(abs(a.balance_over_initial).max(),1e-8)

    def test_local_Au_error_is_retained(self):
        d=pd.read_csv(R/'static_MR6_comparison.csv');d=d[d.cells==320]
        self.assertTrue((d.radial_gold_peak_K/d.MR6_gold_peak_K<.4).all())
        self.assertTrue((d.balance_over_absorbed<1e-7).all())
        self.assertEqual(len(pd.read_csv(R/'static_domain_shape.csv')),3)

if __name__=='__main__':unittest.main()
