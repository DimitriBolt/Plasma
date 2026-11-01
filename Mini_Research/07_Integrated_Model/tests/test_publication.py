"""The public data subset has explicit coverage and actual readable references."""
import json,unittest
from src.upstream import ROOT
from src.data_access import compact,metadata
from src.reference_verify import check_references

@unittest.skipUnless(compact(),'Compact publication profile only')
class PublicationTests(unittest.TestCase):
    def test_lossless_array_inventory_and_current_physical_powers(self):
        rows=check_references();self.assertEqual(len(rows),4)
        self.assertTrue(all(r['passed'] for r in rows))

    def test_three_amplitudes_and_censoring_are_preserved(self):
        rows=check_references();self.assertEqual({r['Pa_bar'] for r in rows},{1.,1.15,1.3})
        high=[r for r in rows if r['Pa_bar']==1.3]
        self.assertTrue(all(not r['completed'] and r['status']=='gas_coefficient_upper_boundary' for r in high))

    def test_every_admitted_run_keeps_its_metadata(self):
        import pandas as pd
        d=pd.read_csv(ROOT/'results/all_runs.csv');self.assertEqual(len(d),115)
        for name in d.run_id:self.assertEqual(metadata(name)['summary']['run_id'],name)

if __name__=='__main__':unittest.main()
