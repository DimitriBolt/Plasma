"""Bounded fresh-process reproduction, usable from a clean publication archive."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from .upstream import ROOT,REPO
from .run_experiments import run_case,run_key,write_json

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=ROOT/'results');args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    saved=pd.read_csv(ROOT/'results/summary.csv',float_precision='round_trip');rows=[]
    cases=[(1.3,'air','coupled','M1_variable',1.),(1.3,'argon','one_way','M1_variable',1.),
           (1.15,'argon','coupled','M2_sat',1.),(1.3,'argon','coupled','M2_sat',1.)]
    for task in cases:
        key=run_key(*task);ref=saved[saved.run_id==key].iloc[0]
        m,z,sol=run_case(*task)
        diff=m['observed_Tmax_K']/ref.observed_Tmax_K-1
        assert abs(diff)<1e-8 and bool(ref.completed)==m['completed'],key
        rows.append(dict(run_id=key,published_Tmax_K=ref.observed_Tmax_K,fresh_Tmax_K=m['observed_Tmax_K'],
                         relative_difference=diff,fresh_t_Tmax_s=m['t_Tmax_s'],published_t_Tmax_s=ref.t_Tmax_s,
                         completed=m['completed'],energy_residual_scaled=m['energy_residual_scaled']))
        print(key,diff,m['completed'],flush=True)
    pd.DataFrame(rows).to_csv(args.output_dir/'clean_reproduction.csv',index=False)
    inventory=json.loads((ROOT/'results/input_inventory.json').read_text())
    for item in inventory['files']:
        assert hashlib.sha256((REPO/item['path']).read_bytes()).hexdigest()==item['sha256'],item['path']
    hashes={rel:hashlib.sha256((ROOT/rel).read_bytes()).hexdigest() for rel in
            ['results/summary.csv','results/species_coefficients.json','report/04_Temperature_Barrier.pdf']}
    write_json(args.output_dir/'clean_snapshot_verification.json',dict(
        isolated_checkout=str(REPO),python=platform.python_version(),numpy=np.__version__,
        bounded_reproduction_cases=len(rows),max_relative_peak_difference=max(abs(r['relative_difference']) for r in rows),
        all_statuses_match=True,unchanged_upstream_files=len(inventory['files']),scientific_artifact_hashes=hashes,
        scope='Clean source/data/PDF copy, before adding this verification journal; no old results regenerated'))

if __name__=='__main__': main()
