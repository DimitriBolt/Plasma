"""Verify a bounded correction run; saved meshes are input, not remesher targets."""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from .upstream import ROOT
from .final_corrections import check_original_inputs

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('clean_repository',type=Path)
    args=parser.parse_args();other=args.clean_repository/'Mini_Research/06_AuNP_Bubble_Coupling'
    r=ROOT/'results';s=other/'results';files,sources=check_original_inputs()
    arrays=0;npzs=0
    for p in sorted(r.rglob('*.npz')):
        with np.load(p,allow_pickle=False) as a,np.load(other/p.relative_to(ROOT),allow_pickle=False) as b:
            if set(a.files)!=set(b.files):raise RuntimeError('Source/field schema differs: '+p.name)
            for key in a.files:
                # These are preserved archived arrays/re-exported sources, not
                # fresh remesher outputs. Immutable input equality is justified.
                if not np.array_equal(a[key],b[key]):raise RuntimeError('Saved input/export differs '+str(p)+':'+key)
                arrays+=1
            npzs+=1
    claims=pd.read_csv(s/'handoff_numeric_audit.csv')
    fields=pd.read_csv(s/'final_field_verification.csv')
    replay=pd.read_csv(s/'saved_mesh_replay.csv')
    if len(claims)!=59 or not claims.passed.all():raise RuntimeError('Clean numeric claim audit failed')
    if len(fields)!=5 or fields.independently_integrated_water_residual_over_E.max()>1e-10 or fields.wall_power_field_relative.max()>1e-10:
        raise RuntimeError('Clean archived-field quadrature failed')
    cols=[c for c in replay if c.endswith('_maxnorm_relative')]
    worst=float(replay[cols].to_numpy().max())
    if len(replay)!=3 or worst>2e-5:raise RuntimeError('Clean archived-mesh observables failed')
    a=pd.read_csv(r/'gas_feedback_estimate.csv');b=pd.read_csv(s/'gas_feedback_estimate.csv')
    for col in ['estimated_omitted_liquid_return_J','estimated_return_fraction','max_extra_surface_water_K',
                'conservative_prescribed_trace_return_bound_J','already_counted_MR4_incremental_gas_loss_J']:
        if not np.allclose(a[col],b[col],rtol=1e-6,atol=0):raise RuntimeError('Saved-trace analytic estimate differs: '+col)
    jpgs={p.name:p.read_bytes()==(other/'figures'/p.name).read_bytes() for p in (ROOT/'figures').glob('*.jpg')}
    pdf=Path('report/06_AuNP_Bubble_Coupling.pdf')
    text=lambda p:subprocess.check_output(['pdftotext','-layout',str(p),'-']).decode()
    first,second=text(ROOT/pdf),text(other/pdf)
    if first.split()!=second.split():raise RuntimeError('Recompiled scientific text differs')
    tests=(s/'clean_environment_tests.txt').read_text()
    if '31 passed' not in tests:raise RuntimeError('Final clean test suite did not pass')
    record=dict(original_MR6_snapshot='faa2828f554517a0e2e43d9ae42cff37d682fbd0',
        clean_repository=str(args.clean_repository),fresh_locked_environment=True,
        scope='bounded final-correction reproduction; no new sweep or moving-boundary model',
        original_numeric_files_unchanged=files,source_archives_reexported_from_original_schema=sources,
        preserved_or_reexported_NPZ_count=npzs,matching_archived_array_count=arrays,
        scientific_claims_checked=59,independent_stored_field_checks=5,saved_mesh_replay_cases=3,
        mesh_policy='archived nodes/connectivity are input; regenerated mesh identity is never required',
        saved_mesh_observable_maxnorm_tolerance=2e-5,clean_worst_observed_maxnorm_difference=worst,
        saved_trace_feedback_cases=len(b),feedback_reproduction_relative_tolerance=1e-6,
        same_host_JPG_byte_matches=jpgs,extracted_PDF_text_exact=first==second,
        normalized_PDF_text_equal=True,clean_test_output=tests,
        not_repeated=['original parameter grids','original mechanical/gas trajectory production','original convergence sweeps'],
        physical_validation=False)
    (r/'clean_reproduction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(f'Clean corrections: 31 tests, 59 claims, 5 field checks, 3 archived-mesh replays; {npzs} NPZ / {arrays} archived arrays match')

if __name__=='__main__':main()
