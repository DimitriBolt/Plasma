"""Compare a freshly recalculated exported repository with this production tree."""
import argparse,hashlib,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from .upstream import ROOT

def main():
    parser=argparse.ArgumentParser();parser.add_argument('clean_root',type=Path);args=parser.parse_args()
    other=args.clean_root/'Mini_Research/06_AuNP_Bubble_Coupling'
    csvs=[];arrays=[];figures=[]
    for p in sorted((ROOT/'results').glob('*.csv')):
        q=other/'results'/p.name
        if not q.exists():raise RuntimeError('Missing clean result '+str(q))
        a=pd.read_csv(p,float_precision='round_trip');b=pd.read_csv(q,float_precision='round_trip')
        pd.testing.assert_frame_equal(a,b,check_exact=True)
        csvs.append(p.name)
    for p in sorted((ROOT/'results').rglob('*.npz')):
        q=other/p.relative_to(ROOT)
        with np.load(p,allow_pickle=False) as a,np.load(q,allow_pickle=False) as b:
            if set(a.files)!=set(b.files):raise RuntimeError('NPZ fields differ '+str(p))
            for k in a.files:
                if a[k].dtype.kind in 'fc':same=np.array_equal(a[k],b[k],equal_nan=True)
                else:same=np.array_equal(a[k],b[k])
                if not same:raise RuntimeError('NPZ array differs: '+str(p)+' '+k)
                arrays.append(str(p.relative_to(ROOT))+':'+k)
    for p in sorted((ROOT/'figures').glob('*.jpg')):
        q=other/'figures'/p.name
        if p.read_bytes()!=q.read_bytes():raise RuntimeError('Figure bytes differ '+p.name)
        figures.append(p.name)
    pdf=Path('report/06_AuNP_Bubble_Coupling.pdf')
    txt=lambda p:subprocess.check_output(['pdftotext','-layout',str(p),'-'])
    if txt(ROOT/pdf)!=txt(other/pdf):raise RuntimeError('Recompiled PDF text differs')
    checks={}
    for name in ['parameters.json','design_extrema.json','validation_summary.json','optical_contract.json','local_flux_metrics.json']:
        if json.loads((ROOT/'results'/name).read_text())!=json.loads((other/'results'/name).read_text()):
            raise RuntimeError('Scientific JSON differs: '+name)
        checks[name]=True
    result=dict(clean_repository=str(args.clean_root),fresh_locked_environment=True,
                exact_csv_count=len(csvs),exact_csv_files=csvs,exact_NPZ_array_count=len(arrays),
                exact_JPG_count=len(figures),exact_JPG_files=figures,scientific_json_matches=checks,
                extracted_PDF_text_exact=True,
                tests_output=(other/'results/clean_environment_tests.txt').read_text(),
                scope='same-machine fresh-environment full recomputation; not independent physical validation')
    (ROOT/'results/clean_reproduction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(len(csvs),'CSV,',len(arrays),'NPZ arrays,',len(figures),'JPG, scientific JSON and PDF text reproduce exactly')

if __name__=='__main__':main()
