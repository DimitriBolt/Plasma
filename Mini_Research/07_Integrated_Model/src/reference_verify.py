"""Current compact-data integrity and constitutive checks, distinct from history."""
from dataclasses import fields
import hashlib,json
import numpy as np
import pandas as pd
from .upstream import ROOT,T0
from .model import Model,Config
from .particle import Particle
from .data_access import metadata,data_path

def check_references():
    R=ROOT/'results';scope=json.loads((R/'publication_scope.json').read_text())
    inventory=json.loads((R/'public_data_inventory.json').read_text());rows=[]
    assert set(scope['selected_runs'])=={x['run_id'] for x in inventory}
    for entry in inventory:
        with np.load(data_path(entry['run_id'])) as z:
            assert set(z.files)=={x['key'] for x in entry['arrays']}
            for item in entry['arrays']:
                value=z[item['key']];digest=hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
                assert digest==item['sha256'] and list(value.shape)==item['shape']
                if not item['selection']:assert digest==item['source_array_sha256']
    for name in scope['full_accepted_state_runs']:
        md=metadata(name);m=Model(Config(**md['config']),Particle(**{k:v for k,v in md['particle'].items() if k in {f.name for f in fields(Particle)}}))
        with np.load(data_path(name)) as z:
            t=z['accepted_t_s'];y=z['accepted_y'];assert y.shape==(m.dim,len(t))
            ids=np.unique(np.linspace(0,len(t)-1,min(70,len(t))).astype(int));errs=[]
            for j in ids:
                i=np.searchsorted(z['t_s'],t[j]);assert z['t_s'][i]==t[j]
                vals=np.array(m.powers(t[j],y[:,j])[:4])
                archived=np.array([z[k][i] for k in ['Q_l_to_g_W','Q_Au_to_l_W','Q_bath_W','work_by_gas_W']])
                errs.append(abs(vals-archived));assert abs(y[2,j]*T0-z['T_g_K'][i])<1e-9
            err=np.max(errs,axis=0);scale=np.array([max(abs(z[k])) for k in ['Q_l_to_g_W','Q_Au_to_l_W','Q_bath_W','work_by_gas_W']])
            passed=bool(np.all((err/np.maximum(scale,1e-15)<1e-8)|(err<1e-17)))
            assert passed,name
            rows.append(dict(run_id=name,Pa_bar=m.cfg.amplitude_bar,sampled_accepted_states=len(ids),
                 max_power_error_W=float(err.max()),passed=passed,completed=md['summary']['completed'],
                 status=md['summary']['status'],scope='current shipped full-state references only'))
    return rows

def main():
    rows=check_references();pd.DataFrame(rows).to_csv(ROOT/'results/reference_checks.csv',index=False,float_format='%.17g')
    print(json.dumps(dict(shipped_arrays_verified=True,current_baseline_power_replays=rows),indent=2))

if __name__=='__main__':main()
