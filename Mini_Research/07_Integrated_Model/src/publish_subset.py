"""Lossless, explicitly selected data for the compact public package.

Read the complete local MR7 archive; never modify it. Every retained value is
copied without interpolation, rounding or refitting. All run metadata and summary
tables are public, while high-dimensional sweep outputs stay local.
"""
import argparse,hashlib,json,shutil
from pathlib import Path
import numpy as np
from .upstream import ROOT

def digest(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True,type=Path);a=p.parse_args()
    src=a.source.resolve()/'results/raw';out=ROOT/'results/references';out.mkdir(parents=True,exist_ok=True)
    md=ROOT/'results/metadata';md.mkdir(parents=True,exist_ok=True)
    for f in src.glob('*.json'):shutil.copy2(f,md/f.name)
    baseline=['base_1.00','base_1.15','base_1.30','extra_high_1.30_BDF_base']
    pairs=['stable_base_1.00_F10_L3_q3','stable_base_1.15_F10_L3_q3','stable_extra_high_1.30_BDF_base_F10_L3_q3']
    particle=['sens_passive_near_base','sens_passive_far_base','sens_fixed_base']
    dry=['dry_air_base','dry_argon_base'];names=baseline+pairs+particle+dry+['laser_F10_1.00','static_reference_320']
    record=[]
    for name in names:
        f=src/(name+'.npz')
        with np.load(f) as z:
            fields={};selections={}
            if name in baseline or name in pairs or name=='laser_F10_1.00':
                fields={k:z[k] for k in z.files if z[k].ndim<=1}
                if name in baseline:fields['accepted_y']=z['accepted_y']
                if name in pairs:fields['paired_delta_ledgers_J']=z['paired_delta_ledgers_J']
            elif name in particle:fields={k:z[k] for k in ['t_s','r_np_m','gap_m']}
            elif name in dry:fields={k:z[k] for k in ['t_s','R_m','U_m_s','T_g_K','p_g_Pa','peak_candidates_s','peak_candidates_K']}
            else:fields={k:z[k] for k in ['time_s','gas_energy_J','gold_rise_K']}
            if name=='laser_F10_1.00':
                pos=np.unique(np.minimum(np.searchsorted(z['field_t_s'],np.array([10,14,16,18,33])*1e-6),len(z['field_t_s'])-1))
                fields.update(T_l_rise_K=z['T_l_rise_K'][:,pos],field_t_s=z['field_t_s'][pos],field_indices=z['field_indices'][pos])
                selections={k:dict(axis=1 if k=='T_l_rise_K' else 0,positions=pos.tolist()) for k in ['T_l_rise_K','field_t_s','field_indices']}
            entries=[]
            for k,v in fields.items():
                selection=selections.get(k)
                expected=np.take(z[k],selection['positions'],axis=selection['axis']) if selection else z[k]
                assert np.array_equal(v,expected,equal_nan=True)
                entries.append(dict(key=k,shape=list(v.shape),dtype=str(v.dtype),sha256=digest(v),source_array_sha256=digest(z[k]),selection=selection))
            np.savez_compressed(out/f.name,**fields)
        record.append(dict(run_id=name,source_raw_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
            role='full accepted baseline states' if name in baseline else ('primary nonlinear scalar response' if name in pairs else 'figure or formal-extreme trace'),
            arrays=entries,bytes=(out/f.name).stat().st_size))
    scope=dict(profile='compact-publication-v1',local_full_snapshot='c467678d2f48db51ceafe09ea4cb5040b6ecc293',
        local_full_snapshot_is_not_published=True,metadata_runs=len(list(md.glob('*.json'))),
        selected_runs=names,full_accepted_state_runs=baseline,reference_bytes=sum(x['bytes'] for x in record),
        exact_values=True,interpolation_or_rounding=False,
        public_scope='All source, report, figures, parameter/summary/convergence/applicability tables, all metadata, four accepted baseline states at three amplitudes and selected exact scalar/profiling traces.',
        omitted='All other per-sweep NPZ fields, nonessential solver matrices, and local prototype/production raw files. Full raw archives remain local; expensive scripts can regenerate them.',
        verification_scope='Historical full-archive validation is retained separately. Compact tests and reference checks explicitly cover only shipped arrays.')
    (ROOT/'results/publication_scope.json').write_text(json.dumps(scope,indent=2)+'\n')
    (ROOT/'results/public_data_inventory.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in scope.items() if k not in ['selected_runs']},indent=2))

if __name__=='__main__':main()
