"""Publication storage: retain all scalar samples and exact selected liquid fields.

No field value is interpolated, rounded or clipped. All liquid-cell extrema,
pulse/thermal-event neighborhoods and two uniform diagnostic grids are retained.
No-laser baselines keep every accepted solver state for variational replay.
"""
import json,hashlib,re
import numpy as np
from .upstream import ROOT

def digest(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def main():
    rows=[]
    for path in sorted((ROOT/'results/raw').glob('*.npz')):
        if path.name.startswith('static_'):continue
        with np.load(path) as z:
            if 'storage_version' in z:continue
            data={k:z[k] for k in z.files}
        if 'T_l_rise_K' not in data:continue
        t=data['t_s'];n=len(t)
        ids=list(np.searchsorted(t,np.linspace(0,t[-1],1400)))+list(np.searchsorted(t,np.linspace(0,t[-1],1300)))
        for key in ['temperature_stationary_times_s','radius_stationary_times_s']:
            for event in data.get(key,[]):
                ids+=list(np.searchsorted(t,event+np.r_[0.,-np.geomspace(1e-13,1e-6,35),np.geomspace(1e-13,1e-6,35)]))
        meta=json.loads(path.with_suffix('.json').read_text())
        center=meta['summary']['pulse_center_s'];duration=meta['summary']['duration_s']
        pulse=np.r_[np.linspace(center-duration,center+duration,90),center+np.geomspace(1e-13,5e-6,80)]
        ids+=list(np.searchsorted(t,pulse))
        fieldkeys=[k for k in ['T_l_rise_K','dTl_dF','paired_delta_Tl_K'] if k in data]
        for key in fieldkeys:
            ids+=list(np.argmax(data[key],axis=1))+list(np.argmin(data[key],axis=1))
        for key in ['T_Au_rise_K','Q_l_to_g_W','Q_Au_to_l_W','Q_bath_W','paired_delta_T_K','T_g_K','R_m']:
            if key in data:ids += [int(np.argmax(data[key])),int(np.argmin(data[key]))]
        ids=np.unique(np.clip(ids,0,n-1)).astype(int)
        scalar_before={k:digest(v) for k,v in data.items() if np.asarray(v).ndim==1}
        for key in fieldkeys:data[key]=data[key][:,ids]
        # Accepted state duplicates all physical arrays and is unnecessary outside
        # the no-laser reference used for reconstruction/replay.
        keep_reference=path.stem.startswith(('base_','N0_')) or path.stem.endswith('_base')
        if not keep_reference:data.pop('accepted_y',None)
        data['field_indices']=ids;data['field_t_s']=t[ids];data['storage_version']=np.array(2)
        after={k:digest(data[k]) for k in scalar_before}
        if scalar_before!=after:raise AssertionError('Scalar history changed')
        before_bytes=path.stat().st_size;np.savez_compressed(path,**data)
        rows.append(dict(file=path.name,original_times=n,stored_field_times=len(ids),
             before_bytes=before_bytes,after_bytes=path.stat().st_size,scalar_arrays_unchanged=True,
             retained_field_values='exact subset, float64',accepted_state_kept=keep_reference))
    (ROOT/'results/storage_compaction.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('Compacted',len(rows),'files; scalar histories unchanged; exact field snapshots retained.')

if __name__=='__main__':main()
