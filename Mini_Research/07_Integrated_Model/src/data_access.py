"""Explicit compact-reference versus full-local archive selection."""
import os,json
from .upstream import ROOT

def compact():
    return (ROOT/'results/publication_scope.json').exists() and os.environ.get('MR7_DATA_MODE')!='full'

def data_path(name,suffix='npz'):
    folder=('metadata' if suffix=='json' else 'references') if compact() else 'raw'
    p=ROOT/'results'/folder/(name+'.'+suffix)
    if not p.exists():raise FileNotFoundError(f'{name} is not shipped in this data profile. The full replay regenerates results/raw; set MR7_DATA_MODE=full to use it.')
    return p

def metadata(name):return json.loads(data_path(name,'json').read_text())

def available_dynamic_runs():
    if compact():return [p.stem for p in (ROOT/'results/references').glob('*.npz') if (ROOT/'results/metadata'/(p.stem+'.json')).exists()]
    return [p.stem for p in (ROOT/'results/raw').glob('*.json')]
