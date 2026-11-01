"""Write or verify the exact public file inventory and SHA256 manifest."""
import argparse,hashlib,json
from pathlib import Path
from .upstream import ROOT

def admitted(p):
    q=p.relative_to(ROOT)
    return not any(x in {'.venv','tmp','__pycache__','local_archive'} for x in q.parts) and q.parts[:2]!=('results','raw') and not (
        p.suffix in {'.pyc','.aux','.log','.out','.fls','.fdb_latexmk','.gz'} or p.name=='SHA256SUMS')

def main():
    a=argparse.ArgumentParser();a.add_argument('--write',action='store_true');args=a.parse_args()
    if args.write:
        names=sorted({str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and admitted(p)}|{'FILE_LIST.txt','SHA256SUMS'})
        (ROOT/'FILE_LIST.txt').write_text('\n'.join(names)+'\n')
        lines=[hashlib.sha256((ROOT/n).read_bytes()).hexdigest()+'  '+n for n in names if n!='SHA256SUMS']
        (ROOT/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
    lines=(ROOT/'SHA256SUMS').read_text().splitlines();bad=[]
    for line in lines:
        digest,name=line.split('  ',1);p=ROOT/name
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:bad.append(name)
    listed=set((ROOT/'FILE_LIST.txt').read_text().splitlines())
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and admitted(p)}|{'SHA256SUMS'}
    if listed!=actual:raise AssertionError('Inventory mismatch: '+str(listed^actual))
    if bad:raise AssertionError('Checksum mismatch: '+str(bad))
    print(json.dumps(dict(files=len(listed),checksums=len(lines),all_match=True,
        bytes=sum((ROOT/n).stat().st_size for n in listed))))

if __name__=='__main__':main()
