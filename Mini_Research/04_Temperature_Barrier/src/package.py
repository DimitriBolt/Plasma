"""Deterministic publication inventory; never enumerates user files outside MR4."""
import hashlib
from .upstream import ROOT

def files():
    excluded={'.venv','__pycache__','.pytest_cache'}
    suffixes={'.aux','.log','.out','.fls','.fdb_latexmk','.gz'}
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and
                  not excluded.intersection(p.relative_to(ROOT).parts) and
                  p.suffix not in suffixes and p.name not in {'SHA256SUMS','FILE_LIST.txt'})

def main():
    paths=files();names=sorted([str(p.relative_to(ROOT)) for p in paths]+['FILE_LIST.txt','SHA256SUMS'])
    (ROOT/'FILE_LIST.txt').write_text('\n'.join(names)+'\n')
    paths=sorted(paths+[ROOT/'FILE_LIST.txt'])
    (ROOT/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n' for p in paths))
    print(f'Publication files: {len(names)}; hashes: {len(paths)}')

if __name__=='__main__': main()
