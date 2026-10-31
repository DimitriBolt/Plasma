"""Write the complete package list and SHA-256 manifest; omit runtime caches."""
import hashlib
from .plasmonics import ROOT


def files():
    suffixes={'.py','.md','.txt','.csv','.json','.npz','.tex','.pdf','.jpg','.ipynb','.yml'}
    return sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
        if p.is_file() and not any(part.startswith('.') or part=='__pycache__' for part in p.relative_to(ROOT).parts)
        and (p.suffix in suffixes or p.name in {'FILE_LIST.txt','SHA256SUMS'}))


def run():
    names=sorted(set(files()+['.gitignore','.python-version','FILE_LIST.txt','SHA256SUMS']))
    required={'data/Johnson.yml','data/provenance.json','src/plasmonics.py',
              'report/05_Plasmonic_Heating.tex','report/05_Plasmonic_Heating.pdf',
              'results/review_handoff.md','requirements-lock.txt','tests/test_plasmonics.py'}
    if not required.issubset(names):
        raise ValueError(f'Missing required package files: {required-set(names)}')
    (ROOT/'FILE_LIST.txt').write_text('\n'.join(names)+'\n')
    lines=[hashlib.sha256((ROOT/name).read_bytes()).hexdigest()+'  '+name for name in names if name!='SHA256SUMS']
    (ROOT/'SHA256SUMS').write_text('\n'.join(lines)+'\n')
    print(f'{len(names)} files; {len(lines)} SHA-256 entries')


if __name__=='__main__':run()
