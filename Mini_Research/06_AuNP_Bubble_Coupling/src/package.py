"""Explicit complete artifact inventory, excluding caches and TeX intermediates."""
from pathlib import Path
import hashlib
from .upstream import ROOT

def files():
    allowed={'.py','.md','.csv','.json','.npz','.jpg','.pdf','.tex','.txt','.ipynb'}
    special={'.gitattributes','.gitignore','.python-version'}
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and not any(
        x in p.parts for x in ['__pycache__','.pytest_cache','.venv']) and (p.suffix in allowed or p.name in special)
        and p.name not in ('FILE_LIST.txt','SHA256SUMS'))

def main():
    all_files=files(); names=[str(p.relative_to(ROOT)) for p in all_files]+['FILE_LIST.txt','SHA256SUMS']
    required=['src/coupling.py','report/06_AuNP_Bubble_Coupling.pdf','results/reference.npz',
              'results/review_handoff.md','results/clean_reproduction.json','results/input_inventory.json']
    if any(x not in names for x in required):raise RuntimeError('Missing required package artifact')
    (ROOT/'FILE_LIST.txt').write_text('\n'.join(sorted(names))+'\n')
    all_files.append(ROOT/'FILE_LIST.txt')
    (ROOT/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n'
                                        for p in sorted(all_files)))
    print(len(names),'files;',len(all_files),'checksums')
if __name__=='__main__':main()
