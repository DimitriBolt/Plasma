"""Manifest only deliverable files; preserve unrelated user files."""
import hashlib
from .model import ROOT


def main():
    excluded={'.aux','.log','.out','.fdb_latexmk','.fls','.gz.tmp','.pyc'}
    files=sorted(p for p in ROOT.rglob('*') if p.is_file()
        and not any(x.startswith('.') for x in p.relative_to(ROOT).parts if x!='.gitignore')
        and not any(x=='__pycache__' for x in p.relative_to(ROOT).parts)
        and p.suffix not in excluded and not p.name.endswith('.synctex.gz')
        and p.name not in ('SHA256SUMS','FILE_LIST.txt'))
    listing=sorted([str(p.relative_to(ROOT)) for p in files]+['FILE_LIST.txt','SHA256SUMS','.gitignore'])
    # Include the local ignore policy in the integrity set, once only.
    files=sorted(set(files)|{ROOT/'.gitignore'})
    (ROOT/'FILE_LIST.txt').write_text('\n'.join(sorted(set(listing)))+'\n')
    files.append(ROOT/'FILE_LIST.txt')
    (ROOT/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n' for p in sorted(files)))


if __name__=='__main__': main()
