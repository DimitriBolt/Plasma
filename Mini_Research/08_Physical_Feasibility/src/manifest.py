"""Hash only MR8 and Final_Report scientific/supporting files; exclude auxiliaries."""
from pathlib import Path
import hashlib
from .upstream import ROOT,REPO

EXCLUDE={'.aux','.log','.out','.toc','.fls','.fdb_latexmk','.gz','.pyc'}

def main():
    for directory in [ROOT,REPO/'Final_Report']:
        files=sorted(p for p in directory.rglob('*') if p.is_file() and p.suffix not in EXCLUDE
           and not any(x in p.parts for x in ['.venv','__pycache__','.pytest_cache'])
           and p.name not in ['SHA256SUMS','FILE_LIST.txt'])
        # The bounded computation log is evidence, unlike compiler logs.
        log=directory/'results/computation.log'
        if log.exists():files=sorted(set(files+[log]))
        (directory/'FILE_LIST.txt').write_text('\n'.join(str(p.relative_to(directory)) for p in files)+
            '\nFILE_LIST.txt\nSHA256SUMS\n')
        files=sorted(files+[directory/'FILE_LIST.txt'])
        (directory/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(directory))+'\n' for p in files))
        print(directory.name,len(files),'manifest entries')

if __name__=='__main__':main()
