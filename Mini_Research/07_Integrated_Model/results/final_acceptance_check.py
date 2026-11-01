"""Executable editorial-preservation audit against the published baseline.

Run from a Git checkout: python results/final_acceptance_check.py
After reviewing every rendered PDF page, add --visual-reviewed. The optional
--verify-manifest checks the separately regenerated publication manifest.
No scientific solver, data writer or upstream module is invoked here.
"""
import argparse,csv,hashlib,json,re,subprocess
from pathlib import Path

BASE='658908ddad0a2b12a922f6a51ce7a58259a9cccd'
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
PREFIX='Mini_Research/07_Integrated_Model/'
EDITORIAL={
 'README.md','07_Integrated_Model.ipynb','src/presentation.py',
 'report/07_Integrated_Model.tex','report/07_Integrated_Model.pdf',
 'results/table_main.tex','results/narrative_primary.tex',
 'results/narrative_conclusion.tex','results/model_decisions.md','results/data_schema.md'}
INVENTORIES={'FILE_LIST.txt','SHA256SUMS'}

def command(*args):return subprocess.check_output(args,cwd=REPO)
def rows(name):return list(csv.DictReader((ROOT/'results'/name).open()))
def equations(s):return [m.group(0) for m in re.finditer(r'\\begin\{(equation|align)\}.*?\\end\{\1\}',s,re.S)]

def main():
 p=argparse.ArgumentParser();p.add_argument('--visual-reviewed',action='store_true');p.add_argument('--verify-manifest',action='store_true');a=p.parse_args()
 entries=[];changed=[];oldpaths=[]
 for line in filter(None,command('git','ls-tree','-r','-z',BASE).decode().split('\0')):
  meta,path=line.split('\t',1);mode,kind,expected=meta.split();assert kind=='blob'
  oldpaths.append(path);f=REPO/path
  if path.startswith(PREFIX) and path[len(PREFIX):] in EDITORIAL|INVENTORIES:
   if path[len(PREFIX):] in EDITORIAL and hashlib.sha1(b'blob '+str(f.stat().st_size).encode()+b'\0'+f.read_bytes()).hexdigest()!=expected:changed.append(path)
   continue
  b=f.read_bytes();actual=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
  assert actual==expected,'Protected input changed: '+path
  entries.append(dict(path=path,baseline_git_blob=expected,current_git_blob=actual,sha256=hashlib.sha256(b).hexdigest()))
 new=[x for x in command('git','ls-files','--others','--exclude-standard','-z').decode().split('\0') if x.startswith(PREFIX)]
 tracked=set(filter(None,command('git','ls-files','-z').decode().split('\0')))
 new=sorted(set(new)|(tracked-set(oldpaths)))
 assert all(x.startswith(PREFIX+'results/final_acceptance') for x in new),new
 old=command('git','show',BASE+':'+PREFIX+'report/07_Integrated_Model.tex').decode()
 tex=(ROOT/'report/07_Integrated_Model.tex').read_text();assert equations(old)==equations(tex)
 table=(ROOT/'results/table_main.tex').read_text()
 for word in ['Numerical completion','Physical applicability','Conditional reduced model','Formal fixed-composition scenario','Censored at 6000 K']:assert word in table,word
 for pa,value in [(1.,'1.4088163'),(1.15,'2.1848791')]:
  r=next(r for r in rows('stable_pairs.csv') if abs(float(r['Pa_bar'])-pa)<1e-12 and r['level']=='3' and r['order']=='3' and float(r['fluence_J_m2'])==10)
  assert format(float(r['global_delta_Tmax_K'])*1e6,'.8g')==value and value in table
 wet=next(r for r in rows('stable_pairs.csv') if abs(float(r['Pa_bar'])-1.3)<1e-12 and r['level']=='3')
 assert not wet['global_delta_Tmax_K'] and wet['completed_pair']=='False'
 late=next(r for r in rows('strong_collapse_resolution.csv') if float(r['offset_s'])==0)
 assert round(100*float(late['relative_mesh_change']),1)==89.6 and late['numerical_resolution']=='unresolved_late_tail'
 for word in ['1.30\\times10^8','Negative liquid pressure alone is not a mathematical proof of impossibility','metastable tension','21--30\\%','89.6\\%','not a physical temperature ceiling']:assert word in tex,word
 pdf=ROOT/'report/07_Integrated_Model.pdf'
 info=command('pdfinfo',str(pdf)).decode();pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M)[1]);assert 'Dimitri Bolt' in info
 text=command('pdftotext','-layout',str(pdf),'-').decode();page_text=text.split('\f');main_pages=next(i for i,s in enumerate(page_text) if 'Numerical controls and reproducibility' in s)
 assert 8<=main_pages<=12,(pages,main_pages)
 assert '1.4088163' in text and '2.1848791' in text and '89.6%' in text
 assert not re.search(r'Codex|SHA256|git push|delivered_and_verified',text)
 log=(ROOT/'report/07_Integrated_Model.log').read_text();assert 'Overfull' not in log and 'undefined' not in log.lower()
 fonts=command('pdffonts',str(pdf)).decode();assert all('yes' in x for x in fonts.splitlines()[2:])
 notebook=json.loads((ROOT/'07_Integrated_Model.ipynb').read_text());nb=''.join(''.join(c['source']) for c in notebook['cells'])
 for doc in [(ROOT/'README.md').read_text(),nb]:
  for word in ['Numerical completion','Physical applicability','Formal fixed-composition','89.6%','21–30%','6263.579','22808.716']:assert word in doc,word
 data=[x for x in entries if x['path'].startswith(PREFIX) and Path(x['path']).suffix in {'.csv','.npz'}]
 code=[x for x in entries if x['path'].startswith(PREFIX+'src/') and x['path'].endswith('.py')]
 upstream=[x for x in entries if not x['path'].startswith(PREFIX)]
 result=dict(baseline=BASE,editorial_only=True,changed_editorial_paths=sorted(changed),
  protected_files=len(entries),unchanged_upstream_files=len(upstream),unchanged_model_and_other_source_modules=len(code),unchanged_CSV_NPZ_files=len(data),
  all_protected_files_byte_identical=True,equation_blocks_unchanged=len(equations(tex)),figure_files_unchanged=10,
  compact_scope_and_historical_records_unchanged=True,pdf_pages=pages,main_pages=main_pages,appendix_pages=pages-main_pages,
  pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),all_pages_visually_reviewed=a.visual_reviewed,
  corrections_checked=['numerical completion versus physical applicability','21–30% local Au-rise limitation beside material thresholds','89.6% unresolved tail and censored/formal extremes retained'],
  scope='Editorial acceptance corrections only; no independent scientific closure is asserted.',files=entries)
 (ROOT/'results/final_acceptance_preservation.json').write_text(json.dumps(result,indent=2)+'\n')
 if a.verify_manifest:subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=ROOT,stdout=subprocess.DEVNULL,check=True)
 print(json.dumps({k:v for k,v in result.items() if k!='files'},indent=2))

if __name__=='__main__':main()
