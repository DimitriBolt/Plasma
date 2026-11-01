"""Saved-data report assertions and PDF metadata/text/font checks."""
import hashlib,json,re,subprocess
import pandas as pd
from .upstream import ROOT

def main():
    R=ROOT/'results';pdf=ROOT/'report/07_Integrated_Model.pdf'
    run=lambda *x:subprocess.check_output(x,text=True)
    info=run('pdfinfo',str(pdf));text=run('pdftotext','-layout',str(pdf),'-')
    fonts=run('pdffonts',str(pdf));lp=ROOT/'report/07_Integrated_Model.log'
    if lp.exists():log=lp.read_text()
    else:
        original=json.loads((R/'full_archive_final_verification.json').read_text())
        assert hashlib.sha256(pdf.read_bytes()).hexdigest()==original['pdf_sha256']
        assert not original['final_latex_undefined_references'] and not original['overfull_boxes']
        log='Verified original PDF and archived final compilation diagnostics.'
    pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M)[1])
    assert pages==13 and 'Author:          Dimitri Bolt' in info
    assert 'Overfull' not in log and 'undefined' not in log.lower()
    assert not re.search(r'Codex|SHA256|review_handoff|not_delivered|delivered_and_verified|git push',text)
    assert 'Mini-Research 7' in text and '6000' in text and 'unresolved' in text
    assert all('yes' in line for line in fonts.splitlines()[2:])
    s=pd.read_csv(R/'stable_pairs.csv',float_precision='round_trip')
    maintex=(R/'table_main.tex').read_text()
    for r in s[(s.level==3)&(s.order==3)&(s.fluence_J_m2==10)].itertuples():
        if r.completed_pair:assert f'{r.global_delta_Tmax_K*1e6:.4g}' in maintex
    assert len(list((ROOT/'figures').glob('*.jpg')))==10
    assert not list((ROOT/'figures').glob('*.png')) and not list((ROOT/'figures').glob('*.svg'))
    notebook=json.loads((ROOT/'07_Integrated_Model.ipynb').read_text())
    assert len({c['id'] for c in notebook['cells']})==len(notebook['cells'])
    testlog=(R/'tests.txt').read_text();assert testlog.rstrip().endswith('OK')
    count=int(re.search(r'Ran (\d+) tests',testlog)[1])
    comp=json.loads((R/'storage_compaction.json').read_text());assert len(comp)==115 and all(x['scalar_arrays_unchanged'] for x in comp)
    preservation=json.loads((R/'input_preservation.json').read_text());assert preservation['verified_files']==1209 and preservation['all_unchanged']
    result=dict(pdf_pages=pages,main_pages=11,appendix_pages=2,pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
        tests=count,historical_full_archive_constitutive_replays=115,unchanged_upstream_files=1209,
        figures_JPG=10,report_value_assertions_passed=True,fonts_embedded=True,
        final_latex_undefined_references=False,overfull_boxes=False,
        visual_review='all 13 pages reviewed; changed pages rechecked after final presentation edits',
        physical_validation='not established; reduced-model and formal results only')
    if (R/'clean_archive_check.json').exists():
        clean=json.loads((R/'clean_archive_check.json').read_text());assert clean['all_passed'] and clean['tests_passed']==31
        result['clean_archive_snapshot']=clean['scientific_snapshot_sha'];result['clean_archive_all_checks_passed']=True
    if (R/'publication_scope.json').exists():
        result['profile']='compact-publication-v1';result['current_shipped_baseline_replays']=4
        result['historical_clean_archive_record_scope']='full local snapshot; not the current compact export'
    (R/'final_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
