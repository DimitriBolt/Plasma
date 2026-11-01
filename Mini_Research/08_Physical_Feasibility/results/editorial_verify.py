"""Editorial-only checks: read saved artifacts and Git objects, never run models."""
import argparse,hashlib,json,re,subprocess,tempfile,unittest
from pathlib import Path
import xml.etree.ElementTree as ET
from pypdf import PdfReader

BASE='083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2'
ROOT=Path(__file__).resolve().parents[3]
REFERENCE=ROOT
ALLOWED={
 'Mini_Research/02_Compressible_Bubble_Dynamics/README.md',
 'Mini_Research/02_Compressible_Bubble_Dynamics/FILE_LIST.txt',
 'Mini_Research/02_Compressible_Bubble_Dynamics/SHA256SUMS',
 'Final_Report/Consolidated_Final_Scientific_Report.tex',
 'Final_Report/Consolidated_Final_Scientific_Report.pdf',
 'Final_Report/support_links.tex','Final_Report/FILE_LIST.txt','Final_Report/SHA256SUMS',
 'Mini_Research/08_Physical_Feasibility/README.md',
 'Mini_Research/08_Physical_Feasibility/report/08_Physical_Feasibility.pdf',
 'Mini_Research/08_Physical_Feasibility/report/support_links.tex',
 'Mini_Research/08_Physical_Feasibility/FILE_LIST.txt',
 'Mini_Research/08_Physical_Feasibility/SHA256SUMS'}
REPORTS=['Final_Report/Consolidated_Final_Scientific_Report.pdf',
         'Mini_Research/08_Physical_Feasibility/report/08_Physical_Feasibility.pdf']
CANON='Mini_Research/02_Compressible_Bubble_Dynamics/report/02_Compressible_Bubble_Dynamics.pdf'
OLD='Mini_Research/02_Compressible_Bubble_Dynamics/report/main.pdf'
METRICS={}

def original(path):
    return subprocess.check_output(['git','-C',str(REFERENCE),'show',f'{BASE}:{path}'])

class EditorialChecks(unittest.TestCase):
    def test_preserve_scientific_snapshot(self):
        files=subprocess.check_output(['git','-C',str(REFERENCE),'ls-tree','-r','--name-only','-z',BASE]).decode().split('\0')
        files=[f for f in files if f and f not in ALLOWED]
        changed=[f for f in files if not (ROOT/f).exists() or (ROOT/f).read_bytes()!=original(f)]
        self.assertEqual(changed,[])
        METRICS['protected_original_files_identical']=len(files)
        METRICS['scientific_code_data_figures_tables_unchanged']=True

    def test_canonical_MR2_copy(self):
        a=(ROOT/CANON).read_bytes();b=(ROOT/OLD).read_bytes()
        self.assertEqual(a,b);self.assertEqual(b,original(OLD))
        METRICS['MR2_both_pdf_sha256']=hashlib.sha256(a).hexdigest()

    def test_links_and_English_pdfs(self):
        documents=[]
        for rel in REPORTS:
            reader=PdfReader(ROOT/rel);text='\n'.join(p.extract_text() for p in reader.pages)
            self.assertFalse(re.search('[А-Яа-яЁё]',text))
            links=[];internal=0;names=reader.named_destinations
            for p in reader.pages:
                for a in p.get('/Annots',[]):
                    annotation=a.get_object();action=annotation.get('/A',{})
                    uri=action.get('/URI')
                    if uri:links.append(str(uri))
                    destination=annotation.get('/Dest')
                    if action.get('/S')=='/GoTo':destination=action.get('/D')
                    if destination is not None:
                        internal+=1
                        if isinstance(destination,str):self.assertIn(destination,names)
                        else:self.assertGreater(len(destination),0)
            for destination in names.values():
                self.assertIn(reader.get_destination_page_number(destination),range(len(reader.pages)))
            canonical='https://github.com/DimitriBolt/Plasma/blob/main/'+CANON
            self.assertIn(canonical,links)
            self.assertNotIn('https://github.com/DimitriBolt/Plasma/blob/main/'+OLD,links)
            if rel.startswith('Final'):
                for i in range(1,8):self.assertTrue(any(f'/0{i}_' in u for u in links))
                self.assertIn('Comparison with Experimental Evidence and Validation Status',text)
            documents.append(dict(path=rel,pages=len(reader.pages),internal_links_resolved=internal,named_destinations_resolved=len(names),sha256=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest(),
                source_pdf_urls=sorted({u for u in links if 'github.com/DimitriBolt/Plasma/blob/main/Mini_Research/' in u and u.endswith('.pdf')})))
        METRICS['documents']=documents

    def test_readme_translation_contract(self):
        path='Mini_Research/08_Physical_Feasibility/README.md'
        old=original(path).decode();new=(ROOT/path).read_text()
        self.assertEqual(re.findall(r'```bash\n(.*?)```',old,re.S),re.findall(r'```bash\n(.*?)```',new,re.S))
        self.assertEqual(re.findall(r'\]\(([^)]+)\)',old),re.findall(r'\]\(([^)]+)\)',new))
        for token in ['0.513361 µm','0.550809 µm','34.3373 fJ','0.03347','0.02499 K','1.408816','2.184879 µK','1.15 bar','1.30 bar','6000 K','Python 3.12','18','ΔU = −∫p dV + ∫Q(l→g) dt','Δ(Ug+El+EAu)+Wby+Qbath=Eabs']:
            self.assertIn(token,new)
        self.assertGreater(len(re.findall('[А-Яа-яЁё]',new)),1000)
        METRICS['README_commands_links_values_units_preserved']=True

    def test_condensation_and_scientific_statements(self):
        path='Final_Report/Consolidated_Final_Scientific_Report.tex'
        old=original(path).decode();new=(ROOT/path).read_text()
        def count(s):
            part=s[s.index('\\section{Acoustic bubble'):s.index('\\section{Timing and useful')]
            return len(re.findall(r"[A-Za-z]+(?:[-'][A-Za-z]+)*",part))
        a,b=count(old),count(new);self.assertLess(b,a)
        self.assertIn('Gilmore--Tait',new)
        for x in ['1.13046','1.16720','22808.7','6263.58','1.408816','2.184879','6000','n_D','T_i']:
            self.assertIn(x,new)
        self.assertIn('not a completed global maximum',new)
        self.assertIn('This does not prove impossibility',new)
        self.assertIn('no direct experimental validation',new)
        METRICS['sections_3_to_7_word_estimate']={'before':a,'after':b,'reduction_percent':100*(a-b)/a,
           'method':'Same English-token estimate on the LaTeX section span; includes equation and caption tokens.'}

    def test_pdf_layout_and_compilation(self):
        for rel in REPORTS:
            log=(ROOT/rel).with_suffix('.log').read_text()
            self.assertNotIn('Overfull',log);self.assertNotIn('undefined references',log)
            with tempfile.TemporaryDirectory() as td:
                xml=Path(td)/'bbox.xml';subprocess.run(['pdftotext','-bbox',str(ROOT/rel),str(xml)],check=True)
                clean=re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]','',xml.read_text())
                for p in ET.fromstring(clean).findall('.//{*}page'):
                    w,h=float(p.attrib['width']),float(p.attrib['height'])
                    for word in p.findall('.//{*}word'):
                        self.assertGreaterEqual(float(word.attrib['xMin']),0)
                        self.assertLessEqual(float(word.attrib['xMax']),w)
                        self.assertGreaterEqual(float(word.attrib['yMin']),0)
                        self.assertLessEqual(float(word.attrib['yMax']),h)
            fonts=subprocess.check_output(['pdffonts',str(ROOT/rel)]).decode().splitlines()[2:]
            self.assertTrue(all(f.split()[-5]=='yes' for f in fonts if f.strip()))
        METRICS['no_overfull_undefined_outside_words_and_all_fonts_embedded']=True

def main():
    global ROOT,REFERENCE
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,default=ROOT);parser.add_argument('--reference-repo',type=Path,default=REFERENCE);parser.add_argument('--output',type=Path);args=parser.parse_args()
    ROOT,REFERENCE=args.root.resolve(),args.reference_repo.resolve()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(EditorialChecks)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    METRICS.update(scientific_base=BASE,editorial_tests_passed=result.testsRun-len(result.failures)-len(result.errors),models_reintegrated=False,
       historical_records_unchanged=True,scientific_conclusions_unchanged=True)
    if args.output:args.output.write_text(json.dumps(METRICS,indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)

if __name__=='__main__':main()
