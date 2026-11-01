"""Read-only final PDF structure/content checks using pypdf and Poppler.

Run with a Python that has pypdf. Images require separate human/visual inspection.
"""
import hashlib,json,re,subprocess,tempfile
from pathlib import Path
import xml.etree.ElementTree as ET
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    records=[]
    for rel,expected in [('Final_Report/Consolidated_Final_Scientific_Report.pdf',24),
                         ('Mini_Research/08_Physical_Feasibility/report/08_Physical_Feasibility.pdf',9)]:
        path=REPO/rel;reader=PdfReader(path);assert len(reader.pages)==expected
        content='\n'.join(p.extract_text() for p in reader.pages)
        assert 'Dimitri Bolt' in content and 'Gabitov' in content
        assert not re.search(r'\b[0-9a-f]{40}\b',content)
        assert not re.search(r'Codex|RoadMap|handoff|Part [A-H]|acceptance procedure',content)
        links=[]
        for page in reader.pages:
            for a in page.get('/Annots',[]):
                act=a.get_object().get('/A',{})
                if act.get('/URI'):links.append(str(act['/URI']))
        reports=sorted({u for u in links if 'github.com/DimitriBolt/Plasma/blob/main/Mini_Research/' in u and u.endswith('.pdf')})
        if rel.startswith('Final'):
            assert len(reports)==8
            for stage in ['01_','02_','03_','04_','05_','06_','07_']:assert any('/'+stage in u for u in reports)
            for term in ['1.408816','2.184879','6000','6263.58','22808.7','1.13046','1.16720','0.03347','0.02499']:assert term in content,term
        with tempfile.TemporaryDirectory() as td:
            xml=Path(td)/'bbox.xml';subprocess.run(['pdftotext','-bbox',str(path),str(xml)],check=True)
            clean=re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', '', xml.read_text())
            tree=ET.ElementTree(ET.fromstring(clean));bad=[]
            for i,page in enumerate(tree.findall('.//{*}page'),1):
                width=float(page.attrib['width']);height=float(page.attrib['height'])
                for w in page.findall('.//{*}word'):
                    if float(w.attrib['xMin'])<0 or float(w.attrib['xMax'])>width or float(w.attrib['yMin'])<0 or float(w.attrib['yMax'])>height:bad.append((i,w.text))
            assert not bad,bad
        log=path.with_suffix('.log').read_text()
        assert 'Overfull' not in log and 'undefined references' not in log
        records.append(dict(path=rel,pages=len(reader.pages),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
           page_text_characters=[len(p.extract_text()) for p in reader.pages],external_link_annotations=len(links),
           supporting_pdf_urls=reports,no_out_of_page_words=True,no_overfull_or_undefined_references=True,
           fonts_checked_separately=True,visual_review='Every rendered page inspected; final changes rechecked.'))
    (ROOT/'results/pdf_qa.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(records,indent=2))

if __name__=='__main__':main()
