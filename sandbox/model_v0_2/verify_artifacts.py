"""Check final PDF structure, equation labels, figures and compilation logs."""
from pathlib import Path
import re,json,hashlib
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[2]
DOCS=ROOT/'docs/modelo_v0_2'
PDFS=ROOT/'output/pdf/modelo_v0_2'

def main():
    summary=[]
    for md in sorted(DOCS.glob('[ABCD]_*_v0_2.md')):
        text=md.read_text(encoding='utf-8')
        pdf=PDFS/(md.stem+'.pdf')
        reader=PdfReader(pdf)
        content='\n'.join(page.extract_text() or '' for page in reader.pages)
        compact=re.sub(r'\s+','',content)
        tags=re.findall(r'\\tag\{([^}]+)\}',text)
        missing=[tag for tag in tags if '('+tag+')' not in compact]
        assert not missing,(md.name,missing)
        assert len(tags)==len(set(tags)),md.name
        source_figures=re.findall(r'!\[[^\]]*\]\(([^)]+)\)',text)
        embedded=sum(len(page.images) for page in reader.pages)
        assert embedded==len(source_figures),(md.name,embedded,source_figures)
        for fig in source_figures:assert (DOCS/fig).is_file()
        assert '\ufffd' not in content,pdf.name
        log=ROOT/'tmp/pdfs/modelo_v0_2'/(md.stem+'.log')
        problems=[]
        if log.is_file():
            problems=re.findall(r'^.*(?:Overfull|Underfull|Missing character|Undefined control|undefined references).*$',log.read_text(encoding='utf-8',errors='replace'),re.M)
            assert not problems,(md.name,problems)
        summary.append({'document':md.name,'pdf':pdf.name,'pages':len(reader.pages),'equation_labels':len(tags),
                        'missing_equation_labels':missing,'embedded_figures':embedded,'compile_layout_warnings':problems,
                        'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()})
    (DOCS/'verificaciones/QA_estructura.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
