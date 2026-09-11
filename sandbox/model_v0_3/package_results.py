"""Package revision 0.3 with hashes and the completed visual/structural QA."""
from pathlib import Path
import hashlib,json,zipfile,subprocess,datetime,re

ROOT=Path(__file__).resolve().parents[2]
DOCS=ROOT/'docs/modelo_v0_3'
DATA=DOCS/'verificaciones'
OUTPUT=ROOT/'output/pdf/modelo_v0_3'
MANIFEST=DATA/'manifiesto_v0_3.json'
ZIP=OUTPUT/'pysnspd_documentos_v0_3.zip'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    qa=json.loads((DATA/'QA_estructura.json').read_text(encoding='utf-8'))
    visual=json.loads((DATA/'QA_visual.json').read_text(encoding='utf-8'))
    assert visual['status']=='approved'
    assert len(qa)==5
    for entry in qa:
        pdf=OUTPUT/entry['pdf']
        assert digest(pdf)==entry['pdf_sha256']
        v=next(x for x in visual['documents'] if x['pdf']==pdf.name)
        assert v['pdf_sha256']==digest(pdf) and v['pages_reviewed']==list(range(1,entry['pages']+1))
    files=sorted([p for p in DOCS.rglob('*') if p.is_file() and p!=MANIFEST and '__pycache__' not in p.parts]
        +list((ROOT/'sandbox/model_v0_3').glob('*.py'))+list(OUTPUT.glob('[ABCDE]*.pdf')))
    source_specs=[
        ('observaciones_v0_3',Path('C:/Users/joaqu/.codex/attachments/e06cc86c-982e-480e-a9de-315cf912cb67/pasted-text.txt')),
        ('memoria_02.pdf',ROOT/'tmp/pdfs/modelo_v0_2/sources/memoria_02.pdf'),
        ('Allmaras_Thesis_Final.pdf',ROOT/'tmp/pdfs/modelo_v0_3/sources/Allmaras_Thesis_Final.pdf')]
    provenance=[{'name':name,'sha256':digest(path),'redistributed':False} for name,path in source_specs if path.is_file()]
    figs=set()
    for md in DOCS.glob('[ABCDE]*.md'):
        figs.update(re.findall(r'!\[[^\]]*\]\(([^)]+)\)',md.read_text(encoding='utf-8')))
    assert len(figs)==17
    manifest={
        'revision':'0.3','date':'2026-09-08','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'repository_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'scope':'Documentos, figuras y pruebas ligeras; sin modificación del solver ni ejecución PRE/SS/photon',
        'pdf_pages':sum(item['pages'] for item in qa),'equation_labels':sum(item['equation_labels'] for item in qa),
        'distinct_figures':len(figs),'pedagogical_topics':5,'learning_state':'E01–E05 abiertos; respuestas y evaluación pendientes',
        'remote_root':'geminga:/home/jdiaz/pysnspd','source_documents':provenance,
        'primary_reference_urls':['https://arxiv.org/pdf/1611.06060','https://doi.org/10.1103/PhysRevB.101.094507',
          'https://arxiv.org/pdf/1909.00992v1','https://arxiv.org/pdf/2501.13791v3',
          'https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf','https://arxiv.org/pdf/1805.00130v2'],
        'files':[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':digest(p)} for p in files]}
    MANIFEST.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    with zipfile.ZipFile(ZIP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[MANIFEST]:z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(ZIP) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(files)+1
        for entry in manifest['files']:
            assert hashlib.sha256(z.read(entry['path'])).hexdigest()==entry['sha256']
    print(json.dumps({'zip':str(ZIP),'bytes':ZIP.stat().st_size,'files':len(files)+1,
      'pdf_pages':manifest['pdf_pages'],'figures':len(figs),'sha256':digest(ZIP)},indent=2))

if __name__=='__main__':main()
