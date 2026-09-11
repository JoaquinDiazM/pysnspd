"""Package the editable sources, PDFs, figures and small checks for revision 0.2."""
from pathlib import Path
import hashlib,json,zipfile,subprocess,datetime

ROOT=Path(__file__).resolve().parents[2]
DOCS=ROOT/'docs/modelo_v0_2'
DATA=DOCS/'verificaciones'
OUTPUT=ROOT/'output/pdf/modelo_v0_2'
MANIFEST=DATA/'manifiesto_v0_2.json'
ZIP=OUTPUT/'pysnspd_documentos_v0_2.zip'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def payload():
    return sorted([p for p in DOCS.rglob('*') if p.is_file() and p!=MANIFEST and '__pycache__' not in p.parts]
        +list((ROOT/'sandbox/model_v0_2').glob('*.py'))+list(OUTPUT.glob('[ABCD]*.pdf')))

def main():
    files=payload()
    qa=json.loads((DATA/'QA_estructura.json').read_text(encoding='utf-8'))
    originals=[]
    downloads=Path('C:/Users/joaqu/Downloads')
    for stem in ('A_microscopia_y_funcional','B_cinetica_energia_y_temperaturas','C_condensado_corriente_y_senal','D_sintesis_plan_y_verificaciones'):
        for suffix in ('.md','.pdf'):
            path=downloads/(stem+suffix)
            if path.is_file():originals.append({'name':path.name,'sha256':digest(path),'redistributed':False})
    memory=ROOT/'tmp/pdfs/modelo_v0_2/sources/memoria_02.pdf'
    if memory.is_file():originals.append({'name':'memoria_02.pdf','source':'geminga:/home/jdiaz/memoria/main/memoria_02.pdf','sha256':digest(memory),'redistributed':False})
    manifest={
        'revision':'0.2','date':'2026-09-08','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'repository_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'scope':'Documentos, pruebas ligeras y figuras; sin cambios de solver ni ejecución PRE/SS/photon',
        'pdf_pages':sum(item['pages'] for item in qa),'distinct_figures':len(list((DOCS/'figuras').glob('*.png'))),
        'pdf_visual_review':'50 páginas finales inspeccionadas; ecuaciones y figuras legibles, sin recortes ni superposiciones',
        'source_documents':originals,
        'primary_reference_urls':['https://arxiv.org/pdf/1611.06060','https://doi.org/10.1103/PhysRevB.101.094507','https://arxiv.org/pdf/1909.00992v1'],
        'files':[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':digest(p)} for p in files]
    }
    MANIFEST.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    with zipfile.ZipFile(ZIP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files+[MANIFEST]:z.write(p,p.relative_to(ROOT).as_posix())
    with zipfile.ZipFile(ZIP) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(files)+1
        for entry in manifest['files']:
            assert hashlib.sha256(z.read(entry['path'])).hexdigest()==entry['sha256']
    print(json.dumps({'zip':str(ZIP),'bytes':ZIP.stat().st_size,'files':len(files)+1,'pdf_pages':manifest['pdf_pages'],'distinct_figures':manifest['distinct_figures'],'sha256':digest(ZIP)},indent=2))

if __name__=='__main__':main()
