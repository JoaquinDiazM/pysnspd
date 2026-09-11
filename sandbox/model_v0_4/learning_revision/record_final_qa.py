"""Record completed visual inspections and bind the E delivery to its inputs."""
from pathlib import Path
import datetime
import hashlib
import json
import re
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[3]
DOCS=ROOT/'docs/modelo_v0_4'
WORK=Path(__file__).resolve().parent
PDF=ROOT/'output/pdf/modelo_v0_4/E_cuaderno_de_aprendizaje_v0_4.pdf'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
digest=sha(PDF)
assert len(PdfReader(PDF).pages)==26
root_review={'pdf_sha256':digest,'reviewer':'root','pages_actually_viewed':list(range(19,27)),
 'method':'Inspección visual directa de los ocho PNG finales de Poppler a 1500 px, en contactos de dos páginas.',
 'defects':[], 'notes':['Ecuaciones, DOS, comparación de alpha2F y lambda legibles y sin solapamientos.',
 'La última página conserva espacio para observaciones de aprendizaje.']}
(WORK/'qa_final_19_26.json').write_text(json.dumps(root_review,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
reviews=[json.loads((WORK/f'qa_final_{name}.json').read_text(encoding='utf-8')) for name in ['01_09','10_18','19_26']]
seen=[]
for review in reviews:
    assert review.get('sha256',review.get('pdf_sha256')).lower()==digest
    assert not review['defects']
    seen.extend(review['pages_actually_viewed'])
assert sorted(seen)==list(range(1,27)) and len(seen)==len(set(seen))
qa={'revision':'0.4','notebook_revision':'E-r02','pdf':PDF.relative_to(ROOT).as_posix(),
 'pdf_sha256':digest,'pages_reviewed':sorted(seen),'status':'approved',
 'method':'Todas las páginas finales renderizadas con Poppler e inspeccionadas directamente a 1500 px; revisión repartida entre tres revisores. Se corrigieron el índice residual, la tabla partida y saltos que producían espacios excesivos.',
 'findings':[], 'scope':'Calidad documental y verificación de ejemplos pedagógicos; no certifica aprendizaje ni valida un material o el detector.',
 'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(DOCS/'verificaciones/QA_visual_E.json').write_text(json.dumps(qa,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
md=DOCS/'E_cuaderno_de_aprendizaje_v0_4.md'
source=md.read_text(encoding='utf-8')
figs=[DOCS/p for p in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',source)]
files=[md,PDF,DOCS/'E_VERSIONADO.md',DOCS/'E_registro_de_revisiones.json',
       DOCS/'latex'/f'{md.stem}.tex',DOCS/'latex/header_learning.tex',
       DOCS/'verificaciones/QA_visual_E.json',*figs,*sorted(WORK.glob('*.py'))]
manifest={'model_edition':'0.4','notebook_revision':'E-r02','pages':26,'lessons_active':8,
 'activities_pending':24,'equations':len(re.findall(r'\\tag\{',source)),
 'figures_in_document':len(figs),'canonical_source':md.relative_to(ROOT).as_posix(),
 'files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in files]}
(DOCS/'verificaciones/manifiesto_E_r02.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps({'status':'approved','pdf_sha256':digest,'pages_reviewed':len(seen),'figures':len(figs)}))
