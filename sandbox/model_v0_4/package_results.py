"""Package model edition 0.4 and the current independently versioned notebook."""
from pathlib import Path
import hashlib,json,zipfile,subprocess,datetime,re,shutil

ROOT=Path(__file__).resolve().parents[2]
DOCS=ROOT/'docs/modelo_v0_4'
DATA=DOCS/'verificaciones'
OUTPUT=ROOT/'output/pdf/modelo_v0_4'
MANIFEST=DATA/'manifiesto_v0_4.json'
ZIP=OUTPUT/'pysnspd_documentos_v0_4.zip'
LEARNING=ROOT/'sandbox/model_v0_4/learning_revision'
REGISTRY=DOCS/'E_registro_de_revisiones.json'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def learning_metadata():
    registry=json.loads(REGISTRY.read_text(encoding='utf-8'))
    source=(DOCS/'E_cuaderno_de_aprendizaje_v0_4.md').read_text(encoding='utf-8')
    headings=re.findall(r'^# E\.\d+\..*?— (E\d+)\s*$',source,re.M)
    assert headings==list(registry['classes']), 'La ruta de clases no coincide con el registro'
    answers=re.findall(r'^\*\*Respuesta (E\d+\.(?:A)?\d+):\*\*[ \t]*([^\n]*)',source,re.M)
    assert len(answers)==len({name for name,_ in answers}), 'Identificadores de respuesta repetidos'
    assert {name.split('.')[0] for name,_ in answers}==set(registry['classes'])
    pending=sum(text.strip()=='_Escribir aquí._' for _,text in answers)
    assert pending==registry['activities_pending'], 'Respuestas pendientes y registro difieren'
    assert not set(registry['classes']).intersection(registry['retired'])
    return registry, len(answers), pending

def preserve_previous_delivery(registry):
    """Archive existing outputs once by content hash before replacing either file.

    Historical ZIPs live outside DOCS and are not nested into the new bundle.
    The initial 0.4 manifest predates notebook_revision; its supplied registry
    identifies that delivery as the previous notebook snapshot.
    """
    if not MANIFEST.is_file():
        assert not ZIP.exists(), 'Hay un ZIP previo sin manifiesto que identifique su entrega'
        return None
    old=json.loads(MANIFEST.read_text(encoding='utf-8'))
    old_revision=old.get('notebook_revision',Path(registry['previous_snapshot']).name)
    assert re.fullmatch(r'E-r\d+',old_revision), old_revision
    snapshot=Path(old_revision)/('entrega-'+digest(MANIFEST)[:12])
    manifest_copy=DOCS/'historial'/snapshot/MANIFEST.name
    zip_copy=OUTPUT/'historial'/snapshot/ZIP.name
    result={'notebook_revision':old_revision}
    for original,target,key in [(MANIFEST,manifest_copy,'manifest'),(ZIP,zip_copy,'zip')]:
        if not original.is_file():
            continue
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            assert digest(target)==digest(original), f'El historial no se sobrescribe: {target}'
        else:
            shutil.copy2(original,target)
        result[key]={'path':target.relative_to(ROOT).as_posix(),'sha256':digest(target)}
    return result

def main():
    registry,activity_count,pending=learning_metadata()
    qa=json.loads((DATA/'QA_estructura.json').read_text(encoding='utf-8'))
    visual={'revision':registry['model_edition'],'notebook_revision':registry['notebook_revision'],
            'status':'approved','documents':[]}
    for letter in 'ABCDE':
        review=json.loads((DATA/f'QA_visual_{letter}.json').read_text(encoding='utf-8'))
        if letter=='A':
            assert review['resultado']=='aprobado'
            row={'pdf':Path(review['pdf']).name,'pdf_sha256':review['sha256_pdf'],
                 'pages_reviewed':[p['pagina'] for p in review['paginas'] if p['visto'] and p['resultado']=='aprobado'],
                 'method':review['metodo']}
        elif letter=='C':
            assert review['status']=='aprobado' and review['all_pages_viewed']
            row={'pdf':review['document'],'pdf_sha256':review['sha256'],
                 'pages_reviewed':review['pages_inspected'],'method':review['method']}
        else:
            assert review['status']=='approved'
            row={'pdf':Path(review['pdf']).name,'pdf_sha256':review['pdf_sha256'],
                 'pages_reviewed':review['pages_reviewed'],'method':review['method']}
        row['record']=f'QA_visual_{letter}.json'
        visual['documents'].append(row)
    (DATA/'QA_visual.json').write_text(json.dumps(visual,indent=2,ensure_ascii=False),encoding='utf-8')
    assert len(qa)==5
    for entry in qa:
        pdf=OUTPUT/entry['pdf']
        assert digest(pdf)==entry['pdf_sha256']
        v=next(x for x in visual['documents'] if x['pdf']==pdf.name)
        assert v['pdf_sha256']==digest(pdf) and v['pages_reviewed']==list(range(1,entry['pages']+1))
    previous_delivery=preserve_previous_delivery(registry)
    files=sorted(set([p for p in DOCS.rglob('*') if p.is_file() and p!=MANIFEST and '__pycache__' not in p.parts]
        +list((ROOT/'sandbox/model_v0_4').glob('*.py'))
        +[p for p in LEARNING.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix in {'.py','.md','.json'}]
        +list(OUTPUT.glob('[ABCDE]*.pdf'))))
    source_specs=[
        ('observaciones_v0_3_adjuntas',Path('C:/Users/joaqu/.codex/attachments/e06cc86c-982e-480e-a9de-315cf912cb67/pasted-text.txt')),
        ('memoria_02.pdf',ROOT/'tmp/pdfs/modelo_v0_2/sources/memoria_02.pdf'),
        ('Allmaras_Thesis_Final.pdf',ROOT/'tmp/pdfs/modelo_v0_3/sources/Allmaras_Thesis_Final.pdf')]
    provenance=[{'name':name,'sha256':digest(path),'redistributed':False} for name,path in source_specs if path.is_file()]
    figs=set()
    for entry in qa:
        md=DOCS/(Path(entry['pdf']).stem+'.md')
        figs.update(re.findall(r'!\[[^\]]*\]\(([^)]+)\)',md.read_text(encoding='utf-8')))
    assert figs and all((DOCS/name).is_file() for name in figs), 'Falta una figura citada'
    learning_state='; '.join(f"{code}: {item['learning_state']} (r{item['revision']})"
                            for code,item in registry['classes'].items())
    learning_state+=f"; retiradas: {', '.join(registry['retired'])}; {pending} respuestas pendientes"
    manifest={
        'revision':registry['model_edition'],'notebook_revision':registry['notebook_revision'],
        'date':registry['date'],'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'repository_base':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'scope':'Documentos, figuras y pruebas ligeras; sin modificación del solver ni ejecución PRE/SS/photon',
        'pdf_pages':sum(item['pages'] for item in qa),'equation_labels':sum(item['equation_labels'] for item in qa),
        'distinct_figures':len(figs),'pedagogical_topics':len(registry['classes']),
        'activities_total':activity_count,'activities_pending':pending,
        'class_revisions':registry['classes'],'retired_classes':registry['retired'],
        'learning_state':learning_state,'previous_delivery_snapshot':previous_delivery,
        'model_admission':'Candidato experimental restringido; no se promueve a producción. Fuerza de núcleo sensible, símbolo principal negativo en estados probados, entrada fonónica inadmisible y tiempos cinéticos sin identificar.',
        'independent_diagnostics':['checks_a_v04.py','checks_b_v04.py','checks_c_v04.py','checks_d_v04.py',
          'learning_revision/figures_foundations.py','learning_revision/figures_e04.py',
          'learning_revision/checks_modes.py','learning_revision/figures_dft_dfpt.py'],
        'learning_build':['learning_revision/build_learning.py'],
        'historical_learning_assembly':{'path':'sandbox/model_v0_4/learning_revision/assemble_learning.py',
          'purpose':'Montaje editorial de E-r02 desde fragmentos; reconstruir sólo sobre una copia para preservar la edición personal del Markdown.'},
        'historical_learning_script':{'path':'sandbox/model_v0_4/learning_v04.py',
          'purpose':'Reproduce la entrega E-r01; no ejecutar para regenerar las figuras actuales.'},
        'material_source':{'path_on_geminga':'/home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat',
            'sha256':'e94f11273e6c42b7c8fad978b78d119166e241bec1614ee77e4176906fa8f5df',
            'redistributed':False,'admitted_for_absolute_rates':False},
        'remote_root':'geminga:/home/jdiaz/pysnspd',
        'delivery_status':{'notebook':'E-r02; comprobar verificacion_entrega.json externo para la integridad local/remota vigente.',
          'documents_A_D':'Edición 0.4; corrección editorial en D para referir el cuaderno E-r02. Ecuaciones y resultados conservados.',
          'software_release':'v1.0.0 conserva la implementación de la memoria; el modelo actualizado es documentación.'},
        'source_documents':provenance,
        'primary_reference_urls':['https://arxiv.org/pdf/1611.06060','https://doi.org/10.1103/PhysRevB.101.094507',
          'https://arxiv.org/pdf/1909.00992v1','https://arxiv.org/pdf/2501.13791v3',
          'https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf',
          'https://ocw.mit.edu/courses/18-325-topics-in-applied-mathematics-waves-and-imaging-fall-2015/834e98a77f64597a422c205b123cd6d7_MIT18_325F15_Appendix_A.pdf',
          'https://topocondmat.org/test/w1_topointro/0d.html',
          'https://www.quantum-espresso.org/Doc/ph_user_guide/node19.html',
          'https://arxiv.org/abs/cond-mat/0012092'],
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
