"""Compile maintained mathematical TeX sources; never regenerate their content.

No simulation is executed. Builds, renders, logs and replaced PDF history stay
in scratch. Publication is explicit and only allowed after the stage-5 closure.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess

REPO=Path('/home/jdiaz/pysnspd')
SOURCE=REPO/'sandbox/stage5_math_documents'
DEFAULT_BUILD=Path('/home/jdiaz/scratch/stage5_math_documents_20261011')
DOCUMENTS=(
 ('Matematica_etapa5_principales_2D','Principal.pdf'),
 ('Matematica_etapa5_auxiliares','Auxiliar.pdf'),
 ('Matematica_etapa5_precomputados','Precomputo.pdf'),
 ('Matematica_etapa5_parametros','Parametros.pdf'),
)
PHYSICAL_SOURCES=(
 'pysnspd/experimental/effective_bcs_functional.py',
 'pysnspd/experimental/effective_bcs_kinetics.py',
 'pysnspd/experimental/effective_bcs_kinetics_numba.py',
 'pysnspd/experimental/effective_bcs_ports.py',
 'pysnspd/experimental/effective_bcs_preparation.py',
 'pysnspd/experimental/effective_coupled.py',
 'pysnspd/experimental/effective_material.py',
 'pysnspd/experimental/electrical_ports.py',
 'pysnspd/experimental/trigger_window.py',
 'docs/implementation/MODELO_VIGENTE.md',
 'docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md',
)
def digest(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def preserve_pdf(path,history):
 if not path.exists():
  return
 sha=digest(path)
 target=history/(path.stem+'_'+sha[:16]+'.pdf')
 if not target.exists():
  shutil.copy2(path,target)
 if digest(target)!=sha:
  raise RuntimeError('Historical PDF verification failed: '+str(path))

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--build-root',type=Path,default=DEFAULT_BUILD)
 parser.add_argument('--documents',nargs='+',choices=tuple(short for _,short in DOCUMENTS),
  help='Compile only these exact maintained PDF filenames; default: all documents')
 parser.add_argument('--publish',action='store_true',
  help='After visual QA, copy compiled PDFs to the exact maintained Matematica filenames')
 parser.add_argument('--stage5-status',choices=('in_progress','complete'),default='in_progress')
 args=parser.parse_args()
 build=args.build_root.resolve()
 if not build.is_relative_to(Path('/home/jdiaz/scratch')):
  raise ValueError('Builds and renders must remain under /home/jdiaz/scratch')
 if args.publish and args.stage5_status!='complete':
  raise ValueError('Do not publish mathematical PDFs before complete stage-5 closure')
 build.mkdir(parents=True,exist_ok=True)
 compiler='/usr/bin/pdflatex'
 environment=os.environ.copy()
 for variable,leaf in (('TEXMFVAR','texlive-var'),('TEXMFCONFIG','texlive-config'),('TEXMFHOME','texmf')):
  environment.setdefault(variable,str(Path('/home/jdiaz/scratch/tmp/pysnspd')/leaf))
 receipt=dict(schema='pysnspd.mathematical_documents.v2',date='2026-10-12',
  model='Effective neutral native BCS/shared radial gradient/KWT/kinetics/CM',
  stage5_status_declared_by_caller=args.stage5_status,compiled_in_scratch=str(build),
  physical_simulations_launched=0,tex_sources_rewritten=False,published=False,
  pdf_operation_marker='Caller follows its PDF artifact workflow; not asserted by this script',
  physical_source_sha256={p:digest(REPO/p)for p in PHYSICAL_SOURCES if(REPO/p).exists()},
  outputs=[])
 for name,short_name in DOCUMENTS:
  if args.documents is not None and short_name not in args.documents:
   continue
  source=SOURCE/(name+'.tex')
  before=digest(source)
  for attempt in (1,2):
   command=[compiler,'-interaction=nonstopmode','-halt-on-error',
     '-output-directory='+str(build),str(source)]
   proc=subprocess.run(command,cwd=build,env=environment,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=120)
   (build/(name+f'.compile_{attempt}.log')).write_bytes(proc.stdout)
   if proc.returncode:
    print(proc.stdout.decode('utf8',errors='replace')[-6000:])
    raise RuntimeError('Compilation failed: '+name)
  if digest(source)!=before:
   raise RuntimeError('Compilation unexpectedly changed maintained TeX source')
  pdf=build/(name+'.pdf')
  log=(build/(name+'.log')).read_text(errors='replace')
  warnings=[line for line in log.splitlines()if'Overfull' in line or'Underfull' in line]
  receipt['outputs'].append(dict(name=name,short_name=short_name,scratch_pdf=str(pdf),
   bytes=pdf.stat().st_size,sha256=digest(pdf),tex_sha256=before,layout_warnings=warnings))
 if args.publish:
  history=build/'previous_final_pdfs'
  history.mkdir(exist_ok=True)
  mathdir=REPO/'Matematica'
  mathdir.mkdir(exist_ok=True)
  for record in receipt['outputs']:
   targets=(mathdir/record['short_name'],)
   for target in targets:
    preserve_pdf(target,history)
   for target in targets:
    shutil.copy2(record['scratch_pdf'],target)
    if digest(target)!=record['sha256']:
     raise RuntimeError('Published PDF verification failed: '+str(target))
   record['published_paths']=[str(target)for target in targets]
  receipt['published']=True
  (SOURCE/'mathematics_document_receipt.json').write_text(
   json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
 (build/'mathematics_document_receipt.json').write_text(
  json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
 print(json.dumps(receipt,indent=2,ensure_ascii=False))

if __name__=='__main__':
 main()
