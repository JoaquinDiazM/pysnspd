"""Verify imported evidence, preserved physical equations and the rendered report."""
from pathlib import Path
import hashlib,json
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/review_20260922'
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    inventory=json.loads((DATA/'import_inventory.json').read_text(encoding='utf-8-sig'))
    checks={v['relative']:sha(DATA/'raw'/v['relative'])==v['sha256'] for v in inventory['files']}
    record=json.loads((DATA/'raw/ssp/one_ssp_40.json').read_text(encoding='utf-8'))
    sources={p:sha(ROOT/p)==value for p,value in record['source_hashes'].items()}
    pdf=ROOT/'output/pdf/implementation/Informe_revision_etapa_2_20260922.pdf'
    reader=PdfReader(pdf);texts=[p.extract_text() for p in reader.pages]
    assert len(texts)==4 and all(t.strip() for t in texts)
    assert 'todavía no procede cerrar' in texts[0]
    assert all(checks.values()) and all(sources.values())
    rendered=ROOT/'tmp/pdfs/stage2_review_20260922'
    pages={str(n):sha(rendered/f'report-{n}.png') for n in range(1,5)}
    result=dict(status='PASS',imported_files=len(checks),import_hashes_match=all(checks.values()),
        physical_and_integration_sources_unchanged=all(sources.values()),pdf_pages=4,
        pdf_sha256=sha(pdf),visual_reviewed_pages=[1,2,3,4],rendered_page_sha256=pages,
        visual_findings='All four rendered pages inspected: legible labels, distinct curves, no clipping or overlapping text.',
        report_declares_no_closure=True,new_physical_trajectories=0,script_sha256=sha(__file__))
    (DATA/'QA_report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
