"""Concise read-only status for the active Geminga command notebook."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/review_20260922'
if __name__=='__main__':
    record=json.loads((DATA/'independent_audit.json').read_text(encoding='utf-8'))
    print('Etapa 2: no cerrada. No hay admisión espacial ni de producción.')
    for statement in record['conclusions'][:3]:print('- '+statement)
    print('Informe: output/pdf/implementation/Informe_revision_etapa_2_20260922.pdf')
    print('Único plan pendiente: docs/implementation/stage2/review_20260922/manual_time_plan.json')
