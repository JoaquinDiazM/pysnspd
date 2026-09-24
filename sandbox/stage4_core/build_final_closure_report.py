"""Build the final physical stage-4 report from archived results, without solves."""
from pathlib import Path
import hashlib
import json
from html import escape
from matplotlib.font_manager import findfont
from PIL import Image
from pypdf import PdfReader
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'docs/implementation/stage4/closure_20260924'
PDF = ROOT / 'output/pdf/implementation/Informe_cierre_etapa_4.pdf'


def main():
    content = json.loads((DATA/'report_content.json').read_text(encoding='utf8'))
    for name, query in [('DV', 'DejaVu Sans'), ('DVB', 'DejaVu Sans:weight=bold')]:
        pdfmetrics.registerFont(TTFont(name, findfont(query)))
    pdfmetrics.registerFontFamily('DV', normal='DV', bold='DVB', italic='DV', boldItalic='DVB')
    W, H = 595.276, 841.89
    L, WIDTH = 38., 519.276
    body = ParagraphStyle('body', fontName='DV', fontSize=9.2, leading=12.4,
                          textColor=HexColor('#183247'))
    small = ParagraphStyle('small', parent=body, fontSize=8.1, leading=10.8)
    title = ParagraphStyle('title', parent=body, fontName='DVB', fontSize=18, leading=22)
    PDF.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(PDF), pagesize=(W, H), pageCompression=1)
    c.setTitle('Etapa 4 - Cierre de desarrollo y resultados físicos')
    c.setAuthor('pySNSPD')
    bottoms, image_hashes = [], {}
    markdown = ['# Etapa 4: cierre de desarrollo y resultados físicos', '',
                '24 de septiembre de 2026. Sin fotón. Producción y v1.0.0 sin cambios.', '']

    def paragraph(text, y, style=body):
        q = Paragraph(text, style)
        _, h = q.wrap(WIDTH, 1000)
        q.drawOn(c, L, y-h)
        return y-h

    for number, page in enumerate(content['pages'], 1):
        c.setFillColor(HexColor('#135b70')); c.rect(0, H-28, W, 28, fill=1, stroke=0)
        c.setFillColor(HexColor('#ffffff')); c.setFont('DVB', 9)
        c.drawString(L, H-18, 'pySNSPD | ETAPA 4 | INFORME FINAL DE DESARROLLO')
        y = paragraph(page['title'], H-48, title)-8
        y = paragraph(page['lead'], y)-10
        markdown += ['## '+page['title'], '', page['lead'], '']
        if page.get('figure'):
            path = DATA/'figures'/page['figure']
            image_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
            with Image.open(path) as im:
                height = min(page.get('figure_height', 367), WIDTH*im.height/im.width)
                width = height*im.width/im.height
            c.drawImage(str(path), L+(WIDTH-width)/2, y-height, width, height, mask='auto')
            y -= height+7
            y = paragraph(page['caption'], y, small)-12
            markdown += [f"![{page['title']}](figures/{page['figure']})", '', page['caption'], '']
        for text in page.get('paragraphs', []):
            y = paragraph(text, y)-9
            markdown += [text, '']
        if page.get('table'):
            rows = [[Paragraph(escape(str(cell)), small) for cell in row] for row in page['table']]
            widths = page.get('column_widths', [WIDTH/len(rows[0])]*len(rows[0]))
            table = Table(rows, colWidths=widths)
            table.setStyle(TableStyle([
                ('BACKGROUND', (0,0), (-1,0), HexColor('#e3eef1')),
                ('VALIGN', (0,0), (-1,-1), 'TOP'),
                ('BOTTOMPADDING', (0,0), (-1,-1), 6),
                ('TOPPADDING', (0,0), (-1,-1), 6),
                ('LINEBELOW', (0,0), (-1,0), .7, HexColor('#7a9ba7')),
                ('LINEBELOW', (0,1), (-1,-1), .3, HexColor('#d8e2e6'))]))
            _, height = table.wrap(WIDTH, H)
            table.drawOn(c, L, y-height)
            y -= height+10
            markdown += ['| '+' | '.join(map(str,page['table'][0]))+' |',
                         '| '+' | '.join(['---']*len(rows[0]))+' |']
            markdown += ['| '+' | '.join(map(str,row))+' |' for row in page['table'][1:]]
            markdown += ['']
        for text in page.get('after_table', []):
            y = paragraph(text, y, small)-7
            markdown += [text, '']
        if y < 42:
            raise ValueError(f'Page {number} overflows: bottom={y:.1f}')
        bottoms.append(y)
        c.setFillColor(HexColor('#627681')); c.setFont('DV', 7.5)
        c.drawString(L, 23, '24 septiembre 2026 | Controles sin fotón | No es una validación de latencia')
        c.drawRightString(W-L, 23, f'{number} / {len(content["pages"])}')
        c.showPage()
    c.save()
    (DATA/'Informe_cierre_etapa_4.md').write_text('\n'.join(markdown), encoding='utf8')
    reader = PdfReader(PDF)
    assert len(reader.pages) == len(content['pages'])
    context = dict(pdf_path=str(PDF.relative_to(ROOT)),
        pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest(), pages=len(reader.pages),
        page_content_bottoms_pt=bottoms, physics_solves_performed=0,
        stage4_development_closed=True, full_contract_completed=False,
        photon_dynamics_admitted=False, production_changed=False,
        figure_sha256=image_hashes, visual_inspection='PENDING',
        builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (DATA/'report_qa.json').write_text(json.dumps(context, indent=2)+'\n', encoding='utf8')
    print(json.dumps(dict(pdf=str(PDF), pages=len(reader.pages), bottoms=bottoms)))


if __name__ == '__main__':
    main()
