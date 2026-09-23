"""Render the results-only stage-2 results report from its reviewed content record.

This assembles existing results; it does not calculate a physical trajectory.
Run from any directory with reportlab installed. Fonts are system Arial on
Windows or DejaVu Sans on Linux. Every page is explicitly partitioned for QA.
"""
from __future__ import annotations

import json
import re
from functools import partial
from pathlib import Path
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak,
)

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'docs/implementation/stage2/closure_prep_20260922'
OUTPUT = ROOT / 'output/pdf/implementation/Informe_validacion_dos_celdas_etapa_2_20260922.pdf'
NAVY = colors.HexColor('#16334b')
TEAL = colors.HexColor('#007c83')
INK = colors.HexColor('#202d39')
GRAY = colors.HexColor('#596775')
WIDTH = A4[0] - 34*mm


def fonts():
    candidates = [
        (Path('C:/Windows/Fonts/arial.ttf'), Path('C:/Windows/Fonts/arialbd.ttf')),
        (Path('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'),
         Path('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf')),
    ]
    try:
        import matplotlib
        base = Path(matplotlib.get_data_path()) / 'fonts/ttf'
        candidates.append((base/'DejaVuSans.ttf', base/'DejaVuSans-Bold.ttf'))
    except ImportError:
        pass
    for normal, bold in candidates:
        if normal.exists() and bold.exists():
            pdfmetrics.registerFont(TTFont('Body', normal))
            pdfmetrics.registerFont(TTFont('BodyBold', bold))
            pdfmetrics.registerFontFamily('Body', normal='Body', bold='BodyBold')
            return
    raise RuntimeError('Install Arial or DejaVu Sans to render the report.')


def main():
    content = json.loads((DATA/'report_content.json').read_text(encoding='utf-8'))
    fonts()
    styles = {
        'title': ParagraphStyle('title', fontName='BodyBold', fontSize=23,
                                leading=27, textColor=NAVY, spaceAfter=10),
        'heading': ParagraphStyle('heading', fontName='BodyBold', fontSize=15,
                                  leading=19, textColor=NAVY, spaceAfter=8),
        'subheading': ParagraphStyle('subheading', fontName='BodyBold', fontSize=10.5,
                                     leading=14, textColor=TEAL, spaceBefore=7, spaceAfter=5),
        'body': ParagraphStyle('body', fontName='Body', fontSize=9.5,
                               leading=13.5, textColor=INK, spaceAfter=7),
        'small': ParagraphStyle('small', fontName='Body', fontSize=8,
                                leading=10.7, textColor=GRAY, spaceAfter=5),
        'cell': ParagraphStyle('cell', fontName='Body', fontSize=8.3,
                               leading=11, textColor=INK),
        'tablehead': ParagraphStyle('tablehead', fontName='BodyBold', fontSize=8.3,
                                    leading=11, textColor=colors.white),
    }
    story, markdown = [], ['# '+content['title'], '', content['subtitle'], '']
    def p(text, style='body'):
        # Use font-supported baseline glyphs for mathematical scripts.
        superscripts = str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹', '-0123456789')
        text = re.sub(r'[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+',
                      lambda m: '<super>'+m.group().translate(superscripts)+'</super>', text)
        subscripts = str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789')
        text = re.sub(r'[₀₁₂₃₄₅₆₇₈₉]+',
                      lambda m: '<sub>'+m.group().translate(subscripts)+'</sub>', text)
        return Paragraph(text, styles[style])
    for number, page in enumerate(content['pages']):
        if number:
            story.append(PageBreak())
        story.append(p(page['title'], 'title' if number == 0 else 'heading'))
        markdown += ['## '+page['title'], '']
        if number == 0:
            story.append(p(content['subtitle'], 'small'))
        for block in page['blocks']:
            kind = block['type']
            if kind in ('body', 'small', 'subheading'):
                story.append(p(block['text'], kind))
                markdown += [block['text'], '']
            elif kind == 'table':
                rows = block['rows']
                rendered = [[p(escape(str(cell)), 'tablehead' if i == 0 else 'cell')
                             for cell in row] for i, row in enumerate(rows)]
                widths = block.get('widths', [1/len(rows[0])]*len(rows[0]))
                table = Table(rendered, colWidths=[WIDTH*w for w in widths], repeatRows=1)
                table.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), NAVY),
                    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor('#edf4f7'), colors.white]),
                    ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                    ('LEFTPADDING', (0, 0), (-1, -1), 7),
                    ('RIGHTPADDING', (0, 0), (-1, -1), 7),
                    ('TOPPADDING', (0, 0), (-1, -1), 7),
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 7),
                ]))
                story += [table, Spacer(1, 8)]
                markdown += ['| '+' | '.join(map(str, rows[0]))+' |',
                             '| '+' | '.join(['---']*len(rows[0]))+' |']
                markdown += ['| '+' | '.join(map(str, row))+' |' for row in rows[1:]]
                markdown += ['']
            elif kind == 'image':
                from PIL import Image as PILImage
                path = ROOT/block['path']
                with PILImage.open(path) as im:
                    w, h = im.size
                display_width = WIDTH*block.get('width_fraction', 1)
                story.append(Image(str(path), width=display_width, height=display_width*h/w))
                story.append(p(block['caption'], 'small'))
                markdown += [f"![{block['caption']}]({path.relative_to(DATA).as_posix()})", '']
            else:
                raise ValueError('Unknown content block: '+kind)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    def footer(canvas, doc):
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(.7)
        canvas.line(17*mm, 15*mm, A4[0]-17*mm, 15*mm)
        canvas.setFont('Body', 8)
        canvas.setFillColor(GRAY)
        canvas.drawString(17*mm, 10.5*mm, 'pySNSPD · Etapa 2 - Celdas acopladas · 22 septiembre 2026')
        canvas.drawRightString(A4[0]-17*mm, 10.5*mm, str(doc.page))
    doc = SimpleDocTemplate(str(OUTPUT), pagesize=A4, rightMargin=17*mm,
                            leftMargin=17*mm, topMargin=16*mm, bottomMargin=22*mm,
                            title=content['title'], author='pySNSPD',
                            pageCompression=1)
    doc.build(story, onFirstPage=footer, onLaterPages=footer,
              canvasmaker=partial(Canvas, invariant=1))
    (DATA/'Informe_validacion_dos_celdas_etapa_2_20260922.md').write_text('\n'.join(markdown).rstrip()+'\n', encoding='utf-8')
    print(OUTPUT)


if __name__ == '__main__':
    main()
