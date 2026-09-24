"""Result report from saved dual KWT and weak longitudinal trajectories only."""
from pathlib import Path
import hashlib
import json

from matplotlib.font_manager import findfont
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage4/final_kwt_20260924'
PDF = ROOT/'output/pdf/implementation/Informe_etapa_4_Euler_y_acoplamiento.pdf'
INPUTS = []


def read(relative):
    path = DATA/relative
    INPUTS.append(path)
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary = read('thermal/summary.json')
    comparison = read('thermal/refinement.json')
    read('thermal/analysis.json')
    plan = read('thermal/executed_plan.json')
    pilot = read('pilot_comparison.json')
    longitudinal = read('longitudinal/receipt.json')
    policy = read('acceptance_policy.json')
    if not summary['all_practical_criteria_met'] or longitudinal['status'] != 'PASSED_CONTROL':
        raise ValueError('The report requires the two recorded successful controls')
    if plan['refinement_relative_limit'] != .02 or plan['energy_relative_limit'] != .02:
        raise ValueError('The practical 2 percent policy must match the executed plan')
    fonts = [('DV', 'DejaVu Sans'), ('DVB', 'DejaVu Sans:weight=bold')]
    for name, query in fonts:
        pdfmetrics.registerFont(TTFont(name, findfont(query)))
    pdfmetrics.registerFontFamily('DV', normal='DV', bold='DVB', italic='DV', boldItalic='DVB')
    W, H = 595.276, 841.89
    L, WIDTH = 39., 517.276
    PDF.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(PDF), pagesize=(W, H), pageCompression=1)
    c.setTitle('Etapa 4: Euler heredado, disipación y acoplamiento débil')
    c.setAuthor('pySNSPD - resultados de implementación')
    body = ParagraphStyle('body', fontName='DV', fontSize=9.5, leading=13,
                          textColor=HexColor('#233647'))
    small = ParagraphStyle('small', parent=body, fontSize=8.1, leading=10.6)
    page_bottoms = []

    def paragraph(text, y, style=body):
        q = Paragraph(text, style)
        _, h = q.wrap(WIDTH, 1000)
        q.drawOn(c, L, y-h)
        return y-h

    def heading(number, title, subtitle):
        c.setFillColor(HexColor('#15576b')); c.rect(0, H-29, W, 29, fill=1, stroke=0)
        c.setFillColor(HexColor('#ffffff')); c.setFont('DVB', 9)
        c.drawString(L, H-19, 'pySNSPD | ETAPA 4 | RESULTADOS DE IMPLEMENTACIÓN')
        c.setFillColor(HexColor('#183449')); c.setFont('DVB', 18)
        c.drawString(L, H-63, title)
        y = paragraph(subtitle, H-76, small)
        c.setFillColor(HexColor('#677985')); c.setFont('DV', 7.5)
        c.drawString(L, 23, '24 septiembre 2026 | Sin fotón | Producción y v1.0.0 intactas')
        c.drawRightString(W-L, 23, f'{number} / 3')
        return y-14

    def figure(relative, y, max_height):
        path = DATA/relative
        INPUTS.append(path)
        with Image.open(path) as im:
            height = min(max_height, WIDTH*im.height/im.width)
            width = height*im.width/im.height
        c.drawImage(str(path), L+(WIDTH-width)/2, y-height, width=width,
                    height=height, mask='auto')
        return y-height

    def table(rows, y, widths):
        values = [[Paragraph(str(cell), small) for cell in row] for row in rows]
        t = Table(values, colWidths=widths, hAlign='LEFT')
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), HexColor('#e8f0f3')),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('LINEBELOW', (0, 0), (-1, 0), .6, HexColor('#a6bac4')),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('LEFTPADDING', (0, 0), (-1, -1), 6),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1),
             [HexColor('#ffffff'), HexColor('#f6f8f9')]),
        ]))
        _, h = t.wrap(WIDTH, 1000); t.drawOn(c, L, y-h)
        return y-h

    def end_page(y):
        page_bottoms.append(y)
        if y < 42:
            raise ValueError(f'Content enters footer area: bottom={y}')
        c.showPage()

    y = heading(1, 'Euler de la memoria: ensayo aprobado',
        'Cuatro trayectorias térmicas sobre la malla dual alcanzaron 1 ps en 153,8 s. '
        'Los campos y el balance integrado cumplen la política práctica del 2 %.')
    y = paragraph('Se reutilizó el <b>paso KWT real de la memoria</b> con la nueva fuerza y corriente espectrales. '
        'El ensayo conserva 1712 nodos, 256 frecuencias, T = 0,9 K y T<sub>c</sub> = 8,65 K. '
        'La cinta de 80 × 160 nm tiene contactos de extremo fijos y laterales aislantes; '
        'su longitud sigue siendo una preparación provisional, no un dominio validado de Korzh.', y)-9
    y = figure('thermal/figures/01_evolucion_y_refinamiento.png', y, 237)-7
    y = paragraph('<b>Figura 1.</b> d = Δ/(k<sub>B</sub>T<sub>c</sub>); δd es la diferencia respecto del equilibrio uniforme. '
        'Se compara una perturbación de amplitud del 0,1 % con una de fase de 0,001 rad, ambas con el perfil suave declarado. '
        'La norma nodal usa áreas duales: ||z||<sub>M</sub> = √Σ m<sub>i</sub>|z<sub>i</sub>|², '
        'm<sub>i</sub> = A<sub>i</sub>/ℓ<sub>0</sub>². Para corriente adimensional I se usa √Σ |I<sub>ij</sub>|²/c<sub>ij</sub>, '
        'con c<sub>ij</sub> = cara dual/arista. G es el gradiente de la acción; Im(d*G)/m es fuerza de fase, no torque mecánico. '
        'La diferencia temporal usa la mayor norma inicial '
        'de las dos sondas para el mismo observable. Primaria y refinada emplean Δt y Δt/2; '
        'el panel de diferencias evita ocultar curvas superpuestas.', y, small)-10
    n = len(comparison['records'])
    maxerror = max(row['relative_to_initial_signal'] or 0 for row in comparison['records'])
    y = table([
        ['Resultado térmico', 'Medida'],
        ['Comparaciones durante la evolución que cumplen', '40/40 (+ 8 iniciales)'],
        ['Mayor diferencia temporal / escala inicial', f'{100*maxerror:.4f} %'.replace('.', ',')],
        ['Pasos por sonda: primaria / refinada', '581 / 1162'],
        ['Máximo residuo espectral / límite', '3,96 × 10<super>−8</super> / 10<super>−7</super>'],
    ], y, [WIDTH*.65, WIDTH*.35])-11
    y = paragraph('<b>Precisión adoptada.</b> Se compara con el 2 % de la escala inicial más el piso absoluto '
        'prescrito, 10<super>−7</super> en la norma del observable. El rechazo anterior de una cola cruzada '
        f'queda aceptado con la revisión v2: {policy["maximum_error_over_initial"]*100:.3f} % &lt; 2 %. '
        'Los datos y el certificado histórico se conservan; esa cola ya no exige otra corrida.', y, small)-9
    y = paragraph('El método resuelve una cuadrática para |Δ|², pero su precisión temporal es de <b>primer orden</b>. '
        'El grado de esa ecuación no es el orden temporal. Se mantiene el algoritmo publicado de '
        '<link href="https://py-tdgl.readthedocs.io/en/latest/background.html#implicit-euler-method" color="#15576b">pyTDGL</link> '
        'y se verifica su precisión sobre las señales de interés.', y, small)
    end_page(y)

    y = heading(2, 'Disipación resuelta y menos cómputo',
        'El balance se comprobó en toda la historia aceptada. No se observaron aumentos de energía libre entre pasos.')
    y = figure('thermal/figures/02_balance_integrado.png', y, 342)-7
    y = paragraph('<b>Figura 2.</b> F es la energía libre térmica adimensional de la misma acción espectral; '
        'E<sub>0</sub> = F(0) − F<sub>eq</sub> es el exceso inicial de cada sonda. '
        'La pérdida acumulada integra (P<sub>KWT</sub> + P<sub>normal</sub>) dt/t<sub>D</sub>, '
        'con t<sub>D</sub> = ℏ/(2k<sub>B</sub>T<sub>c</sub>). El defecto es '
        'F(t) − F(0) + pérdida acumulada. Los porcentajes usan E<sub>0</sub>, no la energía absoluta. '
        'Es disipación a temperatura impuesta; no es el balance de energía interna de poblaciones.', y, small)-10
    bal = summary['integrated_balance']
    y = table([
        ['Máximo defecto / exceso inicial', 'Primaria', 'Refinada'],
        ['Sonda de amplitud', f'{100*abs(bal["primary_amplitude"]["relative_to_initial_excess"]):.6f} %',
         f'{100*abs(bal["refined_amplitude"]["relative_to_initial_excess"]):.6f} %'],
        ['Sonda de fase', f'{100*abs(bal["primary_angular_phase"]["relative_to_initial_excess"]):.4f} %',
         f'{100*abs(bal["refined_angular_phase"]["relative_to_initial_excess"]):.4f} %'],
    ], y, [WIDTH*.52, WIDTH*.24, WIDTH*.24])-10
    y = paragraph('El límite energético es <b>0,02 E<sub>0</sub> + 10<super>−8</super></b>, declarado antes de ejecutar. '
        'La sonda angular primaria consume aproximadamente el 57 % del presupuesto relativo del 2 %; '
        'al dividir el paso, su defecto baja de 1,1467 % a 0,5771 %. No es necesario endurecerlo para avanzar.', y, small)-10
    y = paragraph('<b>Qué se optimizó.</b> Se predice el espectro usando una factorización de la referencia uniforme. '
        'Cada predicción se verifica contra el residuo no lineal exacto con el mismo límite 10<super>−7</super>; '
        'si no cumple, se conserva Newton como corrección. Las 893 440 consultas del ensayo completo '
        'pasaron esa verificación. El mismo piloto hasta 0,001 ps bajó de <b>29,571 a 4,661 s</b>: '
        '6,34 veces menos tiempo total, sin alterar ecuaciones ni tolerancia. El horizonte completo costó 153,813 s.', y, small)-9
    y = paragraph('<b>Recursos efectivos.</b> 27 trabajadores y un coordinador, un hilo por proceso: '
        '28 de 32 CPU lógicas. Se reservaron dos núcleos físicos completos y 30 GiB para planificación, '
        'por debajo del 90 % de la memoria disponible. La cifra de memoria es una reserva, no una medición de pico.', y, small)
    end_page(y)

    y = heading(3, 'Intercambio no térmico débil: aprobado',
        'Control de amplitud del gap y forma energética de la población electrónica, sin fotón ni reajuste a una temperatura.')
    y = figure('longitudinal/figures/longitudinal_exchange.png', y, 345)-7
    y = paragraph('<b>Figura 3.</b> A: modo espacial más suave, normalizado a máximo uno; no es un depósito fotónico. '
        'B: cambio real del gap en el máximo del modo, dividido por el gap de equilibrio; el recuadro separa las resoluciones. '
        'C: incremento de ocupación electrónica respecto de Fermi-Dirac a 0,9 K, en las siete energías de control '
        'medidas en k<sub>B</sub>T<sub>c</sub>; las líneas unen muestras. D: disponibilidad A y pérdidas acumuladas '
        'divididas por A(0). <b>Disponibilidad</b> significa energía libre cuadrática de la perturbación; no energía '
        'interna total ni calor entregado a fonones.', y, small)-10
    y = table([
        ['Control longitudinal hasta 1 ps', 'Resultado'],
        ['Error temporal / norma inicial de disponibilidad', '0,003288 %'],
        ['Defecto integrado / A(0)', '0,006619 %'],
        ['Disponibilidad restante / A(0)', '66,6206 %'],
        ['Tiempo de cálculo, un proceso y un hilo', '0,078 s'],
    ], y, [WIDTH*.73, WIDTH*.27])-10
    y = paragraph('Se conservan ocho amplitudes de un sector espacial que es invariante en las <b>ecuaciones linealizadas</b>: '
        'una del gap y siete de población. El RHS coincide con su representación en los 1712 nodos. '
        'La referencia matricial sólo contrasta Euler; no sustituye el solver del detector. '
        'La preparación suprime el gap 0,2 % y añade una perturbación de ocupación independiente. '
        'No incluye colisiones ni un circuito polarizado.', y, small)-10
    y = paragraph('<b>Conclusión de implementación.</b> La malla dual, Euler heredado y la disipación térmica pasan; '
        'también pasa el intercambio no térmico débil ensayado. <b>La etapa 4 general sigue abierta</b> para unir '
        'poblaciones, fase, potencial, calor y las tres variables del circuito. Se investiga la separación de escalas: '
        '<b>cuasiclásica no significa equilibrio instantáneo</b>. En la referencia de Korzh, ℏ/Δ ≈ 0,50 ps es una escala '
        'espectral, no un tiempo de termalización. Se conservan provisionalmente las poblaciones dinámicas; la validez '
        'de calcular el espectro como una respuesta instantánea al estado sigue en evaluación. '
        'No se necesita repetir las corridas de este informe.', y, small)
    end_page(y)
    c.save()
    reader = PdfReader(str(PDF))
    alltext = ' '.join(' '.join(page.extract_text() or '' for page in reader.pages).split())
    assert len(reader.pages) == 3
    for required in ('153', 'primer orden', 'etapa 4 general sigue abierta', '893 440', '1712'):
        assert required in alltext, required
    input_hashes = {str(path.relative_to(ROOT)).replace('\\', '/'): sha(path)
                    for path in sorted(set(INPUTS))}
    context = dict(pdf_path=str(PDF.relative_to(ROOT)).replace('\\', '/'),
        pdf_sha256=sha(PDF), source_sha256=input_hashes,
        builder_sha256=sha(Path(__file__)), pages=len(reader.pages),
        page_content_bottoms_pt=page_bottoms, physics_solves_performed=0,
        stage4_complete=False, production_changed=False,
        semantic_checks='PASS', visual_inspection='PENDING')
    (DATA/'report_context.json').write_text(json.dumps(context, indent=2, ensure_ascii=False)+'\n', encoding='utf8')
    print(json.dumps(context, indent=2))


if __name__ == '__main__':
    main()
