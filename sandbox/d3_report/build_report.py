"""Rebuild the historical D3 report from external cached results and renders.

Inputs remain under ``tmp/pdfs/d3_report``; see the adjacent README for the
required files and the original Windows font dependency. This script only
assembles the report and does not run a detector simulation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tmp" / "pdfs" / "d3_report"
ASSETS = WORK / "assets"
OUTPUT = ROOT / "output" / "pdf" / "Informe_D3_proyeccion_energetica.pdf"
CACHE = WORK / "D3_energy_projection_cache.npz"

PAGE_W, PAGE_H = landscape(A4)
MARGIN = 12 * mm
NAVY = colors.HexColor("#17324D")
BLUE = colors.HexColor("#1769AA")
TEAL = colors.HexColor("#008C95")
GREEN = colors.HexColor("#2E7D32")
AMBER = colors.HexColor("#B26A00")
RED = colors.HexColor("#A33A46")
INK = colors.HexColor("#15202B")
MUTED = colors.HexColor("#536273")
LIGHT_BLUE = colors.HexColor("#EAF3FA")
LIGHT_GREEN = colors.HexColor("#EAF5EC")
LIGHT_AMBER = colors.HexColor("#FFF3DC")
LIGHT_GRAY = colors.HexColor("#F3F5F7")
WHITE = colors.white


def register_fonts() -> None:
    font_dir = Path(r"C:\Windows\Fonts")
    pdfmetrics.registerFont(TTFont("Arial", font_dir / "arial.ttf"))
    pdfmetrics.registerFont(TTFont("Arial-Bold", font_dir / "arialbd.ttf"))
    pdfmetrics.registerFont(TTFont("Arial-Italic", font_dir / "ariali.ttf"))
    pdfmetrics.registerFont(TTFont("Arial-BoldItalic", font_dir / "arialbi.ttf"))


def styles() -> dict[str, ParagraphStyle]:
    return {
        "title": ParagraphStyle(
            "title", fontName="Arial-Bold", fontSize=26, leading=30, textColor=WHITE
        ),
        "subtitle": ParagraphStyle(
            "subtitle", fontName="Arial", fontSize=12, leading=16, textColor=colors.HexColor("#DDEAF4")
        ),
        "h1": ParagraphStyle(
            "h1", fontName="Arial-Bold", fontSize=17, leading=20, textColor=NAVY
        ),
        "h2": ParagraphStyle(
            "h2", fontName="Arial-Bold", fontSize=11.5, leading=14, textColor=NAVY
        ),
        "body": ParagraphStyle(
            "body", fontName="Arial", fontSize=8.8, leading=11.3, textColor=INK
        ),
        "small": ParagraphStyle(
            "small", fontName="Arial", fontSize=7.3, leading=9.3, textColor=INK
        ),
        "tiny": ParagraphStyle(
            "tiny", fontName="Arial", fontSize=6.3, leading=7.8, textColor=INK
        ),
        "center": ParagraphStyle(
            "center", fontName="Arial", fontSize=8.5, leading=11, textColor=INK, alignment=TA_CENTER
        ),
        "equation": ParagraphStyle(
            "equation", fontName="Arial", fontSize=11, leading=15, textColor=NAVY, alignment=TA_CENTER
        ),
        "box_title": ParagraphStyle(
            "box_title", fontName="Arial-Bold", fontSize=9.2, leading=11, textColor=NAVY
        ),
    }


def paragraph(c: canvas.Canvas, text: str, style: ParagraphStyle, x: float, top: float, width: float) -> float:
    item = Paragraph(text, style)
    _, height = item.wrap(width, PAGE_H)
    item.drawOn(c, x, top - height)
    return height


def box(
    c: canvas.Canvas,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    fill: colors.Color,
    stroke: colors.Color,
    title: str,
    body: str,
    st: dict[str, ParagraphStyle],
    title_color: colors.Color | None = None,
) -> None:
    c.setFillColor(fill)
    c.setStrokeColor(stroke)
    c.setLineWidth(0.8)
    c.roundRect(x, y, width, height, 5, fill=1, stroke=1)
    c.setFillColor(stroke)
    c.roundRect(x, y, 4, height, 2, fill=1, stroke=0)
    title_style = st["box_title"].clone("box_title_local")
    title_style.textColor = title_color or stroke
    title_height = paragraph(c, title, title_style, x + 11, y + height - 8, width - 20)
    paragraph(
        c,
        body,
        st["small"],
        x + 11,
        y + height - 11 - title_height,
        width - 20,
    )


def badge(c: canvas.Canvas, x: float, y: float, label: str, color: colors.Color) -> None:
    width = pdfmetrics.stringWidth(label, "Arial-Bold", 7.2) + 12
    c.setFillColor(color)
    c.roundRect(x, y, width, 14, 7, fill=1, stroke=0)
    c.setFillColor(WHITE)
    c.setFont("Arial-Bold", 7.2)
    c.drawCentredString(x + width / 2, y + 4.1, label)


def header(c: canvas.Canvas, page: int, title: str, section: str) -> None:
    c.setFillColor(NAVY)
    c.rect(0, PAGE_H - 43, PAGE_W, 43, fill=1, stroke=0)
    c.setFillColor(WHITE)
    c.setFont("Arial-Bold", 12.5)
    c.drawString(MARGIN, PAGE_H - 27, title)
    c.setFont("Arial", 7.5)
    c.setFillColor(colors.HexColor("#D8E5EF"))
    c.drawRightString(PAGE_W - MARGIN, PAGE_H - 25.5, section)
    c.setStrokeColor(colors.HexColor("#C8D1D9"))
    c.line(MARGIN, 21, PAGE_W - MARGIN, 21)
    c.setFillColor(MUTED)
    c.setFont("Arial", 6.8)
    c.drawString(MARGIN, 10, "pySNSPD · D3 · corrida fotónica completa · 10-08-2026")
    c.drawRightString(PAGE_W - MARGIN, 10, f"Página {page} de 10")


def crop_image(source: Path) -> Path:
    target_dir = WORK / "cropped"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / source.name.replace("page-", source.parent.name + "-")
    if target.exists() and target.stat().st_mtime_ns >= source.stat().st_mtime_ns:
        return target
    image = PILImage.open(source).convert("RGB")
    arr = np.asarray(image)
    mask = np.any(arr < 247, axis=2)
    rows, cols = np.where(mask)
    if rows.size:
        pad = 12
        left = max(0, int(cols.min()) - pad)
        right = min(image.width, int(cols.max()) + pad + 1)
        top = max(0, int(rows.min()) - pad)
        bottom = min(image.height, int(rows.max()) + pad + 1)
        image = image.crop((left, top, right, bottom))
    image.save(target, optimize=True)
    return target


def draw_image_fit(
    c: canvas.Canvas,
    source: Path,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    border: bool = True,
) -> None:
    path = crop_image(source)
    image = PILImage.open(path)
    scale = min(width / image.width, height / image.height)
    draw_w = image.width * scale
    draw_h = image.height * scale
    draw_x = x + 0.5 * (width - draw_w)
    draw_y = y + 0.5 * (height - draw_h)
    if border:
        c.setFillColor(WHITE)
        c.setStrokeColor(colors.HexColor("#CFD7DE"))
        c.roundRect(x, y, width, height, 4, fill=1, stroke=1)
    c.drawImage(str(path), draw_x, draw_y, draw_w, draw_h, preserveAspectRatio=True, mask="auto")


def bullet_list(
    c: canvas.Canvas,
    items: Iterable[str],
    x: float,
    top: float,
    width: float,
    st: dict[str, ParagraphStyle],
    *,
    style: str = "small",
    gap: float = 3,
) -> float:
    cursor = top
    for item in items:
        text = f"<bullet>&bull;</bullet>{item}"
        pstyle = st[style].clone(f"bullet_{cursor}")
        pstyle.leftIndent = 10
        pstyle.firstLineIndent = -7
        used = paragraph(c, text, pstyle, x, cursor, width)
        cursor -= used + gap
    return cursor


def build_metrics() -> dict[str, float]:
    with np.load(CACHE) as z:
        time_ps = z["time_ps"]
        dt_s = z["dt_ps"] * 1.0e-12
        photon_ps = float(z["photon_time_ps"][0])
        post = time_ps >= photon_ps
        post_corr = time_ps >= photon_ps + 0.2
        eV = 1.602176634e-19
        metrics: dict[str, float] = {
            "photon_ps": photon_ps,
            "vmax_ps": float(z["vmax_time_ps"][0]),
            "final_ps": float(z["final_time_ps"][0]),
            "n_intervals": float(time_ps.size),
            "duplicate_times": float(z["dropped_duplicate_time_count"][0]),
        }
        for short in (
            "P_spec",
            "P_delta",
            "P_delta_cond",
            "P_q",
            "P_path",
            "P_J",
            "minus_P_ep",
            "P_diff",
            "Q_ret",
            "residual",
            "P_esc",
        ):
            values = z[f"integrated_{short}_W"]
            metrics[f"E_{short}_signed_eV"] = float(np.sum(values[post] * dt_s[post]) / eV)
            metrics[f"E_{short}_abs_eV"] = float(np.sum(np.abs(values[post]) * dt_s[post]) / eV)
            idx = np.flatnonzero(post)[int(np.argmax(np.abs(values[post])))]
            metrics[f"peak_{short}_W"] = float(values[idx])
            metrics[f"peak_{short}_time_ps"] = float(time_ps[idx])
        for short in ("P_spec", "P_delta", "P_delta_cond", "P_q", "P_path", "Q_ret", "residual"):
            values = z[f"p99_abs_{short}_W_m3"]
            idx = np.flatnonzero(post)[int(np.argmax(values[post]))]
            metrics[f"p99_{short}"] = float(values[idx])
            metrics[f"p99_{short}_time_ps"] = float(time_ps[idx])
        spec = z["integrated_P_spec_W"][post_corr]
        residual = z["integrated_residual_W"][post_corr]
        metrics["residual_spec_corr"] = float(np.corrcoef(spec, residual)[0, 1])
        metrics["residual_spec_rel_rms"] = float(
            np.sqrt(np.mean((residual - spec) ** 2)) / np.sqrt(np.mean(spec**2))
        )
        for key in ("catalog_clipped_fraction", "catalog_Te_clipped_fraction", "catalog_Tph_clipped_fraction", "catalog_delta_clipped_fraction", "catalog_q_clipped_fraction"):
            values = z[key]
            metrics[f"{key}_max"] = float(np.max(values[post]))
        q_clip = z["catalog_q_clipped_fraction"]
        metrics["q_clip_duration_ps"] = float(np.sum(dt_s[post][q_clip[post] > 0]) * 1.0e12)
        active = np.flatnonzero(post & (q_clip > 0))
        metrics["q_clip_last_ps"] = float(time_ps[active[-1]])

        mask = z["central_mask"]
        weights = z["node_area_m2"][mask]
        for label, position in (("peak", 0), ("mid", 1), ("final", 2)):
            metrics[f"{label}_time_ps"] = float(z["selected_times_ps"][position])
            delta = z["selected_delta_over_delta0"][position, mask]
            Te = z["selected_Te_K"][position, mask]
            Tph = z["selected_Tph_K"][position, mask]
            qxi = z["selected_q_xi"][position, mask]
            metrics[f"{label}_delta_min"] = float(np.min(delta))
            metrics[f"{label}_delta_mean"] = float(np.average(delta, weights=weights))
            metrics[f"{label}_Te_max"] = float(np.max(Te))
            metrics[f"{label}_Te_mean"] = float(np.average(Te, weights=weights))
            metrics[f"{label}_Tph_max"] = float(np.max(Tph))
            metrics[f"{label}_qxi_p99"] = float(np.percentile(qxi, 99))
            metrics[f"{label}_clip"] = float(np.mean(z["selected_catalog_clipped"][position, mask]))
        for field in ("delta_over_delta0", "Te_K", "Tph_K", "q_xi", "u_e_J_m3"):
            initial = z[f"initial_{field}"][mask]
            final = z[f"selected_{field}"][-1, mask]
            metrics[f"final_minus_initial_{field}"] = float(np.average(final - initial, weights=weights))
        return metrics


def draw_table(c: canvas.Canvas, data: list[list[str]], x: float, y: float, width: float, row_heights: list[float] | None = None) -> float:
    col_widths = [0.42 * width, 0.25 * width, 0.33 * width]
    table = Table(data, colWidths=col_widths, rowHeights=row_heights)
    table.setStyle(
        TableStyle(
            [
                ("FONTNAME", (0, 0), (-1, 0), "Arial-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Arial"),
                ("FONTSIZE", (0, 0), (-1, -1), 7.0),
                ("LEADING", (0, 0), (-1, -1), 8.5),
                ("BACKGROUND", (0, 0), (-1, 0), NAVY),
                ("TEXTCOLOR", (0, 0), (-1, 0), WHITE),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [WHITE, LIGHT_GRAY]),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#C6CED6")),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    _, height = table.wrap(width, PAGE_H)
    table.drawOn(c, x, y - height)
    return height


def page_one(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    c.setFillColor(NAVY)
    c.rect(0, 0, PAGE_W, PAGE_H, fill=1, stroke=0)
    c.setFillColor(TEAL)
    c.rect(0, 0, 12, PAGE_H, fill=1, stroke=0)
    paragraph(c, "Informe D3: proyección energética de una corrida fotónica", st["title"], 52, PAGE_H - 60, 690)
    paragraph(
        c,
        "Diagnóstico de <i>P</i><sub>Δ</sub>, <i>P</i><sub>q</sub> y de la conveniencia de evolucionar energía electrónica en pySNSPD",
        st["subtitle"],
        54,
        PAGE_H - 132,
        690,
    )
    c.setFillColor(colors.HexColor("#DDEAF4"))
    c.setFont("Arial", 8.5)
    c.drawString(55, PAGE_H - 169, "Corrida: photon_phasecg_I30uA_0p8eV_sigma10nm_t50ps_1500ps_stiffness_01")
    c.drawString(55, PAGE_H - 184, "Ventana diagnosticada: 100 nm centrales · 15 001 snapshots · 1 499.9 ps")

    c.setFillColor(WHITE)
    c.roundRect(50, 292, 742, 92, 7, fill=1, stroke=0)
    badge(c, 68, 347, "OBJETIVO ELEMENTAL", TEAL)
    paragraph(
        c,
        "Determinar si la simulación actual omite la energía asociada a los cambios de <b>Δ</b> y <b>q</b>, y decidir entre dos correcciones: añadir <i>P</i><sub>Δ</sub>/<i>P</i><sub>q</sub> o evolucionar directamente <i>u</i><sub>e</sub>.",
        st["h1"],
        68,
        337,
        690,
    )

    widths = [230, 230, 230]
    xs = [50, 306, 562]
    box(
        c, xs[0], 136, widths[0], 126,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="HECHO 1 · El término omitido es grande durante el transitorio",
        body=f"El pico p99 de |P<sub>spec</sub>| es {m['p99_P_spec']:.2e} W m<super>-3</super>. En el instante de ese pico, su potencia central integrada es {1e9*m['peak_P_spec_W']:.2f} nW.",
        st=st,
    )
    box(
        c, xs[1], 136, widths[1], 126,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="HECHO 2 · PΔ domina sobre Pq en esta trayectoria",
        body=f"La variación absoluta temporal de potencia central equivale a {m['E_P_delta_abs_eV']:.2f} eV para P<sub>Δ</sub> y {m['E_P_q_abs_eV']:.3f} eV para P<sub>q</sub>. La comparación es diagnóstica, no una nueva fuente de energía.",
        st=st,
    )
    box(
        c, xs[2], 136, widths[2], 126,
        fill=LIGHT_AMBER, stroke=AMBER,
        title="AMBIGÜEDAD CLAVE · Pq aún no es activable",
        body=f"El catálogo recorta q hasta en {100*m['catalog_q_clipped_fraction_max']:.1f}% de los nodos centrales durante {m['q_clip_duration_ps']:.1f} ps. Además, u<sub>e</sub>(q) no incluye por sí sola toda la energía del circuito y del superflujo.",
        st=st,
    )

    c.setFillColor(colors.HexColor("#DDEAF4"))
    c.setFont("Arial", 7.5)
    c.drawString(55, 92, "Convención de certeza:")
    badge(c, 150, 86, "HECHO", GREEN)
    c.setFillColor(colors.HexColor("#DDEAF4"))
    c.drawString(210, 91, "sale directamente de datos guardados o de una identidad algebraica")
    badge(c, 470, 86, "INFERENCIA", BLUE)
    c.setFillColor(colors.HexColor("#DDEAF4"))
    c.drawString(550, 91, "requiere aceptar el cierre térmico y el catálogo")
    badge(c, 150, 62, "AMBIGÜEDAD", AMBER)
    c.setFillColor(colors.HexColor("#DDEAF4"))
    c.drawString(239, 67, "no puede resolverse con esta corrida por sí sola")
    c.setFont("Arial-Italic", 7)
    c.drawString(55, 32, "“100% cierto” en este informe significa cierto dentro de los archivos persistidos, las definiciones D3 y la precisión numérica mostrada; no una ley universal del modelo reducido.")
    c.showPage()


def page_two(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    header(c, 2, "Qué problema se está diagnosticando", "Objetivo, variables y ecuaciones")
    paragraph(c, "1. Una contabilidad energética sencilla", st["h1"], MARGIN, PAGE_H - 62, 420)
    paragraph(
        c,
        "El subsistema electrónico no se describe solo por una temperatura. También cambia el <b>gap superconductivo Δ</b> y el <b>momento de superflujo q</b>. Por eso una misma temperatura puede corresponder a energías electrónicas distintas. La pregunta es qué parte de la potencia cambia T<sub>e</sub> y qué parte mueve el fondo superconductivo.",
        st["body"],
        MARGIN,
        PAGE_H - 90,
        420,
    )
    box(
        c, MARGIN, 328, 420, 85,
        fill=LIGHT_BLUE, stroke=BLUE,
        title="Analogía para no especialistas",
        body="La temperatura es el dinero disponible en una cuenta corriente; Δ y q describen activos cuyo valor también cambia. Mirar solo el saldo de la cuenta corriente no reconstruye el patrimonio total.",
        st=st,
    )
    paragraph(c, "Variables que aparecen en las figuras", st["h2"], MARGIN, 307, 420)
    variable_data = [
        ["Símbolo", "Nombre", "Papel en el balance"],
        ["Tₑ", "temperatura electrónica", "energía térmica de excitaciones"],
        ["Tph", "temperatura fonónica", "reservorio de fonones"],
        ["Δ", "gap / condensado", "estado del fondo superconductivo"],
        ["q", "momento de superflujo", "depairing y redistribución de corriente"],
        ["uₑ", "energía electrónica", "uₑ = uqp + ucond"],
        ["Qret", "potencia retenida", "PJ + Pdiff − Pe-ph"],
    ]
    draw_table(c, variable_data, MARGIN, 286, 420)

    right_x = 486
    paragraph(c, "2. Ecuaciones usadas por D3", st["h1"], right_x, PAGE_H - 62, 320)
    equations = [
        "<b>Energía de estado:</b><br/>u<sub>e</sub> = u<sub>qp</sub> + u<sub>cond</sub>",
        "<b>Regla de la cadena:</b><br/>du<sub>e</sub>/dt = C<sub>e</sub> dT<sub>e</sub>/dt + P<sub>Δ</sub> + P<sub>q</sub>",
        "<b>Solver actual:</b><br/>C<sub>e</sub> dT<sub>e</sub>/dt = Q<sub>ret</sub> = P<sub>J</sub> + P<sub>diff</sub> − P<sub>e-ph</sub>",
        "<b>Diagnóstico robusto:</b><br/>P<sub>spec</sub> = [u<sub>e</sub>(T<sub>e</sub><super>n</super>, Δ<super>n+1</super>, q<super>n+1</super>) − u<sub>e</sub>(T<sub>e</sub><super>n</super>, Δ<super>n</super>, q<super>n</super>)] / Δt",
        "<b>Residuo:</b><br/>R<sub>e</sub> = Δu<sub>e</sub>/Δt − Q<sub>ret</sub>",
    ]
    y = 474
    heights = [58, 67, 67, 82, 58]
    for text, height in zip(equations, heights):
        c.setFillColor(LIGHT_GRAY)
        c.setStrokeColor(colors.HexColor("#C9D2DA"))
        c.roundRect(right_x, y - height, 320, height - 7, 5, fill=1, stroke=1)
        paragraph(c, text, st["equation"], right_x + 10, y - 12, 300)
        y -= height
    box(
        c, right_x, 55, 320, 88,
        fill=LIGHT_AMBER, stroke=AMBER,
        title="Qué es robusto y qué depende del camino",
        body="P<sub>spec</sub> es la diferencia total de una función de estado. Su separación finita en P<sub>Δ</sub> y P<sub>q</sub> depende de si se cambia primero Δ o primero q; D3 muestra esa ambigüedad como P<sub>path</sub>.",
        st=st,
    )
    c.showPage()


def page_three(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    header(c, 3, "Evolución temporal completa", "Escalas locales, potencias y energía integrada")
    draw_image_fit(c, ASSETS / "temporal" / "page-1.png", MARGIN, 52, 525, 480)
    x = 570
    paragraph(c, "Lectura cuantitativa", st["h1"], x, PAGE_H - 66, 235)
    data = [
        ["Magnitud", "Resultado", "Instante / alcance"],
        ["p99 |Pspec|", f"{m['p99_P_spec']:.2e}", f"{m['p99_P_spec_time_ps']:.2f} ps"],
        ["p99 |PΔ|", f"{m['p99_P_delta']:.2e}", f"{m['p99_P_delta_time_ps']:.2f} ps"],
        ["p99 |Pq|", f"{m['p99_P_q']:.2e}", f"{m['p99_P_q_time_ps']:.2f} ps"],
        ["∫|Pspec,central|dt", f"{m['E_P_spec_abs_eV']:.3f} eV", "desde 50 ps"],
        ["∫Pspec,central dt", f"{m['E_P_spec_signed_eV']:.3f} eV", "casi reversible"],
        ["corr(Re, Pspec)", f"{m['residual_spec_corr']:.5f}", "t ≥ 50.2 ps"],
        ["error RMS relativo", f"{100*m['residual_spec_rel_rms']:.2f}%", "Re frente a Pspec"],
    ]
    draw_table(c, data, x, PAGE_H - 94, 235)
    box(
        c, x, 160, 235, 102,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="HECHO · El término es transitoriamente relevante",
        body=f"La energía central firmada de P<sub>spec</sub> termina cerca de cero ({m['E_P_spec_signed_eV']:.3f} eV), pero su recorrido absoluto es {m['E_P_spec_abs_eV']:.2f} eV. Hay almacenamiento y devolución; no ausencia de dinámica.",
        st=st,
    )
    box(
        c, x, 52, 235, 96,
        fill=LIGHT_BLUE, stroke=BLUE,
        title="INFERENCIA · El residuo identifica la omisión",
        body="La casi identidad temporal entre R<sub>e</sub> y P<sub>spec</sub> indica que el cierre en T<sub>e</sub> satisface el resto del balance y deja fuera el movimiento espectral. Es evidencia fuerte en esta corrida, no una prueba universal.",
        st=st,
    )
    c.showPage()


def page_figure_with_notes(
    c: canvas.Canvas,
    st: dict[str, ParagraphStyle],
    *,
    page: int,
    title: str,
    section: str,
    image: Path,
    notes: list[tuple[str, str, colors.Color, colors.Color]],
) -> None:
    header(c, page, title, section)
    draw_image_fit(c, image, MARGIN, 142, PAGE_W - 2 * MARGIN, 390)
    count = len(notes)
    gap = 12
    width = (PAGE_W - 2 * MARGIN - gap * (count - 1)) / count
    for i, (note_title, body, stroke, fill) in enumerate(notes):
        box(
            c,
            MARGIN + i * (width + gap),
            45,
            width,
            82,
            fill=fill,
            stroke=stroke,
            title=note_title,
            body=body,
            st=st,
        )
    c.showPage()


def page_seven(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    header(c, 7, "Potencias retenidas y cierre electrónico", "Comparación directa entre el solver actual y el libro mayor")
    draw_image_fit(c, ASSETS / "colormaps" / "page-4.png", MARGIN, 292, 385, 240)
    draw_image_fit(c, ASSETS / "colormaps" / "page-5.png", 435, 292, 372, 240)
    paragraph(c, "Potencias presentes en la ecuación de Tₑ", st["h2"], MARGIN + 8, 279, 360)
    paragraph(c, "Derivada total, residuo y sensibilidad al camino", st["h2"], 443, 279, 350)
    data = [
        ["Canal post-fotón", "Energía firmada", "Lectura"],
        ["Joule", f"+{m['E_P_J_signed_eV']:.2f} eV", "calienta electrones"],
        ["Difusión", f"{m['E_P_diff_signed_eV']:.2f} eV", "extrae energía"],
        ["−e-ph", f"+{m['E_minus_P_ep_signed_eV']:.3f} eV", "neto firmado"],
        ["Qret total", f"+{m['E_Q_ret_signed_eV']:.3f} eV", "gran cancelación interna"],
        ["Residuo Re", f"{m['E_residual_signed_eV']:.3f} eV", "sigue a Pspec"],
    ]
    draw_table(c, data, MARGIN, 250, 385)
    box(
        c, 435, 142, 372, 108,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="HECHO · La cancelación de canales es grande",
        body=f"Desde el fotón, P<sub>J</sub> aporta {m['E_P_J_signed_eV']:.2f} eV y la difusión retira {abs(m['E_P_diff_signed_eV']):.2f} eV. Q<sub>ret</sub> queda en {m['E_Q_ret_signed_eV']:.3f} eV. Por eso una omisión menor que cada canal puede ser grande frente al balance neto.",
        st=st,
    )
    box(
        c, MARGIN, 45, 385, 85,
        fill=LIGHT_BLUE, stroke=BLUE,
        title="INFERENCIA · Evolucionar uₑ es numéricamente atractivo",
        body="La formulación conservativa incorpora P<sub>spec</sub> por construcción y evita restar canales grandes para reconstruir una corrección pequeña. Aun así necesita inversión u<sub>e</sub>→T<sub>e</sub> bien condicionada.",
        st=st,
    )
    box(
        c, 435, 45, 372, 85,
        fill=LIGHT_AMBER, stroke=AMBER,
        title="AMBIGÜEDAD · Posible doble conteo con gTDGL/circuito",
        body="El residuo muestra energía electrónica omitida según el funcional de Simon. No determina por sí solo dónde debe ir la disipación gTDGL ni qué parte de la energía asociada a q pertenece al circuito.",
        st=st,
    )
    c.showPage()


def page_eight(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    header(c, 8, "No linealidad térmica y cobertura del catálogo", "Dónde la evaluación es fiable y dónde queda limitada")
    draw_image_fit(c, ASSETS / "colormaps" / "page-6.png", MARGIN, 150, PAGE_W - 2 * MARGIN, 382)
    widths = [240, 240, 240]
    xs = [MARGIN, MARGIN + 255, MARGIN + 510]
    box(
        c, xs[0], 45, widths[0], 88,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="HECHO · El clipping pertenece solo a q",
        body=f"Máximos fuera de catálogo: q={100*m['catalog_q_clipped_fraction_max']:.1f}%; T<sub>e</sub>=0%; T<sub>ph</sub>=0%; Δ=0%. El clipping de q dura {m['q_clip_duration_ps']:.1f} ps y termina a {m['q_clip_last_ps']:.2f} ps.",
        st=st,
    )
    box(
        c, xs[1], 45, widths[1], 88,
        fill=LIGHT_AMBER, stroke=AMBER,
        title="AMBIGÜEDAD · Pq queda subdeterminado",
        body="Durante el clipping, D3 usa el borde del catálogo: no extrapola. El valor mostrado de P<sub>q</sub> es una evaluación saturada y no una medida cuantitativa completa del régimen de q alto.",
        st=st,
    )
    box(
        c, xs[2], 45, widths[2], 88,
        fill=LIGHT_BLUE, stroke=BLUE,
        title="INFERENCIA · Prioridad técnica inmediata",
        body="Extender/refinar el eje q y repetir D3 es obligatorio antes de activar P<sub>q</sub>. La parte explícita de condensación de P<sub>Δ</sub> no depende de q, pero su parte QP/DOS sí puede heredar este límite.",
        st=st,
    )
    c.showPage()


def page_nine(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    header(c, 9, "Perfiles longitudinales del estado", "Promedio transversal en los 100 nm centrales")
    draw_image_fit(c, ASSETS / "profiles" / "page-1.png", MARGIN, 55, 535, 475)
    x = 585
    paragraph(c, "Tres tiempos de decisión", st["h1"], x, PAGE_H - 66, 220)
    data = [
        ["Tiempo", "Estado central", "Valores"],
        [f"{m['peak_time_ps']:.1f} ps", "máximo VTDGL", f"<Δ>/Δ₀={m['peak_delta_mean']:.3f}; max Te={m['peak_Te_max']:.2f} K"],
        [f"{m['mid_time_ps']:.1f} ps", "recuperación", f"<Δ>/Δ₀={m['mid_delta_mean']:.3f}; <Te>={m['mid_Te_mean']:.3f} K"],
        ["1500 ps", "final", f"<Δ>/Δ₀={m['final_delta_mean']:.3f}; <Te>={m['final_Te_mean']:.3f} K"],
    ]
    # Plain text table avoids interpreting angle brackets as markup.
    data[1][2] = f"mean Δ/Δ0={m['peak_delta_mean']:.3f}; max Te={m['peak_Te_max']:.2f} K"
    data[2][2] = f"mean Δ/Δ0={m['mid_delta_mean']:.3f}; mean Te={m['mid_Te_mean']:.3f} K"
    data[3][2] = f"mean Δ/Δ0={m['final_delta_mean']:.3f}; mean Te={m['final_Te_mean']:.3f} K"
    draw_table(c, data, x, PAGE_H - 98, 220)
    box(
        c, x, 176, 220, 102,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="HECHO · Recuperación fuerte, no identidad exacta",
        body=f"Al final, el promedio central de Δ/Δ₀ difiere del inicial en {m['final_minus_initial_delta_over_delta0']:+.4f}; T<sub>e</sub> difiere en {m['final_minus_initial_Te_K']:+.4f} K y qξ en {m['final_minus_initial_q_xi']:+.4f}.",
        st=st,
    )
    box(
        c, x, 55, 220, 108,
        fill=LIGHT_AMBER, stroke=AMBER,
        title="AMBIGÜEDAD · Recuperación del detector ≠ cierre energético",
        body="Que los perfiles vuelvan cerca del estado inicial no prueba conservación durante el trayecto. Una excursión grande puede devolver casi toda la energía y aun alterar latencia, phase slips o umbral de detección.",
        st=st,
    )
    c.showPage()


def page_ten(c: canvas.Canvas, st: dict[str, ParagraphStyle], m: dict[str, float]) -> None:
    header(c, 10, "Conclusiones y decisión recomendada", "Perfiles energéticos, certezas y próximos pasos")
    draw_image_fit(c, ASSETS / "profiles" / "page-2.png", MARGIN, 276, 380, 255)
    draw_image_fit(c, ASSETS / "profiles" / "page-3.png", 428, 276, 379, 255)
    paragraph(c, "Cambios de energía respecto del estado inicial", st["h2"], MARGIN + 5, 267, 360)
    paragraph(c, "Forma espacial normalizada de PΔ, Pq, Pspec y Re", st["h2"], 433, 267, 365)

    box(
        c, MARGIN, 151, 245, 98,
        fill=LIGHT_GREEN, stroke=GREEN,
        title="CONCLUSIÓN CIERTA EN ESTA CORRIDA",
        body="P<sub>spec</sub> no es despreciable durante el evento; P<sub>Δ</sub> es su contribución dominante; y el solver actual no registra explícitamente ese almacenamiento al evolucionar T<sub>e</sub>.",
        st=st,
    )
    box(
        c, MARGIN + 258, 151, 245, 98,
        fill=LIGHT_BLUE, stroke=BLUE,
        title="DECISIÓN TÉCNICA MEJOR SUSTENTADA",
        body="Prototipar una actualización conservativa en u<sub>e</sub> antes que añadir P<sub>Δ</sub> y P<sub>q</sub> como fuentes independientes. Compararla contra el solver actual sin cambiar el resto del splitting.",
        st=st,
    )
    box(
        c, MARGIN + 516, 151, 245, 98,
        fill=LIGHT_AMBER, stroke=AMBER,
        title="NO DECIDIDO TODAVÍA",
        body="No puede activarse P<sub>q</sub> cuantitativamente hasta ampliar el catálogo y fijar la partición de energía entre DOS electrónica, superflujo, campo y circuito.",
        st=st,
    )

    paragraph(c, "Secuencia mínima recomendada", st["h2"], MARGIN, 135, 360)
    bullet_list(
        c,
        [
            "Extender el eje q hasta cubrir toda la trayectoria y repetir D3 con una prueba de convergencia temporal.",
            "Implementar un prototipo local que actualice u<sub>e</sub> y luego invierta T<sub>e</sub>(u<sub>e</sub>, Δ, q).",
            "Verificar ciclo cerrado, balance global y ausencia de doble conteo con gTDGL/circuito.",
            "Comparar las variables de decisión: latencia, primer phase slip, pico de V<sub>TDGL</sub>, recuperación y corriente umbral.",
        ],
        MARGIN,
        118,
        740,
        st,
        style="small",
        gap=1.5,
    )
    c.showPage()


def build() -> Path:
    register_fonts()
    st = styles()
    metrics = build_metrics()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    c = canvas.Canvas(str(OUTPUT), pagesize=(PAGE_W, PAGE_H), pageCompression=1)
    c.setTitle("Informe D3 - Proyección energética de una corrida fotónica")
    c.setAuthor("pySNSPD - informe para supervisión")
    c.setSubject("Diagnóstico de P_Delta, P_q y actualización conservativa en energía")

    page_one(c, st, metrics)
    page_two(c, st, metrics)
    page_three(c, st, metrics)
    page_figure_with_notes(
        c,
        st,
        page=4,
        title="Estado persistido usado por la proyección",
        section="Δ, q, Te y Tph en tres tiempos",
        image=ASSETS / "colormaps" / "page-1.png",
        notes=[
            (
                "HECHO · Máxima respuesta a 56.6 ps",
                f"En la ventana central: mean Δ/Δ₀={metrics['peak_delta_mean']:.3f}, min={metrics['peak_delta_min']:.1e}; max T<sub>e</sub>={metrics['peak_Te_max']:.2f} K y max T<sub>ph</sub>={metrics['peak_Tph_max']:.2f} K.",
                GREEN,
                LIGHT_GREEN,
            ),
            (
                "HECHO · El fondo casi se recupera",
                f"A 778.3 ps y 1500 ps, mean Δ/Δ₀={metrics['mid_delta_mean']:.3f} y {metrics['final_delta_mean']:.3f}; mean T<sub>e</sub>={metrics['mid_Te_mean']:.3f} y {metrics['final_Te_mean']:.3f} K.",
                GREEN,
                LIGHT_GREEN,
            ),
            (
                "LÍMITE · q alto en el evento",
                f"A 56.6 ps, p99(|q|ξ)={metrics['peak_qxi_p99']:.2f} y {100*metrics['peak_clip']:.1f}% de nodos centrales están recortados por el catálogo.",
                AMBER,
                LIGHT_AMBER,
            ),
        ],
    )
    page_figure_with_notes(
        c,
        st,
        page=5,
        title="Libro mayor electrónico y energía fonónica",
        section="ucond, uqp, ue y uph",
        image=ASSETS / "colormaps" / "page-2.png",
        notes=[
            (
                "HECHO · El colapso del gap cuesta energía",
                "u<sub>cond</sub> se vuelve menos negativo cuando Δ colapsa. Simultáneamente aumentan u<sub>qp</sub> y u<sub>ph</sub>. La suma mostrada es u<sub>e</sub>=u<sub>qp</sub>+u<sub>cond</sub>.",
                GREEN,
                LIGHT_GREEN,
            ),
            (
                "HECHO · El estado final no es idéntico",
                f"El promedio central final de u<sub>e</sub> queda {metrics['final_minus_initial_u_e_J_m3']:+.2f} J m<super>-3</super> respecto del inicial. Es pequeño frente a la excursión, pero resoluble.",
                GREEN,
                LIGHT_GREEN,
            ),
            (
                "AMBIGÜEDAD · El cero energético es parcial",
                "u<sub>e</sub> es la energía electrónica de Simon respecto del metal normal. No representa por sí sola la energía electromagnética, inductiva ni toda la energía libre gTDGL.",
                AMBER,
                LIGHT_AMBER,
            ),
        ],
    )
    page_figure_with_notes(
        c,
        st,
        page=6,
        title="Movimiento del espectro superconductivo",
        section="PΔ, condensación, QP/DOS y Pq",
        image=ASSETS / "colormaps" / "page-3.png",
        notes=[
            (
                "HECHO · Dominio de PΔ",
                f"Picos p99: |P<sub>Δ</sub>|={metrics['p99_P_delta']:.2e}, |P<sub>Δ,cond</sub>|={metrics['p99_P_delta_cond']:.2e} y |P<sub>q</sub>|={metrics['p99_P_q']:.2e} W m<super>-3</super>.",
                GREEN,
                LIGHT_GREEN,
            ),
            (
                "HECHO · El total es más robusto que la separación",
                f"La ambigüedad de camino alcanza {metrics['p99_P_path']:.2e} W m<super>-3</super>. P<sub>spec</sub> es independiente de esa elección finita; P<sub>Δ</sub> y P<sub>q</sub> por separado no.",
                GREEN,
                LIGHT_GREEN,
            ),
            (
                "AMBIGÜEDAD · Pq es espectral, no energía total de corriente",
                "El P<sub>q</sub> mostrado proviene solo de la dependencia de la DOS dentro de u<sub>e</sub>. No incluye automáticamente energía de gradiente, inductancia cinética o trabajo del circuito.",
                AMBER,
                LIGHT_AMBER,
            ),
        ],
    )
    page_seven(c, st, metrics)
    page_eight(c, st, metrics)
    page_nine(c, st, metrics)
    page_ten(c, st, metrics)
    c.save()
    return OUTPUT


if __name__ == "__main__":
    print(build())
