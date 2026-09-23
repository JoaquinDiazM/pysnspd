"""Read A20 Fig2.7(a) embedded RGB pixels; no cascade or fitted physics.

Axes were calibrated on the original 1750x1315 image from PDF page35.
Error envelopes combine visible stroke extent with +/-2px axis calibration
and +/-2px horizontal sampling. They are extraction budgets, not confidence
intervals of the physical model. The endpoint t/tau0=1e-4 is not reported:
the legend covers most of the R50 path there and the edge sliver is not used.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from statistics import median

import pdfplumber
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "tmp/pdfs/modelo_v0_3/sources/Allmaras_Thesis_Final.pdf"
OUTPUT = ROOT / "docs/implementation/stage3_5/research_20260923/cascade_digitization.json"
EXPECTED_SHA = "2554e9a2b047f93f7c7f50a2897887afaf522eb530e057307c5a31ef4e1e0e96"
AXES = {"left_x": 228., "right_x": 1582., "top_y": 125., "bottom_y": 1072.,
        "radius_min_nm": 0., "radius_max_nm": 6.,
        "time_axis_factor_tau0": 1e-4, "tau0_ps": 1870.}
COLORS = {"R50": (41, 0, 102), "R90": (204, 0, 0)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def linear_fit(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xm, ym = sum(xs)/len(xs), sum(ys)/len(ys)
    slope = sum((x-xm)*(y-ym) for x, y in points)/sum((x-xm)**2 for x in xs)
    return slope, ym-slope*xm


def extract(image, time_fraction, color, tolerance):
    target_x = AXES["left_x"] + time_fraction*(AXES["right_x"]-AXES["left_x"])
    if not 0 < time_fraction < .6:
        raise ValueError("Requested time outside the declared unobscured extraction interval")
    columns = range(round(target_x)-2, round(target_x)+3)
    pixels = []
    centers = []
    for x in columns:
        rows = [y for y in range(135, 1060)
                if max(abs(c-v) for c, v in zip(image.getpixel((x, y)), color)) <= tolerance]
        if not rows or max(rows)-min(rows)>10:
            raise ValueError(f"Missing/ambiguous curve at x={x}, color={color}")
        pixels.extend((x, y) for y in rows)
        centers.append((x, median(rows)))
    slope, intercept = linear_fit(centers)
    y = slope*target_x+intercept
    projected_rows = [yy+slope*(target_x-xx) for xx, yy in pixels]
    # Visible line plus y-axis calibration, x calibration and pixel rounding.
    pad_y = 2 + abs(slope)*2 + .5
    ylo, yhi = min(projected_rows)-pad_y, max(projected_rows)+pad_y
    radius = (AXES["bottom_y"]-y)*6/(AXES["bottom_y"]-AXES["top_y"])
    candidates = [(bottom-yy)*6/(bottom-top)
                  for bottom in (AXES["bottom_y"]-2, AXES["bottom_y"]+2)
                  for top in (AXES["top_y"]-2, AXES["top_y"]+2)
                  for yy in (ylo, yhi)]
    return {"radius_nm": radius, "extraction_envelope_nm": [min(candidates), max(candidates)],
            "target_x_px": target_x, "fitted_y_px": y, "local_slope_px_per_px": slope,
            "sampled_centers_px": centers, "stroke_pixels_count": len(pixels),
            "endpoint_extrapolation_px": 0.}


def main():
    source_hash = sha(SOURCE)
    if source_hash != EXPECTED_SHA:
        raise ValueError("Source changed: axes and colors require fresh inspection")
    with pdfplumber.open(SOURCE) as pdf:
        page = pdf.pages[34]
        objects = page.images
        obj = objects[0]
        if obj["name"] != "Im13" or tuple(obj["srcsize"]) != (1750, 1315):
            raise ValueError("Unexpected embedded figure")
        raw = obj["stream"].get_data()
        image = Image.frombytes("RGB", tuple(obj["srcsize"]), raw)
        image_hash = hashlib.sha256(raw).hexdigest()
        vector_counts = {"curves": len(page.curves), "lines": len(page.lines)}
    gaussian_ratio = math.sqrt(math.log(10)/math.log(2))
    rows = []
    maximum_threshold_change = 0.
    for fraction in (.3, .5):
        measures = {key: extract(image, fraction, color, 12)
                    for key, color in COLORS.items()}
        for key, color in COLORS.items():
            for tolerance in (4, 20):
                alternate = extract(image, fraction, color, tolerance)
                maximum_threshold_change = max(maximum_threshold_change,
                    abs(alternate["radius_nm"]-measures[key]["radius_nm"]))
        lo50, hi50 = measures["R50"]["extraction_envelope_nm"]
        lo90, hi90 = measures["R90"]["extraction_envelope_nm"]
        ratio_bounds = [lo90/hi50, hi90/lo50]
        r50, r90 = measures["R50"]["radius_nm"], measures["R90"]["radius_nm"]
        rows.append({"time_axis_fraction": fraction, "time_over_tau0": fraction*1e-4,
                     "time_ps": fraction*1e-4*1870, **measures,
                     "R90_over_R50": r90/r50, "ratio_extraction_envelope": ratio_bounds,
                     "gaussian_s_from_R50_nm": r50/math.sqrt(2*math.log(2)),
                     "gaussian_s_from_R90_nm": r90/math.sqrt(2*math.log(10)),
                     "single_2d_gaussian_ratio_outside_extraction_envelope":
                         not ratio_bounds[0] <= gaussian_ratio <= ratio_bounds[1]})
    result = {
        "schema": "pysnspd.stage3_5.A20_radial_figure_digitization.v1",
        "status": "HISTORICAL_CURVE_EXTRACTION_NOT_KORZH_DATA",
        "new_physics_simulations": 0,
        "source": {"path": str(SOURCE.relative_to(ROOT)).replace("\\", "/"),
                   "sha256": source_hash,
                   "url": "https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf",
                   "PDF_page_1based": 35, "printed_page": 23, "figure": "2.7(a)",
                   "subsystem": "total photon energy", "geometry": "cylindrical",
                   "embedded_image": "Im13", "image_size_px": [1750, 1315],
                   "decoded_RGB_sha256": image_hash, "vector_objects": vector_counts},
        "script_sha256": sha(Path(__file__)),
        "axes": AXES,
        "extraction": {"method": "RGB embedded-image curve pixels plus local linear interpolation",
                       "nominal_RGB_Linf_tolerance": 12, "crosscheck_tolerances": [4, 20],
                       "maximum_threshold_change_nm": maximum_threshold_change,
                       "axis_coordinate_uncertainty_px": 2,
                       "envelope_kind": "conservative reading budget, not statistical or physical interval",
                       "no_curve_hidden_by_legend_used": True,
                       "no_raster_resizing_or_smoothing": True},
        "two_dimensional_gaussian_R90_over_R50": gaussian_ratio,
        "rows": rows,
        "not_extracted": [{"time_over_tau0": 1e-4, "R50_nm": None,
                           "R90_over_R50": None,
                           "reason": "Legend obscures R50 over most late times; tiny visible edge segment not used or extrapolated."}],
        "interpretation": "A single 2D Gaussian cannot match both displayed total-energy percentiles within the stated extraction budget; this neither determines a phonon profile nor rejects a future controlled Gaussian reduction.",
        "not_identified": ["phonon-only percentiles", "total second moment",
                           "central energy density", "Korzh775/1550nm0.9K profiles"],
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT.relative_to(ROOT)),
                      "rows": [{k: row[k] for k in ("time_ps", "R90_over_R50",
                               "ratio_extraction_envelope")} for row in rows],
                      "maximum_threshold_change_nm": maximum_threshold_change},
                     indent=2))


if __name__ == "__main__":
    main()
