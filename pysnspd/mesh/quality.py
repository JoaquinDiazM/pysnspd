"""Mesh-quality diagnostics for pySNSPD finite-volume runs.

The OE7 gTDGL/Poisson solver is sensitive to local finite-volume stiffness
factors such as ``s_ij/(a_i e_ij)``.  This module provides a lightweight audit
that can be run after a PRE-run mesh has been generated, before spending time
on Usadel/phase-space catalogues or stationary dynamics.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import math
import numpy as np
import yaml


@dataclass(frozen=True)
class MeshQualityReport:
    """Serializable mesh-quality report."""

    status: str
    metrics: dict[str, float | int | str]
    warnings: list[str]
    recommendations: list[str]

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "metrics": dict(self.metrics),
            "warnings": list(self.warnings),
            "recommendations": list(self.recommendations),
        }


def build_mesh_quality_report(
    mesh,
    ops,
    *,
    min_angle_warn_deg: float = 25.0,
    min_angle_fail_deg: float = 15.0,
    max_edge_ratio_warn: float = 2.25,
    max_edge_ratio_fail: float = 3.50,
    max_fv_weight_ratio_warn: float = 40.0,
    max_fv_weight_ratio_fail: float = 100.0,
    min_node_area_ratio_warn: float = 0.20,
    min_node_area_ratio_fail: float = 0.08,
    max_dual_over_edge_warn: float = 3.0,
    max_dual_over_edge_fail: float = 8.0,
) -> MeshQualityReport:
    """Audit Delaunay/Voronoi finite-volume mesh quality.

    The most important OE7 stiffness proxy is

        s_ij / (a_i e_ij)

    where ``s_ij`` is the Voronoi dual length, ``a_i`` the node control-volume
    area, and ``e_ij`` the primary Delaunay edge length.  Large isolated values
    mean a tiny phase error can become a large local Laplacian/divergence term.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)[:, :2]
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    edges = np.asarray(ops.edges, dtype=np.int64)
    edge_len = np.asarray(ops.edge_length_m, dtype=float)
    dual = np.asarray(ops.dual_face_length_m, dtype=float)
    node_area = np.asarray(ops.node_area_m2, dtype=float)

    tri_area = _triangle_areas(nodes, triangles)
    tri_min_angle = _triangle_min_angles_deg(nodes, triangles)

    target_h = float(getattr(mesh, "target_spacing_m", np.nanmedian(edge_len)))
    median_edge = float(np.nanmedian(edge_len))
    min_edge = float(np.nanmin(edge_len))
    max_edge = float(np.nanmax(edge_len))
    edge_ratio = max_edge / max(min_edge, 1.0e-300)

    median_area = float(np.nanmedian(node_area))
    min_node_area = float(np.nanmin(node_area))
    max_node_area = float(np.nanmax(node_area))
    min_area_ratio = min_node_area / max(median_area, 1.0e-300)
    max_area_ratio = max_node_area / max(median_area, 1.0e-300)

    dual_over_edge = dual / np.maximum(edge_len, 1.0e-300)
    dual_over_edge_max = float(np.nanmax(dual_over_edge))
    dual_over_edge_p99 = float(np.nanpercentile(dual_over_edge, 99.0))

    # Symmetric local FV stiffness proxy: worst of the two endpoint volumes.
    a_i = node_area[np.asarray(ops.edge_i, dtype=np.int64)]
    a_j = node_area[np.asarray(ops.edge_j, dtype=np.int64)]
    fv_w_i = dual / np.maximum(a_i * edge_len, 1.0e-300)
    fv_w_j = dual / np.maximum(a_j * edge_len, 1.0e-300)
    fv_weight = np.maximum(fv_w_i, fv_w_j)
    fv_med = float(np.nanmedian(fv_weight))
    fv_max = float(np.nanmax(fv_weight))
    fv_p99 = float(np.nanpercentile(fv_weight, 99.0))
    fv_ratio = fv_max / max(fv_med, 1.0e-300)

    metrics: dict[str, float | int | str] = {
        "n_nodes": int(nodes.shape[0]),
        "n_triangles": int(triangles.shape[0]),
        "n_edges": int(edges.shape[0]),
        "target_spacing_m": float(target_h),
        "edge_min_m": min_edge,
        "edge_median_m": median_edge,
        "edge_max_m": max_edge,
        "edge_max_over_min": float(edge_ratio),
        "triangle_area_min_m2": float(np.nanmin(tri_area)),
        "triangle_area_median_m2": float(np.nanmedian(tri_area)),
        "triangle_area_max_m2": float(np.nanmax(tri_area)),
        "triangle_min_angle_min_deg": float(np.nanmin(tri_min_angle)),
        "triangle_min_angle_p01_deg": float(np.nanpercentile(tri_min_angle, 1.0)),
        "triangle_min_angle_median_deg": float(np.nanmedian(tri_min_angle)),
        "triangles_below_15deg": int(np.count_nonzero(tri_min_angle < 15.0)),
        "triangles_below_20deg": int(np.count_nonzero(tri_min_angle < 20.0)),
        "triangles_below_25deg": int(np.count_nonzero(tri_min_angle < 25.0)),
        "node_area_min_m2": min_node_area,
        "node_area_median_m2": median_area,
        "node_area_max_m2": max_node_area,
        "node_area_min_over_median": float(min_area_ratio),
        "node_area_max_over_median": float(max_area_ratio),
        "dual_over_edge_max": dual_over_edge_max,
        "dual_over_edge_p99": dual_over_edge_p99,
        "fv_weight_max_1_m2": fv_max,
        "fv_weight_p99_1_m2": fv_p99,
        "fv_weight_median_1_m2": fv_med,
        "fv_weight_max_over_median": float(fv_ratio),
    }

    warnings: list[str] = []
    status = "ok"

    def warn_or_fail(condition_warn: bool, condition_fail: bool, message: str) -> None:
        nonlocal status
        if condition_fail:
            status = "fail"
            warnings.append("FAIL: " + message)
        elif condition_warn and status != "fail":
            status = "warn"
            warnings.append("WARN: " + message)

    warn_or_fail(
        float(metrics["triangle_min_angle_min_deg"]) < min_angle_warn_deg,
        float(metrics["triangle_min_angle_min_deg"]) < min_angle_fail_deg,
        f"minimum triangle angle is {metrics['triangle_min_angle_min_deg']:.3g} deg",
    )
    warn_or_fail(
        edge_ratio > max_edge_ratio_warn,
        edge_ratio > max_edge_ratio_fail,
        f"edge max/min ratio is {edge_ratio:.3g}",
    )
    warn_or_fail(
        fv_ratio > max_fv_weight_ratio_warn,
        fv_ratio > max_fv_weight_ratio_fail,
        f"FV stiffness max/median ratio is {fv_ratio:.3g}",
    )
    warn_or_fail(
        min_area_ratio < min_node_area_ratio_warn,
        min_area_ratio < min_node_area_ratio_fail,
        f"node area min/median ratio is {min_area_ratio:.3g}",
    )
    warn_or_fail(
        dual_over_edge_max > max_dual_over_edge_warn,
        dual_over_edge_max > max_dual_over_edge_fail,
        f"dual/edge max ratio is {dual_over_edge_max:.3g}",
    )

    recommendations = _recommendations(status, metrics)
    metrics["status"] = status
    return MeshQualityReport(
        status=status,
        metrics=metrics,
        warnings=warnings,
        recommendations=recommendations,
    )


def save_mesh_quality_report(report: MeshQualityReport, output_path: str | Path) -> Path:
    """Save a mesh-quality report as YAML."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as f:
        yaml.safe_dump(report.to_dict(), f, sort_keys=False, allow_unicode=True)
    return output


def assert_mesh_quality(report: MeshQualityReport) -> None:
    """Raise if the report has fail status."""
    if report.status == "fail":
        joined = "\n".join(report.warnings + report.recommendations)
        raise ValueError("Mesh-quality audit failed:\n" + joined)


def _recommendations(status: str, metrics: dict[str, float | int | str]) -> list[str]:
    rec: list[str] = []
    if status == "ok":
        rec.append("Mesh quality is acceptable for OE7 stationary relaxation.")
        return rec
    rec.append("For OE7 SS, prefer jitter_fraction <= 0.05 while debugging boundary/current conservation.")
    rec.append("Use boundary_guard_layers >= 2 to keep top/bottom and terminal stencils regular.")
    rec.append("If FV stiffness remains high, reduce target_spacing_m or disable interior jitter for this diagnostic run.")
    if float(metrics.get("triangle_min_angle_min_deg", 90.0)) < 20.0:
        rec.append("Poor minimum triangle angle: regenerate with lower jitter or stronger structured protection.")
    if float(metrics.get("fv_weight_max_over_median", 1.0)) > 40.0:
        rec.append("Large FV stiffness outlier: inspect the node/edge with max s/(a*ell); it can trigger Q and j_s spikes.")
    return rec


def _triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]
    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def _triangle_min_angles_deg(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    p = nodes[triangles]
    a = np.linalg.norm(p[:, 1] - p[:, 2], axis=1)
    b = np.linalg.norm(p[:, 0] - p[:, 2], axis=1)
    c = np.linalg.norm(p[:, 0] - p[:, 1], axis=1)
    eps = 1.0e-300
    A = np.degrees(np.arccos(np.clip((b*b + c*c - a*a) / np.maximum(2*b*c, eps), -1.0, 1.0)))
    B = np.degrees(np.arccos(np.clip((a*a + c*c - b*b) / np.maximum(2*a*c, eps), -1.0, 1.0)))
    C = 180.0 - A - B
    return np.minimum.reduce([A, B, C])
