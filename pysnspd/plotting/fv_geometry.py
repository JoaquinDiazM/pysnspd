"""Finite-volume geometry diagnostics for pySNSPD OE7.

This module is deliberately diagnostic-only: it does not change the operators
used by the stationary solver.  The goal is to compare the *current* control
volumes used by ``FVOperators`` against the circumcentric/Voronoi dual geometry
that is implied by the Delaunay finite-volume stencil.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
import yaml


_NM = 1.0e9
_EPS = 1.0e-300


def compute_fv_geometry_audit(mesh, edge_data=None, ops=None) -> dict[str, Any]:
    """Return arrays and summaries for the FV geometry audit.

    Parameters
    ----------
    mesh:
        pySNSPD mesh object with ``nodes`` and ``triangles``.
    edge_data:
        Optional edge-data object.  If ``ops`` is missing, this is used to build
        the current FV operators with the repository implementation.
    ops:
        Optional ``FVOperators`` object.  When provided, its node areas and
        edge dual lengths are treated as the *current solver geometry*.

    Notes
    -----
    The reference area is an independent circumcentric/Voronoi control-volume
    estimate built directly from triangle circumcenters and edge midpoints.  In
    a perfectly pyTDGL-like geometry, the ratio

        reference_circumcentric_area / current_operator_area

    should be close to one away from boundary/clipping artifacts.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] < 2:
        raise ValueError("mesh.nodes must have shape (n_nodes, >=2).")
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError("mesh.triangles must have shape (n_triangles, 3).")

    if ops is None:
        if edge_data is None:
            raise ValueError("edge_data is required when ops is not provided.")
        from pysnspd.gtdgl.operators import build_fv_operators

        ops = build_fv_operators(mesh, edge_data)

    edges = np.asarray(ops.edges, dtype=np.int64)
    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    edge_length = np.asarray(ops.edge_length_m, dtype=float)
    current_area = np.asarray(ops.node_area_m2, dtype=float)
    current_dual = np.asarray(ops.dual_face_length_m, dtype=float)

    bary_area = barycentric_node_areas(nodes, triangles)
    circum = triangle_circumcenters(nodes, triangles)
    edge_triangles = edge_triangles_from_mesh(triangles, edges)
    reference_dual = circumcenter_dual_lengths(
        nodes=nodes,
        edges=edges,
        edge_length=edge_length,
        circumcenters=circum,
        edge_triangles=edge_triangles,
    )
    reference_area_signed, reference_area_abs = circumcentric_node_areas(nodes, triangles, circum)

    area_ratio_signed = reference_area_signed / np.maximum(current_area, _EPS)
    area_ratio_abs = reference_area_abs / np.maximum(current_area, _EPS)
    dual_ratio = reference_dual / np.maximum(current_dual, _EPS)
    edge_weight_current = current_dual / np.maximum(edge_length, _EPS)
    edge_weight_reference = reference_dual / np.maximum(edge_length, _EPS)

    node_diag_current = incident_sum(edge_i, edge_j, edge_weight_current, n_nodes=nodes.shape[0])
    node_diag_reference = incident_sum(edge_i, edge_j, edge_weight_reference, n_nodes=nodes.shape[0])
    laplace_scale_current = node_diag_current / np.maximum(current_area, _EPS)
    laplace_scale_reference_on_ref_area = node_diag_reference / np.maximum(np.abs(reference_area_signed), _EPS)

    summary = {
        "n_nodes": int(nodes.shape[0]),
        "n_triangles": int(triangles.shape[0]),
        "n_edges": int(edges.shape[0]),
        "control_volume_policy": {
            "current_operator_area": "ops.node_area_m2",
            "reference_area": "circumcentric/Voronoi cell estimate from triangle circumcenters",
            "diagnostic_ratio": "reference_area / current_operator_area",
        },
        "current_area_m2": percentiles(current_area),
        "barycentric_area_m2": percentiles(bary_area),
        "reference_circumcentric_area_signed_m2": percentiles(reference_area_signed),
        "reference_circumcentric_area_abs_m2": percentiles(reference_area_abs),
        "reference_to_current_area_ratio_signed": percentiles(area_ratio_signed),
        "reference_abs_to_current_area_ratio": percentiles(area_ratio_abs),
        "negative_reference_area_nodes": int(np.count_nonzero(reference_area_signed <= 0.0)),
        "negative_reference_area_fraction": float(np.mean(reference_area_signed <= 0.0)),
        "edge_weight_current_s_over_l": percentiles(edge_weight_current),
        "edge_weight_reference_s_over_l": percentiles(edge_weight_reference),
        "reference_to_current_dual_length_ratio": percentiles(dual_ratio),
        "laplace_scale_current_1_m2": percentiles(laplace_scale_current),
        "laplace_scale_reference_on_reference_area_1_m2": percentiles(laplace_scale_reference_on_ref_area),
        "laplace_scale_current_over_median": percentiles(
            laplace_scale_current / np.maximum(np.nanmedian(laplace_scale_current), _EPS)
        ),
    }

    return {
        "nodes": nodes,
        "triangles": triangles,
        "edges": edges,
        "edge_i": edge_i,
        "edge_j": edge_j,
        "edge_length_m": edge_length,
        "current_area_m2": current_area,
        "barycentric_area_m2": bary_area,
        "reference_area_signed_m2": reference_area_signed,
        "reference_area_abs_m2": reference_area_abs,
        "area_ratio_signed": area_ratio_signed,
        "area_ratio_abs": area_ratio_abs,
        "current_dual_length_m": current_dual,
        "reference_dual_length_m": reference_dual,
        "dual_ratio": dual_ratio,
        "edge_weight_current": edge_weight_current,
        "edge_weight_reference": edge_weight_reference,
        "laplace_scale_current_1_m2": laplace_scale_current,
        "laplace_scale_reference_on_ref_area_1_m2": laplace_scale_reference_on_ref_area,
        "summary": summary,
    }


def write_fv_geometry_audit_yaml(audit: dict[str, Any], output_path: str | Path) -> Path:
    """Write the scalar FV geometry audit summary as YAML."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(audit["summary"], f, sort_keys=False)
    return output_path


def plot_fv_geometry_audit(
    mesh,
    edge_data=None,
    ops=None,
    output_path: str | Path = "fv_geometry_audit.png",
    *,
    dpi: int = 480,
) -> Path:
    """Plot FV geometry diagnostics for the current solver operators."""
    audit = compute_fv_geometry_audit(mesh, edge_data=edge_data, ops=ops)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    nodes = audit["nodes"]
    triangles = audit["triangles"]
    edges = audit["edges"]
    x_nm = _NM * nodes[:, 0]
    y_nm = _NM * nodes[:, 1]
    triang = mtri.Triangulation(x_nm, y_nm, triangles)

    area_ratio_signed = np.asarray(audit["area_ratio_signed"], dtype=float)
    edge_weight_current = np.asarray(audit["edge_weight_current"], dtype=float)
    laplace_scale = np.asarray(audit["laplace_scale_current_1_m2"], dtype=float)
    dual_ratio = np.asarray(audit["dual_ratio"], dtype=float)

    log_area_ratio = signed_log10(area_ratio_signed)
    log_edge_weight = np.log10(np.maximum(edge_weight_current, _EPS))
    log_laplace_scale = np.log10(
        np.maximum(laplace_scale / np.maximum(np.nanmedian(laplace_scale), _EPS), _EPS)
    )
    log_dual_ratio = np.log10(np.maximum(dual_ratio, _EPS))

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.8), constrained_layout=True)
    fig.suptitle("OE7 FV geometry audit: current operators vs circumcentric dual reference", y=1.015)

    ax = axes[0, 0]
    vmax = robust_symmetric_limit(log_area_ratio, default=1.0)
    im = ax.tripcolor(triang, log_area_ratio, shading="gouraud", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    mark_negative_reference_nodes(ax, x_nm, y_nm, audit["reference_area_signed_m2"])
    ax.set_title(r"node control volume: signed $\log_{10}(A_{circ}/A_{op})$")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(im, ax=ax, shrink=0.86)
    cb.set_label(r"signed $\log_{10}$ ratio")

    ax = axes[0, 1]
    im = edge_line_collection(
        ax,
        x_nm,
        y_nm,
        edges,
        log_edge_weight,
        cmap="viridis",
        symmetric=False,
    )
    ax.set_title(r"current edge weight $\log_{10}(s_{ij}/e_{ij})$")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(im, ax=ax, shrink=0.86)
    cb.set_label(r"$\log_{10}(s_{ij}/e_{ij})$")

    ax = axes[1, 0]
    vmax = robust_symmetric_limit(log_laplace_scale, default=1.0)
    im = ax.tripcolor(triang, log_laplace_scale, shading="gouraud", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_title(r"node stiffness scale: $\log_{10}[(\sum s/e)/A_{op}]$ relative to median")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    cb = fig.colorbar(im, ax=ax, shrink=0.86)
    cb.set_label(r"relative $\log_{10}$ scale")

    ax = axes[1, 1]
    finite_area = np.isfinite(log_area_ratio)
    finite_weight = np.isfinite(log_edge_weight)
    finite_dual = np.isfinite(log_dual_ratio)
    ax.hist(log_area_ratio[finite_area], bins=60, histtype="step", linewidth=1.8, label=r"$A_{circ}/A_{op}$")
    ax.hist(log_edge_weight[finite_weight], bins=60, histtype="step", linewidth=1.8, label=r"$s_{ij}/e_{ij}$")
    ax.hist(log_dual_ratio[finite_dual], bins=60, histtype="step", linewidth=1.8, label=r"$s_{circ}/s_{op}$")
    ax.axvline(0.0, linewidth=1.0, alpha=0.65)
    ax.set_title("geometry-ratio histograms")
    ax.set_xlabel(r"$\log_{10}$ value")
    ax.set_ylabel("count")
    ax.legend(frameon=False)
    ax.grid(False)

    for ax in axes.flat:
        ax.grid(False)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output_path


def barycentric_node_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return barycentric lumped nodal areas."""
    area = np.zeros(nodes.shape[0], dtype=float)
    tri_area = np.abs(triangle_signed_areas(nodes, triangles))
    for local in range(3):
        np.add.at(area, triangles[:, local], tri_area / 3.0)
    return area


def circumcentric_node_areas(
    nodes: np.ndarray,
    triangles: np.ndarray,
    circumcenters: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return signed and absolute circumcentric/Voronoi nodal areas.

    Each triangle contributes the polygon
    ``[vertex, midpoint(vertex,next), circumcenter, midpoint(vertex,prev)]``
    to the corresponding vertex.  The signed contribution is oriented by the
    parent triangle orientation, so obtuse/ill-conditioned cells can expose
    negative or extreme values instead of being silently hidden.
    """
    pts = np.asarray(nodes, dtype=float)[:, :2]
    tris = np.asarray(triangles, dtype=np.int64)
    if circumcenters is None:
        circumcenters = triangle_circumcenters(nodes, triangles)

    signed = np.zeros(pts.shape[0], dtype=float)
    absolute = np.zeros(pts.shape[0], dtype=float)

    tri_signed = triangle_signed_areas(nodes, triangles)
    for t, tri in enumerate(tris):
        sign = 1.0 if tri_signed[t] >= 0.0 else -1.0
        cc = circumcenters[t]
        for local in range(3):
            i = int(tri[local])
            j_next = int(tri[(local + 1) % 3])
            j_prev = int(tri[(local - 1) % 3])
            p_i = pts[i]
            m_next = 0.5 * (pts[i] + pts[j_next])
            m_prev = 0.5 * (pts[i] + pts[j_prev])
            poly = np.vstack([p_i, m_next, cc, m_prev])
            a = polygon_signed_area(poly) * sign
            signed[i] += a
            absolute[i] += abs(a)
    return signed, absolute


def triangle_signed_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return signed triangle areas."""
    p = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(triangles, dtype=np.int64)
    p0 = p[tri[:, 0]]
    p1 = p[tri[:, 1]]
    p2 = p[tri[:, 2]]
    return 0.5 * ((p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0]))


def triangle_circumcenters(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return triangle circumcenters."""
    pts = np.asarray(nodes, dtype=float)[:, :2]
    out = np.zeros((triangles.shape[0], 2), dtype=float)
    for k, tri in enumerate(np.asarray(triangles, dtype=np.int64)):
        out[k] = circumcenter(pts[tri])
    return out


def circumcenter(p: np.ndarray) -> np.ndarray:
    """Return the circumcenter of a 2D triangle; centroid fallback if singular."""
    a, b, c = p[0], p[1], p[2]
    d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]))
    if abs(d) < _EPS:
        return np.mean(p, axis=0)
    a2 = float(np.dot(a, a))
    b2 = float(np.dot(b, b))
    c2 = float(np.dot(c, c))
    ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d
    uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d
    return np.array([ux, uy], dtype=float)


def edge_triangles_from_mesh(triangles: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Return up to two adjacent triangle indices for each edge."""
    edge_lookup = {tuple(sorted(map(int, e))): k for k, e in enumerate(np.asarray(edges, dtype=np.int64))}
    adj: list[list[int]] = [[] for _ in range(len(edge_lookup))]
    for t, tri in enumerate(np.asarray(triangles, dtype=np.int64)):
        for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            k = edge_lookup.get(tuple(sorted((int(u), int(v)))))
            if k is not None:
                adj[k].append(int(t))
    out = -np.ones((len(adj), 2), dtype=np.int64)
    for k, ts in enumerate(adj):
        if ts:
            out[k, : min(2, len(ts))] = ts[:2]
    return out


def circumcenter_dual_lengths(
    *,
    nodes: np.ndarray,
    edges: np.ndarray,
    edge_length: np.ndarray,
    circumcenters: np.ndarray,
    edge_triangles: np.ndarray,
) -> np.ndarray:
    """Return circumcentric dual-face lengths for each primal edge."""
    pts = np.asarray(nodes, dtype=float)[:, :2]
    dual = np.zeros(edges.shape[0], dtype=float)
    for k, (u, v) in enumerate(np.asarray(edges, dtype=np.int64)):
        ts = np.asarray(edge_triangles[k], dtype=np.int64)
        ts = ts[ts >= 0]
        mid = 0.5 * (pts[int(u)] + pts[int(v)])
        if ts.size == 1:
            sij = float(np.linalg.norm(circumcenters[int(ts[0])] - mid))
        elif ts.size >= 2:
            sij = float(np.linalg.norm(circumcenters[int(ts[0])] - circumcenters[int(ts[1])]))
        else:
            sij = 0.5 * float(edge_length[k])
        if not np.isfinite(sij) or sij <= 0.0:
            sij = 0.5 * float(edge_length[k])
        dual[k] = sij
    return dual


def incident_sum(edge_i: np.ndarray, edge_j: np.ndarray, values: np.ndarray, *, n_nodes: int) -> np.ndarray:
    """Sum edge values on incident nodes."""
    out = np.zeros(int(n_nodes), dtype=float)
    np.add.at(out, edge_i, values)
    np.add.at(out, edge_j, values)
    return out


def polygon_signed_area(poly: np.ndarray) -> float:
    """Signed area of a 2D polygon."""
    x = poly[:, 0]
    y = poly[:, 1]
    return float(0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1)))


def signed_log10(values: np.ndarray) -> np.ndarray:
    """Signed log10(abs(values)) with sign retained."""
    v = np.asarray(values, dtype=float)
    out = np.zeros_like(v, dtype=float)
    mask = np.isfinite(v) & (v != 0.0)
    out[mask] = np.sign(v[mask]) * np.log10(np.maximum(np.abs(v[mask]), _EPS))
    out[~np.isfinite(v)] = np.nan
    return out


def robust_symmetric_limit(values: np.ndarray, *, default: float = 1.0) -> float:
    """Robust symmetric color limit."""
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v)
    if not np.any(finite):
        return float(default)
    lim = float(np.nanpercentile(np.abs(v[finite]), 99.0))
    return max(lim, float(default) * 1.0e-6)


def edge_line_collection(
    ax,
    x_nm: np.ndarray,
    y_nm: np.ndarray,
    edges: np.ndarray,
    values: np.ndarray,
    *,
    cmap: str,
    symmetric: bool,
):
    """Draw edge-colored line collection."""
    segments = np.stack(
        [
            np.column_stack([x_nm[edges[:, 0]], y_nm[edges[:, 0]]]),
            np.column_stack([x_nm[edges[:, 1]], y_nm[edges[:, 1]]]),
        ],
        axis=1,
    )
    finite = np.isfinite(values)
    if symmetric:
        vmax = robust_symmetric_limit(values, default=1.0)
        vmin = -vmax
    else:
        vmin = float(np.nanpercentile(values[finite], 1.0)) if np.any(finite) else -1.0
        vmax = float(np.nanpercentile(values[finite], 99.0)) if np.any(finite) else 1.0
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin, vmax = -1.0, 1.0
    lc = LineCollection(segments, cmap=cmap, linewidths=0.55, alpha=0.95)
    lc.set_array(values)
    lc.set_clim(vmin, vmax)
    ax.add_collection(lc)
    ax.autoscale()
    return lc


def mark_negative_reference_nodes(ax, x_nm: np.ndarray, y_nm: np.ndarray, reference_area_signed: np.ndarray) -> None:
    """Overlay markers on nodes with non-positive signed circumcentric area."""
    neg = np.asarray(reference_area_signed, dtype=float) <= 0.0
    if np.any(neg):
        ax.scatter(x_nm[neg], y_nm[neg], s=4.0, marker="x", linewidths=0.45, alpha=0.9)


def percentiles(values: np.ndarray) -> dict[str, float]:
    """Return stable scalar percentiles as plain Python floats."""
    v = np.asarray(values, dtype=float)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        return {key: float("nan") for key in ("min", "p01", "p05", "p50", "p95", "p99", "max")}
    qs = np.percentile(finite, [0, 1, 5, 50, 95, 99, 100])
    return {
        "min": float(qs[0]),
        "p01": float(qs[1]),
        "p05": float(qs[2]),
        "p50": float(qs[3]),
        "p95": float(qs[4]),
        "p99": float(qs[5]),
        "max": float(qs[6]),
    }
