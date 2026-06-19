from pathlib import Path

import numpy as np

from pysnspd.config import validate_config
from pysnspd.mesh.delaunay import (
    generate_rectangular_delaunay_mesh,
    mesh_summary,
    save_mesh_npz,
    load_mesh_npz,
)
from pysnspd.mesh.edges import (
    assert_edge_data_consistent,
    build_edge_data,
    edge_summary,
    save_edges_npz,
    load_edges_npz,
)
from pysnspd.plotting.figures import (
    plot_boundary_tags,
    plot_mesh_geometry,
)


def _minimal_config(tmp_path: Path) -> dict:
    return {
        "project": {
            "name": "test_project",
            "big_data_root": str(tmp_path / "big_data"),
            "default_run_name": "mesh_test",
        },
        "parallel": {
            "enabled": True,
            "workers": 2,
            "backend": "process",
        },
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "sigma_n_S_m": "4.2e5",
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
        },
        "calibration": {
            "Ic_target_A": 38.8e-6,
            "n_gamma_sweep": 40,
            "gamma_max_fraction": 0.80,
            "D_warn_min_m2_s": 5.0e-5,
            "D_warn_max_m2_s": 5.0e-4,
        },
        "bias": {
            "T_bias_K": 0.9,
            "I_bias_A": 35.0e-6,
        },
        "mesh": {
            "type": "delaunay",
            "length_m": 360.0e-9,
            "target_spacing_m": 20.0e-9,
            "seed": 12345,
        },
        "catalogs": {
            "dos": {
                "n_delta": 5,
                "n_q": 6,
                "n_energy": 50,
                "n_matsubara": 20,
            },
            "phase_space": {
                "n_Te": 5,
                "n_Tph": 5,
                "n_delta": 5,
                "n_q": 6,
                "n_omega": 50,
            },
        },
        "ss_run": {
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "convergence_tol": 1.0e-7,
        },
        "photon_run": {
            "photon_wavelength_m": 1064.0e-9,
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "bubble_radius_m": 10.0e-9,
        },
        "circuit": {
            "R_load_ohm": 50.0,
            "L_bias_H": 1.0e-6,
            "C_rf_F": 1.0e-12,
        },
    }


def test_rectangular_delaunay_mesh_is_reproducible(tmp_path):
    cfg = validate_config(
        _minimal_config(tmp_path),
        require_big_data_root_exists=False,
    )

    mesh_a = generate_rectangular_delaunay_mesh(cfg)
    mesh_b = generate_rectangular_delaunay_mesh(cfg)

    assert np.allclose(mesh_a.nodes, mesh_b.nodes)
    assert np.array_equal(mesh_a.triangles, mesh_b.triangles)

    summary = mesh_summary(mesh_a)
    assert summary["n_nodes"] > 0
    assert summary["n_triangles"] > 0
    assert summary["triangle_area_min_m2"] > 0.0


def test_edge_extraction_and_boundary_tags(tmp_path):
    cfg = validate_config(
        _minimal_config(tmp_path),
        require_big_data_root_exists=False,
    )

    mesh = generate_rectangular_delaunay_mesh(cfg)
    edges = build_edge_data(
        mesh.nodes,
        mesh.triangles,
        length_m=mesh.length_m,
        width_m=mesh.width_m,
    )

    assert_edge_data_consistent(edges)

    summary = edge_summary(edges)
    assert summary["n_edges"] > 0
    assert summary["n_boundary_edges"] > 0
    assert summary["tag_counts"].get("left", 0) > 0
    assert summary["tag_counts"].get("right", 0) > 0
    assert summary["tag_counts"].get("top", 0) > 0
    assert summary["tag_counts"].get("bottom", 0) > 0
    assert summary["tag_counts"].get("boundary_unknown", 0) == 0


def test_mesh_save_load_and_plots(tmp_path):
    cfg = validate_config(
        _minimal_config(tmp_path),
        require_big_data_root_exists=False,
    )

    mesh = generate_rectangular_delaunay_mesh(cfg)
    edges = build_edge_data(
        mesh.nodes,
        mesh.triangles,
        length_m=mesh.length_m,
        width_m=mesh.width_m,
    )

    mesh_path = save_mesh_npz(mesh, tmp_path / "mesh.npz")
    edges_path = save_edges_npz(edges, tmp_path / "edges.npz")

    mesh_loaded = load_mesh_npz(mesh_path)
    edges_loaded = load_edges_npz(edges_path)

    assert np.allclose(mesh.nodes, mesh_loaded.nodes)
    assert np.array_equal(mesh.triangles, mesh_loaded.triangles)
    assert np.array_equal(edges.edges, edges_loaded.edges)
    assert np.array_equal(edges.tags, edges_loaded.tags)

    mesh_plot = plot_mesh_geometry(
        mesh_loaded,
        edges_loaded,
        tmp_path / "mesh_nodes_edges.png",
    )
    tags_plot = plot_boundary_tags(
        mesh_loaded,
        edges_loaded,
        tmp_path / "mesh_boundary_tags.png",
    )

    assert mesh_plot.exists()
    assert tags_plot.exists()


def test_mesh_uses_all_nodes(tmp_path):
    cfg = validate_config(
        _minimal_config(tmp_path),
        require_big_data_root_exists=False,
    )

    mesh = generate_rectangular_delaunay_mesh(cfg)
    used_nodes = np.unique(mesh.triangles.reshape(-1))

    assert used_nodes.size == mesh.n_nodes

    summary = mesh_summary(mesh)
    assert summary["n_unused_nodes"] == 0
    assert summary["area_relative_error"] < 1.0e-12
