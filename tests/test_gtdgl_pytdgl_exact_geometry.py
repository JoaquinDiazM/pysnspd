"""Tests for the exact pyTDGL geometry/meshing compatibility path."""

from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.geometry import box, ensure_unique
from pysnspd.gtdgl.tdgl_compat import Device, Layer, Polygon
from pysnspd.mesh.pytdgl_like import rectangular_boundary_points


def test_box_matches_pytdgl_rectangle_point_count_distribution():
    coords = box(
        width=2.4e-7,
        height=1.2e-7,
        points=101,
        center=(1.2e-7, 0.0),
    )

    perimeter = 2.0 * (2.4e-7 + 1.2e-7)
    x_points = round(101 * 2.4e-7 / perimeter)
    y_points = round(101 * 1.2e-7 / perimeter)

    assert coords.shape == (2 * (x_points + y_points), 2)
    assert np.isclose(coords[:, 0].min(), 0.0)
    assert np.isclose(coords[:, 0].max(), 2.4e-7)
    assert np.isclose(coords[:, 1].min(), -0.6e-7)
    assert np.isclose(coords[:, 1].max(), 0.6e-7)


def test_rectangular_boundary_points_uses_pytdgl_box_coordinate_convention():
    coords = rectangular_boundary_points(
        length_m=2.4e-7,
        width_m=1.2e-7,
        points=101,
    )

    expected = box(
        width=2.4e-7,
        height=1.2e-7,
        points=101,
        center=(1.2e-7, 0.0),
    )

    assert np.allclose(coords, expected)

    # pyTDGL's box() convention can repeat the first coordinate at the end for
    # rectangular boundaries because the side linspace calls include endpoints.
    # That is fine: generate_mesh() calls ensure_unique() before building facets.
    assert np.allclose(coords[0], coords[-1])
    assert ensure_unique(coords).shape[0] < coords.shape[0]


def test_rectangular_boundary_unique_points_still_define_box_extent():
    coords = rectangular_boundary_points(
        length_m=2.4e-7,
        width_m=1.2e-7,
        points=101,
    )

    unique = ensure_unique(coords)

    assert unique.shape[0] < coords.shape[0]
    assert np.isclose(unique[:, 0].min(), 0.0)
    assert np.isclose(unique[:, 0].max(), 2.4e-7)
    assert np.isclose(unique[:, 1].min(), -0.6e-7)
    assert np.isclose(unique[:, 1].max(), 0.6e-7)


def test_minimal_pytdgl_device_copy_and_contains_points():
    film = Polygon(
        "film",
        points=rectangular_boundary_points(
            length_m=2.4e-7,
            width_m=1.2e-7,
            points=101,
        ),
    )

    device = Device(
        "test",
        layer=Layer(
            london_lambda=1.0,
            coherence_length=1.0,
            thickness=1.0,
        ),
        film=film,
        length_units="m",
    )

    points = np.array(
        [
            [1.2e-7, 0.0],
            [4.0e-7, 0.0],
        ]
    )

    mask = device.contains_points(points)

    assert mask.dtype == np.bool_
    assert mask.tolist() == [True, False]

    clone = device.copy(with_mesh=False)

    assert clone is not device
    assert clone.mesh is None
    assert clone.name == device.name
    assert np.allclose(clone.film.points, device.film.points)