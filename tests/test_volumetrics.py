"""Unit tests verifying classical Surveying & Geoinformatics volumetric theorems."""

import numpy as np
from src.volumetric_methods import (
    analytical_frustum_volume,
    prismoidal_correction,
    end_area_volume,
    prismoidal_volume,
    borrow_pit_grid_volume,
    delaunay_tin_volume,
)


def test_prismoidal_correction_identity():
    """
    Classical Surveying Theorem:
    For a two-section rectangular frustum of top (L1 x W1), bottom (L2 x W2), and height D,
    Average End-Area systematically overestimates true prismoidal volume by exactly Cp:
        V_end_area - Cp == V_prismoidal_exact
    """
    L1, W1 = 40.0, 30.0
    L2, W2 = 24.0, 16.0
    D = 12.0

    v_end = 0.5 * D * (L1 * W1 + L2 * W2)
    cp = prismoidal_correction(L1, W1, L2, W2, D)
    v_exact = analytical_frustum_volume(L1, W1, L2, W2, D)

    assert cp > 0.0
    assert np.isclose(v_end - cp, v_exact, rtol=1e-12)


def test_uniform_rectangular_box_exact_across_all_surface_methods():
    """On a flat-bottomed rectangular box, continuous-grid methods agree exactly with L*W*D."""
    xs = np.linspace(0.0, 20.0, 11)
    ys = np.linspace(0.0, 15.0, 11)
    XX, YY = np.meshgrid(xs, ys)
    depth = np.full_like(XX, 8.5)
    expected = 20.0 * 15.0 * 8.5  # 2550.0 m^3

    assert np.isclose(end_area_volume(xs, ys, depth), expected, rtol=1e-12)
    assert np.isclose(prismoidal_volume(xs, ys, depth), expected, rtol=1e-12)
    assert np.isclose(borrow_pit_grid_volume(xs, ys, depth), expected, rtol=1e-12)
    assert np.isclose(delaunay_tin_volume(XX.ravel(), YY.ravel(), depth.ravel()), expected, rtol=1e-12)


def test_simpson_prismoidal_exact_on_quadratic_paraboloid():
    """Prismoidal (Simpson 1/3) is machine-exact for any biquadratic excavation bowl."""
    xs = np.linspace(0.0, 10.0, 11)
    ys = np.linspace(0.0, 10.0, 11)
    XX, YY = np.meshgrid(xs, ys)
    # D(x, y) = x*(10 - x) + y*(10 - y); integral over [0,10]x[0,10]:
    # int_0^10 x(10-x) dx = [5x^2 - x^3/3]_0^10 = 500 - 1000/3 = 500/3
    # Double integral = 10*(500/3) + 10*(500/3) = 10000/3
    depth = XX * (10.0 - XX) + YY * (10.0 - YY)
    expected = 10000.0 / 3.0

    v_prism = prismoidal_volume(xs, ys, depth)
    v_trap = end_area_volume(xs, ys, depth)

    assert np.isclose(v_prism, expected, rtol=1e-12)
    assert abs(v_prism - expected) < abs(v_trap - expected)
