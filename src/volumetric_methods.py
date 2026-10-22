"""
Five Classic Surveying & Geoinformatics Earthwork Volumetric Formulations

Implements and compares five standard volume estimation methods studied in
Surveying & Geoinformatics engineering practice:

  1. End-Area (Trapezoidal Cross-Section) Method
  2. Prismoidal (Composite Simpson's 1/3 Cross-Section) Method
  3. Borrow-Pit Rectangular Grid (1-2-3-4 Corner Weight) Method
  4. Horizontal Contour Planimetric Slicing Method
  5. Triangulated Irregular Network (Delaunay TIN Truncated Prism) Method

Also provides exact analytical validation formulas for rectangular frustum pits
and the classical Surveying Prismoidal Correction (Cp).
"""

from typing import Dict
import numpy as np
from scipy.spatial import Delaunay


def analytical_frustum_volume(L1: float, W1: float, L2: float, W2: float, D: float) -> float:
    """
    Exact volume of a rectangular frustum excavation of depth D with top dimensions
    (L1 x W1) and base dimensions (L2 x W2):
        V = (D / 6) * [A1 + 4 * A_mid + A2]
    """
    A1 = L1 * W1
    A2 = L2 * W2
    A_mid = ((L1 + L2) / 2.0) * ((W1 + W2) / 2.0)
    return (D / 6.0) * (A1 + 4.0 * A_mid + A2)


def prismoidal_correction(L1: float, W1: float, L2: float, W2: float, D: float) -> float:
    """
    Classical Surveying Prismoidal Correction subtracted from Two-Section End-Area Volume:
        Cp = (D / 6) * (L1 - L2) * (W1 - W2)
        V_prismoidal = V_end_area - Cp
    """
    return (D / 6.0) * (L1 - L2) * (W1 - W2)


def end_area_volume(xs: np.ndarray, ys: np.ndarray, depth_grid: np.ndarray) -> float:
    """
    Method 1: Average End-Area (Trapezoidal Cross-Section) Rule.
    Computes vertical cross-sectional cut areas A_j along X via trapezoidal integration,
    then averages adjacent cross-sections separated by dy:
        V = sum_j (dy / 2) * (A_j + A_{j+1})
    """
    dx = float(xs[1] - xs[0])
    dy = float(ys[1] - ys[0])
    areas = np.trapezoid(depth_grid, dx=dx, axis=1)
    return float(np.sum(0.5 * dy * (areas[:-1] + areas[1:])))


def prismoidal_volume(xs: np.ndarray, ys: np.ndarray, depth_grid: np.ndarray) -> float:
    """
    Method 2: Prismoidal / Composite Simpson's 1/3 Cross-Section Method.
    Computes cross-sectional areas using Simpson's 1/3 rule along X (falling back cleanly
    on last panel if station count is even) and integrates cross-sections via Prismoidal
    formula V = (dy / 3) * (A_0 + 4*sum(A_odd) + 2*sum(A_even) + A_M).
    """
    def simpson_1d(vals: np.ndarray, h: float) -> float:
        n = len(vals) - 1  # number of intervals
        if n < 2:
            return float(np.trapezoid(vals, dx=h))
        if n % 2 == 0:
            return float((h / 3.0) * (vals[0] + vals[-1] + 4.0 * np.sum(vals[1:-1:2]) + 2.0 * np.sum(vals[2:-2:2])))
        # Composite Simpson on first (n-1) even intervals + Simpson 3/8 or trapezoidal tail
        even_part = (h / 3.0) * (vals[0] + vals[-2] + 4.0 * np.sum(vals[1:-2:2]) + 2.0 * np.sum(vals[2:-3:2]))
        tail = 0.5 * h * (vals[-2] + vals[-1])
        return float(even_part + tail)

    dx = float(xs[1] - xs[0])
    dy = float(ys[1] - ys[0])
    areas = np.array([simpson_1d(row, dx) for row in depth_grid])
    return simpson_1d(areas, dy)


def borrow_pit_grid_volume(xs: np.ndarray, ys: np.ndarray, depth_grid: np.ndarray) -> float:
    """
    Method 3: Spot-Height Rectangular Grid ('Borrow-Pit' 1-2-3-4 Weight) Method.
    Counts how many active rectangular cells each node participates in (w_ij in {1, 2, 3, 4}):
        V = (dx * dy / 4) * sum(h_1 + 2*h_2 + 3*h_3 + 4*h_4)
    """
    dx = float(xs[1] - xs[0])
    dy = float(ys[1] - ys[0])
    cell_area = dx * dy

    # Count membership of every vertex across all grid cells
    ny, nx = depth_grid.shape
    weights = np.zeros_like(depth_grid, dtype=float)
    weights[:-1, :-1] += 1.0
    weights[:-1, 1:] += 1.0
    weights[1:, :-1] += 1.0
    weights[1:, 1:] += 1.0

    return float((cell_area / 4.0) * np.sum(weights * depth_grid))


def contour_slicing_volume(
    xs: np.ndarray, ys: np.ndarray, depth_grid: np.ndarray, contour_interval: float = 1.0
) -> float:
    """
    Method 4: Horizontal Contour Planimetric Integration.
    Interpolates the survey grid bilinearly onto a high-resolution planimetric surface
    (simulating digital planimetering of contour lines in GIS/CAD) to compute enclosed
    horizontal slice areas A(z) at interval `contour_interval`, then integrates via
    the frustum/cone surveying contour formula V = sum (dz/3)*(A_k + A_{k+1} + sqrt(A_k*A_{k+1})).
    """
    from scipy.interpolate import RegularGridInterpolator

    dmax = float(np.max(depth_grid))
    if dmax <= 0.0:
        return 0.0

    interp = RegularGridInterpolator((ys, xs), depth_grid, method="linear")
    nx_fine, ny_fine = 240, 200
    xf = np.linspace(float(xs[0]), float(xs[-1]), nx_fine, endpoint=False) + (xs[-1] - xs[0]) / (2.0 * nx_fine)
    yf = np.linspace(float(ys[0]), float(ys[-1]), ny_fine, endpoint=False) + (ys[-1] - ys[0]) / (2.0 * ny_fine)
    YYf, XXf = np.meshgrid(yf, xf, indexing="ij")
    fine_d = interp((YYf, XXf))
    pixel_area = float((xs[-1] - xs[0]) * (ys[-1] - ys[0])) / (nx_fine * ny_fine)

    levels = np.append(np.arange(0.0, dmax, contour_interval), dmax)
    slice_areas = np.array([
        float(np.sum(fine_d > 1e-6 if z == 0.0 else fine_d >= z)) * pixel_area
        for z in levels
    ])

    vol = 0.0
    for k in range(len(levels) - 1):
        dz = float(levels[k + 1] - levels[k])
        a1, a2 = slice_areas[k], slice_areas[k + 1]
        vol += (dz / 3.0) * (a1 + a2 + np.sqrt(max(0.0, a1 * a2)))
    return float(vol)


def delaunay_tin_volume(x_flat: np.ndarray, y_flat: np.ndarray, depth_flat: np.ndarray) -> float:
    """
    Method 5: Triangulated Irregular Network (Delaunay TIN Truncated Triangular Prism) Method.
    Constructs a planar 2D Delaunay triangulation over survey (X, Y) coordinates and sums
    exact truncated triangular prism volumes:
        V_TIN = sum_{tri} Area_2D(tri) * (d_A + d_B + d_C) / 3
    """
    pts_2d = np.column_stack([x_flat, y_flat])
    tri = Delaunay(pts_2d)
    simplices = tri.simplices  # (N_tri, 3)

    pA = pts_2d[simplices[:, 0]]
    pB = pts_2d[simplices[:, 1]]
    pC = pts_2d[simplices[:, 2]]

    # Exact 2D triangle area via cross-product determinant
    area_2d = 0.5 * np.abs(
        (pB[:, 0] - pA[:, 0]) * (pC[:, 1] - pA[:, 1])
        - (pC[:, 0] - pA[:, 0]) * (pB[:, 1] - pA[:, 1])
    )
    mean_depth = (
        depth_flat[simplices[:, 0]] + depth_flat[simplices[:, 1]] + depth_flat[simplices[:, 2]]
    ) / 3.0

    return float(np.sum(area_2d * mean_depth))


def evaluate_all_methods(xs: np.ndarray, ys: np.ndarray, depth_grid: np.ndarray) -> Dict[str, float]:
    """Run all five volumetric estimation methods on a structured grid and return m^3 dictionary."""
    XX, YY = np.meshgrid(xs, ys)
    return {
        "End-Area (Trapezoidal)": end_area_volume(xs, ys, depth_grid),
        "Prismoidal (Simpson 1/3)": prismoidal_volume(xs, ys, depth_grid),
        "Borrow-Pit Grid (1-2-3-4)": borrow_pit_grid_volume(xs, ys, depth_grid),
        "Contour Slicing (1.0m CI)": contour_slicing_volume(xs, ys, depth_grid, contour_interval=1.0),
        "Delaunay TIN Surface": delaunay_tin_volume(XX.ravel(), YY.ravel(), depth_grid.ravel()),
    }
