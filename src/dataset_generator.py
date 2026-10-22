"""
Synthetic Engineering Survey Dataset Generator — Abuja Central Area (UTM Zone 32N)

Simulates a realistic pre- and post-excavation total-station & RTK-GNSS survey
over a 60 m x 50 m commercial deep basement excavation site in Central Area,
Abuja, F.C.T. (UTM Zone 32N coordinates), featuring:
  - Gently sloping natural Abuja terrain (~476.4 - 478.1 m orthometric height)
  - Battered side slopes meeting the perimeter crest at a realistic non-zero angle
  - Undulating weathered basement-complex rock floor (~11.6 m peak cut depth)
  - Exact continuous 2D quadrature ground-truth reference volume via scipy.integrate.dblquad
"""

import numpy as np
import pandas as pd
from scipy.integrate import dblquad

# Site origin in UTM Zone 32N (Central Area, Abuja, F.C.T.)
ORIGIN_EASTING = 332400.0   # meters East
ORIGIN_NORTHING = 1001800.0 # meters North
LX = 60.0                   # East-West extent (m)
LY = 50.0                   # North-South extent (m)


def natural_ground_surface(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pre-excavation natural topographic elevation Z_pre(x, y) in meters."""
    return (
        476.40
        + 0.018 * x
        + 0.012 * y
        + 0.22 * np.sin(np.pi * x / LX) * np.cos(np.pi * y / LY)
    )


def cut_depth_field(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Continuous deep-excavation cut depth D(x, y) >= 0 across [0, LX] x [0, LY]
    with realistic non-zero batter slope at the excavation rim (non-periodic boundary).
    """
    u = (2.0 * x - LX) / LX   # in [-1, 1]
    v = (2.0 * y - LY) / LY   # in [-1, 1]

    # Quartic flat-bottom / battered-slope envelope vanishing at |u|=1, |v|=1 with non-zero slope
    px = np.maximum(0.0, 1.0 - u ** 4)
    py = np.maximum(0.0, 1.0 - v ** 4)

    # Gentle asymmetric geologic ramp + weathered Abuja basement undulation
    asym = 1.0 + 0.14 * u - 0.10 * v + 0.08 * (u ** 2 - v ** 2)
    undulation = 0.55 * np.cos(1.5 * np.pi * u) * np.cos(1.5 * np.pi * v) * px * py

    return np.maximum(0.0, 11.20 * px * py * asym + undulation)


def compute_exact_reference_volume() -> float:
    """High-precision adaptive 2D Gauss-Kronrod quadrature ground-truth volume (m^3)."""
    vol, _ = dblquad(
        lambda y, x: float(cut_depth_field(x, y)),
        0.0,
        LX,
        lambda _: 0.0,
        lambda _: LY,
        epsabs=1e-6,
        epsrel=1e-8,
    )
    return vol


def generate_survey_grid(spacing: float = 2.5) -> pd.DataFrame:
    """
    Generate structured survey grid points at interval `spacing` (meters).
    Default 2.5 m spacing yields a 25 x 21 = 525-point survey grid (even number of
    intervals along both Easting and Northing so Composite Simpson's 1/3 applies cleanly).
    """
    nx = int(round(LX / spacing)) + 1
    ny = int(round(LY / spacing)) + 1
    xs = np.linspace(0.0, LX, nx)
    ys = np.linspace(0.0, LY, ny)
    XX, YY = np.meshgrid(xs, ys)

    z_pre = natural_ground_surface(XX, YY)
    depth = cut_depth_field(XX, YY)
    z_post = z_pre - depth

    rows = []
    pt_id = 1
    for j in range(ny):
        for i in range(nx):
            rows.append({
                "Point_ID": f"ABJ-{pt_id:04d}",
                "Grid_Row": j,
                "Grid_Col": i,
                "Local_X_m": round(float(XX[j, i]), 3),
                "Local_Y_m": round(float(YY[j, i]), 3),
                "Easting_m": round(ORIGIN_EASTING + float(XX[j, i]), 3),
                "Northing_m": round(ORIGIN_NORTHING + float(YY[j, i]), 3),
                "Pre_Elev_m": round(float(z_pre[j, i]), 4),
                "Post_Elev_m": round(float(z_post[j, i]), 4),
                "Cut_Depth_m": round(float(depth[j, i]), 4),
            })
            pt_id += 1
    return pd.DataFrame(rows)
