"""
End-to-End Volumetric Comparison & Visualization Pipeline
Author: Shekwomagode Aaron Tonya Solomon
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

from src.dataset_generator import (
    ORIGIN_EASTING,
    ORIGIN_NORTHING,
    compute_exact_reference_volume,
    generate_survey_grid,
)
from src.volumetric_methods import evaluate_all_methods

UNIT_HAULAGE_RATE_USD = 8.50   # Realistic mechanical excavation + haulage rate ($/m^3)
SWELL_FACTOR = 1.22            # Typical weathered Abuja basement/laterite bulking factor


def main() -> None:
    os.makedirs("data", exist_ok=True)
    os.makedirs("outputs", exist_ok=True)

    print("1. Computing exact continuous 2D Gauss-Kronrod reference volume...")
    v_exact = compute_exact_reference_volume()
    print(f"   Exact Reference Cut Volume: {v_exact:,.2f} m^3")

    print("2. Generating 2.5 m RTK-GNSS / Total Station survey grid (Abuja UTM Zone 32N)...")
    df_survey = generate_survey_grid(spacing=2.5)
    df_survey.to_csv("data/abuja_deep_excavation_survey.csv", index=False)
    print(f"   Saved {len(df_survey)} survey stations -> data/abuja_deep_excavation_survey.csv")

    xs = np.sort(df_survey["Local_X_m"].unique())
    ys = np.sort(df_survey["Local_Y_m"].unique())
    depth_grid = df_survey.pivot(index="Local_Y_m", columns="Local_X_m", values="Cut_Depth_m").values
    z_post_grid = df_survey.pivot(index="Local_Y_m", columns="Local_X_m", values="Post_Elev_m").values

    results = evaluate_all_methods(xs, ys, depth_grid)
    summary_rows = []
    for method, vol in results.items():
        diff = vol - v_exact
        pct_err = (diff / v_exact) * 100.0
        loose_vol = vol * SWELL_FACTOR
        cost_usd = vol * UNIT_HAULAGE_RATE_USD
        cost_diff_usd = diff * UNIT_HAULAGE_RATE_USD
        summary_rows.append({
            "Method": method,
            "Bank_Volume_m3": round(vol, 2),
            "Bulked_Haul_Volume_m3": round(loose_vol, 2),
            "Exact_Reference_m3": round(v_exact, 2),
            "Signed_Error_m3": round(diff, 2),
            "Abs_Pct_Error": round(abs(pct_err), 4),
            "Signed_Pct_Error": round(pct_err, 4),
            "Est_Earthwork_Cost_USD": round(cost_usd, 2),
            "Cost_Discrepancy_USD": round(cost_diff_usd, 2),
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv("outputs/method_comparison_summary.csv", index=False)
    print("\n=== Comparative Volumetric Summary (2.5 m Grid Spacing) ===")
    print(df_summary[["Method", "Bank_Volume_m3", "Signed_Error_m3", "Signed_Pct_Error", "Cost_Discrepancy_USD"]].to_string(index=False))

    # 3. Multi-resolution grid convergence experiment (0.5m, 1.0m, 2.0m, 2.5m, 5.0m, 10.0m)
    spacings = [1.0, 2.0, 2.5, 5.0, 7.5, 10.0]
    conv_records = []
    for sp in spacings:
        sub_df = generate_survey_grid(spacing=sp)
        s_xs = np.sort(sub_df["Local_X_m"].unique())
        s_ys = np.sort(sub_df["Local_Y_m"].unique())
        s_d = sub_df.pivot(index="Local_Y_m", columns="Local_X_m", values="Cut_Depth_m").values
        sub_res = evaluate_all_methods(s_xs, s_ys, s_d)
        rec = {"Grid_Spacing_m": sp, "Station_Count": len(sub_df)}
        for m_name, m_vol in sub_res.items():
            rec[m_name] = abs((m_vol - v_exact) / v_exact) * 100.0
        conv_records.append(rec)

    df_conv = pd.DataFrame(conv_records)
    df_conv.to_csv("outputs/grid_convergence_summary.csv", index=False)

    # 4. Generate 4-Panel Publication Dashboard Figure
    plt.rcParams.update({
        "font.family": "sans-serif",
        "axes.edgecolor": "#CBD5E1",
        "axes.linewidth": 0.9,
    })
    fig = plt.figure(figsize=(14.5, 9.6), dpi=180, facecolor="#F8FAFC")
    fig.suptitle(
        "Comparative Analysis of Volumetric Methods for Deep Excavation in Abuja, F.C.T.",
        fontsize=15.5, fontweight="bold", color="#0F172A", y=0.975
    )
    fig.text(
        0.5, 0.942,
        f"UTM Zone 32N Central Area Site (60 m × 50 m | Peak Cut: {depth_grid.max():.2f} m | Exact Reference Volume: {v_exact:,.1f} m³)",
        ha="center", fontsize=10.5, color="#475569"
    )

    XX, YY = np.meshgrid(xs, ys)

    # Panel A: 3D Post-Excavation DEM Surface
    ax1 = fig.add_subplot(2, 2, 1, projection="3d")
    ax1.set_facecolor("#F8FAFC")
    surf = ax1.plot_surface(
        XX, YY, z_post_grid, cmap=cm.viridis, edgecolor="k", linewidth=0.15, alpha=0.94, antialiased=True
    )
    ax1.view_init(elev=30, azim=-132)
    ax1.set_title("A. 3D Post-Excavation Digital Elevation Model", fontsize=11.5, fontweight="bold", color="#16324F", pad=8)
    ax1.set_xlabel("Local Easting ΔE (m)", fontsize=8.5, labelpad=6)
    ax1.set_ylabel("Local Northing ΔN (m)", fontsize=8.5, labelpad=6)
    ax1.set_zlabel("Elevation (m)", fontsize=8.5, labelpad=4)
    ax1.tick_params(labelsize=7.5)
    cb1 = fig.colorbar(surf, ax=ax1, shrink=0.62, pad=0.07)
    cb1.set_label("Elevation (m)", fontsize=8.5)

    # Panel B: Planimetric Cut-Depth Contour Map + Cross-Section Line
    ax2 = fig.add_subplot(2, 2, 2)
    cf = ax2.contourf(xs, ys, depth_grid, levels=12, cmap="YlGnBu")
    cs = ax2.contour(xs, ys, depth_grid, levels=[2, 5, 8, 10, 11], colors="#0F172A", linewidths=0.75)
    ax2.clabel(cs, fmt="%.0f m", fontsize=7.5)
    mid_y_idx = len(ys) // 2
    ax2.axhline(ys[mid_y_idx], color="#DC2626", linestyle="--", linewidth=1.6, label=f"Section A–A' (Y = {ys[mid_y_idx]:.1f} m)")
    ax2.set_title("B. Planimetric Cut-Depth Contours & Survey Grid", fontsize=11.5, fontweight="bold", color="#16324F")
    ax2.set_xlabel("Local Easting X (m)", fontsize=9.5)
    ax2.set_ylabel("Local Northing Y (m)", fontsize=9.5)
    ax2.legend(loc="upper right", fontsize=8.2, framealpha=0.92)
    cb2 = fig.colorbar(cf, ax=ax2)
    cb2.set_label("Cut Depth (m)", fontsize=8.5)

    # Panel C: Signed Volumetric Deviation (m^3) vs Exact Reference
    ax3 = fig.add_subplot(2, 2, 3)
    palette = ["#0284C7", "#059669", "#6366F1", "#D97706", "#7C3AED"]
    methods_short = df_summary["Method"].tolist()
    errs_m3 = df_summary["Signed_Error_m3"].tolist()
    bars = ax3.barh(methods_short, errs_m3, color=palette, height=0.56, edgecolor="#0F172A", linewidth=0.6)
    ax3.axvline(0.0, color="#0F172A", linewidth=1.1, linestyle="--")
    max_neg = abs(min(errs_m3 + [-10.0]))
    max_pos = abs(max(errs_m3 + [10.0]))
    ax3.set_xlim(-max_neg * 1.82, max(max_pos * 1.9, max_neg * 0.35))
    span = max_neg * 1.82
    for bar, val, pct in zip(bars, errs_m3, df_summary["Signed_Pct_Error"]):
        offset = span * 0.03 if val >= 0 else -span * 0.03
        ha = "left" if val >= 0 else "right"
        ax3.text(
            val + offset, bar.get_y() + bar.get_height() / 2.0,
            f"{val:+.1f} m³ ({pct:+.2f}%)",
            va="center", ha=ha, fontsize=8.3, fontweight="bold", color="#0F172A"
        )
    ax3.set_title("C. Signed Volume Error Relative to True Benchmark (2.5 m Grid)", fontsize=11.5, fontweight="bold", color="#16324F")
    ax3.set_xlabel("Volume Discrepancy ΔV (m³)", fontsize=9.5)
    ax3.grid(axis="x", linestyle=":", alpha=0.6)

    # Panel D: Grid Spacing Convergence Sensitivity
    ax4 = fig.add_subplot(2, 2, 4)
    markers = ["o", "s", "^", "D", "P"]
    for m_name, color, mk in zip(methods_short, palette, markers):
        ax4.plot(
            df_conv["Grid_Spacing_m"],
            df_conv[m_name],
            marker=mk, linewidth=2.0, markersize=6, color=color, label=m_name
        )
    ax4.set_title("D. Survey Grid Spacing Sensitivity (1.0 m → 10.0 m Interval)", fontsize=11.5, fontweight="bold", color="#16324F")
    ax4.set_xlabel("Survey Grid Interval Δs (m)", fontsize=9.5)
    ax4.set_ylabel("Absolute Volume Error (%)", fontsize=9.5)
    ax4.grid(True, linestyle=":", alpha=0.65)
    ax4.legend(fontsize=7.8, loc="upper left", framealpha=0.92)

    plt.tight_layout(rect=[0.01, 0.01, 0.99, 0.92])
    out_img = "outputs/volumetric_benchmark_dashboard.png"
    plt.savefig(out_img, dpi=180)
    plt.close(fig)
    print(f"\nSaved publication figure -> {out_img}")


if __name__ == "__main__":
    main()
