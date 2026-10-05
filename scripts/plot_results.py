#!/usr/bin/env python3
"""Figures for the README, read straight from the solver's HDF5 and CSV output.

    python3 scripts/plot_results.py

Inputs:  results/blade/blade.h5, results/blade_te/blade_te.h5, results/verification/*.csv
Outputs: docs/img/*.png
"""

import csv
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
IMG = ROOT / "docs" / "img"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.dpi": 160, "savefig.bbox": "tight"})


def load(path):
    with h5py.File(path, "r") as f:
        xy = f["mesh/nodes"][()] * 1e3  # mm
        tri = mtri.Triangulation(xy[:, 0], xy[:, 1], f["mesh/triangles"][()])
        fields = {k: f["fields"][k][()] for k in f["fields"]}
    return tri, fields


def blade_axes(ax, title):
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")


def plot_field(ax, tri, values, title, cmap, label, norm=None, vmin=None, vmax=None):
    pc = ax.tripcolor(tri, values, shading="gouraud", cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
    blade_axes(ax, title)
    plt.colorbar(pc, ax=ax, label=label, shrink=0.8, pad=0.02)
    return pc


def mark(ax, tri, values, fmt, use_max=True):
    i = int(np.argmax(values) if use_max else np.argmin(values))
    x, y = tri.x[i], tri.y[i]
    # Put the label on the inboard side so it never runs into the axis labels.
    offset = (25, -22) if x < 0.5 * (tri.x.min() + tri.x.max()) else (-70, 20)
    ax.plot(x, y, "o", ms=5, mfc="none", mec="k")
    ax.annotate(fmt.format(values[i]), (x, y), xytext=offset, textcoords="offset points", fontsize=8,
                arrowprops={"arrowstyle": "-", "lw": 0.6})


def blade_figures():
    tri, f = load(ROOT / "results" / "blade" / "blade.h5")
    T = f["temperature"]
    szz = f["stress_zz_nodal"] / 1e6

    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.triplot(tri, lw=0.25, color="0.35")
    blade_axes(ax, f"Mesh: {tri.x.size} nodes, {tri.triangles.shape[0]} linear triangles")
    fig.savefig(IMG / "blade_mesh.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8, 6.4))
    plot_field(axes[0], tri, T, "Temperature", "inferno", "T [K]")
    mark(axes[0], tri, T, "{:.0f} K")
    lim = np.abs(szz).max()
    plot_field(axes[1], tri, szz, "Axial thermal stress (generalized plane strain)", "RdBu_r", r"$\sigma_{zz}$ [MPa]",
               norm=TwoSlopeNorm(0.0, -lim, lim))
    mark(axes[1], tri, szz, "{:.0f} MPa", use_max=False)
    fig.tight_layout()
    fig.savefig(IMG / "blade_results.png")
    plt.close(fig)


def variant_figure():
    cases = [("3 passages (baseline)", ROOT / "results" / "blade" / "blade.h5"),
             ("+ trailing-edge passage", ROOT / "results" / "blade_te" / "blade_te.h5")]
    data = [(name, *load(path)) for name, path in cases]
    tmin = min(f["temperature"].min() for _, _, f in data)
    tmax = max(f["temperature"].max() for _, _, f in data)
    lim = max(np.abs(f["stress_zz_nodal"]).max() for _, _, f in data) / 1e6

    fig, axes = plt.subplots(2, 2, figsize=(12, 5.2))
    for row, (name, tri, f) in enumerate(data):
        plot_field(axes[row, 0], tri, f["temperature"], f"{name}: temperature", "inferno", "T [K]",
                   vmin=tmin, vmax=tmax)
        mark(axes[row, 0], tri, f["temperature"], "{:.0f} K")
        szz = f["stress_zz_nodal"] / 1e6
        plot_field(axes[row, 1], tri, szz, rf"{name}: $\sigma_{{zz}}$", "RdBu_r", "MPa",
                   norm=TwoSlopeNorm(0.0, -lim, lim))
        mark(axes[row, 1], tri, szz, "{:.0f} MPa", use_max=False)
    fig.tight_layout()
    fig.savefig(IMG / "blade_variants.png")
    plt.close(fig)


def read_csv(name):
    with open(ROOT / "results" / "verification" / f"{name}.csv") as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def slope_guide(ax, h, e, order, label):
    """Reference line of the given order through the finest point, offset below the data."""
    ref = e[-1] * 0.5 * (h / h[-1]) ** order
    ax.loglog(h, ref, "k--", lw=0.8)
    ax.text(h[0], ref[0], f" {label}", fontsize=8, va="top")


def convergence_figure():
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.4))

    d = read_csv("heat_manufactured")
    ax = axes[0]
    ax.loglog(d["h"], d["L2_error"], "o-", label=r"$\|T-T_h\|_{L^2}$")
    ax.loglog(d["h"], d["H1_error"], "s-", label=r"$|T-T_h|_{H^1}$")
    slope_guide(ax, d["h"], d["L2_error"], 2, "slope 2")
    slope_guide(ax, d["h"], d["H1_error"], 1, "slope 1")
    ax.set_title("Heat conduction (manufactured)")

    d = read_csv("cantilever")
    ax = axes[1]
    ax.loglog(d["h"], d["disp_L2_rel"], "o-", label="displacement, rel. $L^2$")
    ax.loglog(d["h"], d["tip_rel_err"], "s-", label="tip deflection, rel.")
    slope_guide(ax, d["h"], d["disp_L2_rel"], 2, "slope 2")
    ax.set_title("Timoshenko cantilever")

    d = read_csv("kirsch")
    ax = axes[2]
    ax.loglog(d["h"], d["stress_L2_rel"], "o-", label="stress, rel. $L^2$")
    ax.loglog(d["h"], d["scf_rel_err"], "s-", label="SCF error")
    slope_guide(ax, d["h"], d["stress_L2_rel"], 1, "slope 1")
    ax.set_title(f"Kirsch hole, SCF {d['scf'][-1]:.3f} (exact 3)")

    d = read_csv("cylinder_gps")
    ax = axes[3]
    ax.loglog(d["h"], d["stress_L2_rel"], "o-", label="stress (4 comp.), rel. $L^2$")
    ax.loglog(d["h"], d["szz_inner_err"], "s-", label=r"$\sigma_{zz}$ at bore")
    slope_guide(ax, d["h"], d["stress_L2_rel"], 1, "slope 1")
    ax.set_title("Heated long cylinder (coupled, GPS)")

    for ax in axes:
        ax.set_xlabel("h (max edge length)")
        ax.grid(True, which="both", lw=0.3, alpha=0.5)
        ax.legend(fontsize=7, loc="lower right")
    axes[0].set_ylabel("error")
    fig.tight_layout()
    fig.savefig(IMG / "convergence.png")
    plt.close(fig)


def main():
    IMG.mkdir(parents=True, exist_ok=True)
    convergence_figure()
    blade_figures()
    variant_figure()
    print(f"figures written to {IMG}")


if __name__ == "__main__":
    main()
