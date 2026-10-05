#!/usr/bin/env python3
"""Mesh-independence study of the blade case.

Meshes meshes/blade.geo at several sizes, runs blade-fem on each and reads the
peaks back from the HDF5 results. Writes results/mesh_study/mesh_study.csv and
prints a Markdown table.

    python3 scripts/blade_mesh_study.py --bin build/blade-fem
"""

import argparse
import csv
import subprocess
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SIZES = [0.0016, 0.0008, 0.0004, 0.0002]  # m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=str(ROOT / "build" / "blade-fem"))
    ap.add_argument("--case", default=str(ROOT / "cases" / "blade.json"))
    ap.add_argument("--out", default=str(ROOT / "results" / "mesh_study"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for lc in SIZES:
        tag = f"lc{lc * 1e3:.1f}mm"
        msh = out / f"blade_{tag}.msh"
        subprocess.run(["gmsh", str(ROOT / "meshes" / "blade.geo"), "-2", "-format", "msh22",
                        "-setnumber", "lc", str(lc), "-o", str(msh)], check=True, capture_output=True)
        base = out / tag / "blade"
        subprocess.run([args.bin, args.case, "--mesh", str(msh), "--output", str(base)],
                       check=True, capture_output=True)
        with h5py.File(f"{base}.h5", "r") as f:
            T = f["fields/temperature"][()]
            szz = f["fields/stress_zz"][()]
            vm = f["fields/von_mises"][()]
            rows.append({
                "lc_mm": lc * 1e3,
                "nodes": f["mesh/nodes"].shape[0],
                "triangles": f["mesh/triangles"].shape[0],
                "T_max_K": float(T.max()),
                "T_min_K": float(T.min()),
                "szz_min_MPa": float(szz.min()) / 1e6,
                "szz_max_MPa": float(szz.max()) / 1e6,
                "von_mises_max_MPa": float(vm.max()) / 1e6,
            })

    with open(out / "mesh_study.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print("| lc [mm] | nodes | T min / max [K] | sigma_zz min / max [MPa] | von Mises max [MPa] |")
    print("|---:|---:|---:|---:|---:|")
    for r in rows:
        print(f"| {r['lc_mm']:.1f} | {r['nodes']} | {r['T_min_K']:.1f} / {r['T_max_K']:.1f} "
              f"| {r['szz_min_MPa']:.0f} / {r['szz_max_MPa']:.0f} | {r['von_mises_max_MPa']:.0f} |")


if __name__ == "__main__":
    main()
