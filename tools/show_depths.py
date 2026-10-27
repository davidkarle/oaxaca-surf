"""Print depths (m) around a point from an ESRI ASCII grid, to inspect what the ray tracer sees near a spot.

Usage: python tools/show_depths.py lat lon [half_size_cells] [step_cells] [grid.asc]
Default grid: the GEBCO 2026 subset in bathy/gebco_2026.
"""
import pathlib, sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
asc = pathlib.Path(sys.argv[5]) if len(sys.argv) > 5 else next((ROOT / "bathy" / "gebco_2026").glob("*.asc"))
hdr = {}
with open(asc) as f:
    for _ in range(6):
        k, v = f.readline().split()
        hdr[k.lower()] = float(v)
z = np.loadtxt(asc, skiprows=6)[::-1]
z[z == hdr.get("nodata_value", -32767)] = np.nan
cs = hdr["cellsize"]
lat = hdr["yllcorner"] + (np.arange(int(hdr["nrows"])) + 0.5) * cs
lon = hdr["xllcorner"] + (np.arange(int(hdr["ncols"])) + 0.5) * cs
la, lo = float(sys.argv[1]), float(sys.argv[2])
n = int(sys.argv[3]) if len(sys.argv) > 3 else 8
st = int(sys.argv[4]) if len(sys.argv) > 4 else 1
i, j = int(np.argmin(abs(lat - la))), int(np.argmin(abs(lon - lo)))
cols = range(j - n * st, j + n * st + 1, st)
print(f"{asc.name}, cell {cs * 111:.2f} km, every {st} cell(s)")
print("        " + "".join(f"{lon[c]:8.3f}" for c in cols))
for r in range(i + n * st, i - n * st - 1, -st):
    print(f"{lat[r]:7.3f} " + "".join(f"{-z[r, c]:8.0f}" if z[r, c] < 0 else f"{'land':>8s}" for c in cols))
