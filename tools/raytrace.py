"""Backward wave-ray tracing over GEBCO bathymetry -> per-spot refraction tables.

For each spot and wave period T, rays start at a nearshore point (~10 m depth) off the
break and are traced offshore over the real depth grid until they reach deep water
(h > L0/2) or hit land. Wave action is conserved along rays, so the nearshore energy that
a deep-water swell from direction D delivers is

    Kr(T, D)^2 = (c0 / cn) * sum_rays G(theta0(ray) - D) * dtheta_n

with G the deep-water directional spread (normalised Gaussian). Shoaling is left to the
Komar & Gaughan breaking formula in the app, so Kr is pure refraction/focusing/sheltering
(straight parallel contours give Kr^2 = cos(theta0)/cos(theta_n), the textbook result).

Usage:  python tools/raytrace.py   (writes bathy/refraction.json and prints a summary)
"""
import json, math, pathlib
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

ROOT = pathlib.Path(__file__).resolve().parent.parent
ASC = next((ROOT / "bathy" / "gebco_2026").glob("*.asc"))
OUT = ROOT / "bathy" / "refraction.json"
OUT_WW3 = ROOT / "bathy" / "refraction_ww3.json"   # Kr^2 per GFS-Wave spectral bin, used by tools/fetch_spectra.py
G = 9.81

# Beach coordinates, same as the app's DEFAULT_SPOTS.
SPOTS = {
    "zicatela":   (15.857, -97.062),
    "lapunta":    (15.842, -97.044),
    "puertecito": (15.745, -96.852),
    "plataforma": (15.762, -96.905),
    "barra":      (15.834, -95.952),
}
START_DEPTH = 10.0                                   # m, seaward of the surf zone
PERIODS = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22]
DIRS = list(range(120, 302, 2))                      # deep-water "from" bearing, deg
D_THETA = 0.25                                       # ray fan spacing, deg
DS = 40.0                                            # ray step, m


def spread_deg(T):
    """Assumed deep-water directional spread: wind sea broad, long-period groundswell narrow."""
    return float(np.clip(35 - 2.0 * (T - 6), 8, 35))


def load_grid():
    hdr = {}
    with open(ASC) as f:
        for _ in range(6):
            k, v = f.readline().split()
            hdr[k.lower()] = float(v)
    z = np.loadtxt(ASC, skiprows=6)[::-1]            # row 0 = south
    nrows, ncols, cs = int(hdr["nrows"]), int(hdr["ncols"]), hdr["cellsize"]
    lat = hdr["yllcorner"] + (np.arange(nrows) + 0.5) * cs
    lon = hdr["xllcorner"] + (np.arange(ncols) + 0.5) * cs
    h = gaussian_filter(-z, 0.8)                     # depth (m, +down), light smoothing vs. caustics
    return h, lat, lon, cs


def wavelength(h, T):
    """Fenton & McKee (1990) explicit dispersion approximation."""
    L0 = G * T * T / (2 * math.pi)
    return L0 * np.tanh((2 * math.pi * np.sqrt(np.maximum(h, 0.05) / G) / T) ** 1.5) ** (2 / 3)


WW3_FREQS = [0.035 * 1.07 ** k for k in range(50)]  # GFS-Wave point-output frequency grid
WW3_DIRS = list(range(5, 360, 10))                  # GFS-Wave 10 deg bins, as "from" bearings


def trace_spot(name, beach, h, lat, lon, cs, periods=PERIODS):
    """Trace the ray fan for each period. Returns start info and, per period, (from bearings of rays that
    reach deep water, c0/cn)."""
    lat0 = lat.mean()
    mx = 111320.0 * math.cos(math.radians(lat0)) * cs   # metres per grid cell, x
    my = 110574.0 * cs                                  # metres per grid cell, y

    # start cell: nearest water cell at >= START_DEPTH
    LA, LO = np.meshgrid(lat, lon, indexing="ij")
    dist = np.hypot((LA - beach[0]) * 110574.0, (LO - beach[1]) * 111320.0 * math.cos(math.radians(lat0)))
    dist[h < START_DEPTH] = np.inf
    si, sj = np.unravel_index(np.argmin(dist), dist.shape)
    start = dict(lat=round(float(lat[si]), 4), lon=round(float(lon[sj]), 4),
                 depth=round(float(h[si, sj]), 1), km_from_beach=round(float(dist[si, sj]) / 1000, 2))

    thetas = np.radians(np.arange(0, 360, D_THETA))   # math angle (0 = east, ccw), offshore-going
    traced = []
    for T in periods:
        c = wavelength(h, T) / T
        cy, cx = np.gradient(c, my, mx)                # dc/dy, dc/dx in 1/s
        L0 = G * T * T / (2 * math.pi)
        c0 = L0 / T
        sample = lambda F, i, j: map_coordinates(F, [i, j], order=1, mode="nearest")
        i = np.full(thetas.size, float(si)); j = np.full(thetas.size, float(sj))
        th = thetas.copy()
        alive = np.ones(thetas.size, bool); deep = np.zeros(thetas.size, bool)
        cn = float(c[si, sj])

        def rhs(i, j, th):
            cc = sample(c, i, j)
            return (np.sin(th) * sample(cx, i, j) - np.cos(th) * sample(cy, i, j)) / cc

        for _ in range(8000):
            if not alive.any():
                break
            a = alive
            k1 = rhs(i[a], j[a], th[a])                                  # RK2 (midpoint)
            im = i[a] + 0.5 * DS * np.sin(th[a]) / my
            jm = j[a] + 0.5 * DS * np.cos(th[a]) / mx
            tm = th[a] + 0.5 * DS * k1
            k2 = rhs(im, jm, tm)
            i[a] += DS * np.sin(tm) / my
            j[a] += DS * np.cos(tm) / mx
            th[a] += DS * k2
            out = (i < 0) | (i > h.shape[0] - 1) | (j < 0) | (j > h.shape[1] - 1)
            hd = sample(h, np.clip(i, 0, h.shape[0] - 1), np.clip(j, 0, h.shape[1] - 1))
            reached = alive & ~out & (hd > L0 / 2)
            # leaving through the south edge counts as deep water (it's the trench); other edges / land block
            reached |= alive & (i < 0)
            deep |= reached
            alive &= ~reached & ~out & (hd > 0.5)

        from_bearing = (90 - np.degrees(th[deep])) % 360  # backward ray heads toward where the swell comes from
        traced.append((from_bearing, c0 / cn))
    return start, traced


def spread_table(traced):
    """Kr for a deep-water swell with the assumed Gaussian directional spread (used with model partitions)."""
    table = []
    for T, (from_bearing, ratio) in zip(PERIODS, traced):
        sig = math.radians(spread_deg(T))
        row = []
        for D in DIRS:
            d = np.radians((from_bearing - D + 180) % 360 - 180)
            e = np.sum(np.exp(-0.5 * (d / sig) ** 2)) / (sig * math.sqrt(2 * math.pi)) * math.radians(D_THETA)
            row.append(round(math.sqrt(ratio * e), 3))
        table.append(row)
    return table


def bin_table(traced, width=10.0):
    """Kr^2 for energy spread evenly across each 10 deg WW3 direction bin (used with full spectra:
    the spectrum already carries the real directional spread, so no assumed spread here)."""
    table = []
    for from_bearing, ratio in traced:
        row = []
        for D in WW3_DIRS:
            d = (from_bearing - D + 180) % 360 - 180
            n = np.count_nonzero(np.abs(d) < width / 2)
            row.append(round(ratio * n * D_THETA / width, 4))
        table.append(row)
    return table


def main():
    h, lat, lon, cs = load_grid()
    result = {"source": "GEBCO 2026 15-arc-second grid, " + ASC.name, "start_depth_m": START_DEPTH,
              "periods": PERIODS, "dirs": DIRS, "spots": {}}
    ww3 = {"freqs": [round(f, 5) for f in WW3_FREQS], "dirs": WW3_DIRS, "spots": {}}
    for name, beach in SPOTS.items():
        start, traced = trace_spot(name, beach, h, lat, lon, cs)
        table = spread_table(traced)
        result["spots"][name] = {"start": start, "K": table}
        _, traced_ww3 = trace_spot(name, beach, h, lat, lon, cs, periods=[1 / f for f in WW3_FREQS])
        ww3["spots"][name] = bin_table(traced_ww3)
        K = np.array(table)
        best = [DIRS[int(np.argmax(r))] for r in K]
        print(f"{name:11s} start {start}")
        for T in (8, 12, 16, 20):
            r = K[PERIODS.index(T)]
            print(f"   T={T:2d}s  Kr@180={r[DIRS.index(180)]:.2f} @200={r[DIRS.index(200)]:.2f} "
                  f"@220={r[DIRS.index(220)]:.2f} @240={r[DIRS.index(240)]:.2f}  max {r.max():.2f} from {best[PERIODS.index(T)]}°")
    OUT.write_text(json.dumps(result, separators=(",", ":")))
    print("wrote", OUT, OUT.stat().st_size, "bytes")
    OUT_WW3.write_text(json.dumps(ww3, separators=(",", ":")))
    print("wrote", OUT_WW3, OUT_WW3.stat().st_size, "bytes")

    # embed into the app so it keeps working as a single file
    app = ROOT / "oaxaca-surf.html"
    html = app.read_text(encoding="utf-8")
    a, b = html.index("// BEGIN REFRACTION\n"), html.index("// END REFRACTION")
    html = html[:a] + "// BEGIN REFRACTION\nconst REFRACTION = " + json.dumps(result, separators=(",", ":")) + ";\n" + html[b:]
    app.write_text(html, encoding="utf-8")
    print("embedded into", app.name)


if __name__ == "__main__":
    main()
