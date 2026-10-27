"""Fetch GFS-Wave (WAVEWATCH III) 2D spectra for the stations off Oaxaca and transform them to each spot.

Every frequency/direction bin of the deep-water spectrum is multiplied by that spot's GEBCO ray-traced
refraction factor Kr^2 (bathy/refraction_ww3.json), giving the energy that actually arrives. The result
per spot and hour is an equivalent deep-water height H0 = 4*sqrt(m0') and a breaking-weighted period
Te = (sum E' T^0.8 / sum E')^(1/0.8) (Komar & Gaughan: Hb^2 ~ H0^1.6 T^0.8), which the app turns into
breaking height. Writes spectra/latest.json, which the app loads when served over http.

Usage: python tools/fetch_spectra.py
"""
import datetime as dt, json, math, pathlib, urllib.request
import numpy as np
from ww3_spec import read_spec, band_widths

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "spectra" / "raw"
OUT = ROOT / "spectra" / "latest.json"
NOMADS = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
# Nearest GFS-Wave point-output stations (from tools/find_ww3_points.py), both in deep water.
STATIONS = ["RW-NH2-76", "TPC54"]
# Deep-water point per spot, same as the app's wlat/wlon.
SPOTS = {"zicatela": (15.55, -97.10), "lapunta": (15.55, -97.10), "puertecito": (15.60, -96.86),
         "plataforma": (15.60, -96.90), "barra": (15.55, -95.95)}
LOCAL_UTC_OFFSET = -6          # Oaxaca, no DST since 2022


def exists(url):
    try:
        urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=30)
        return True
    except Exception:
        return False


def latest_cycle():
    now = dt.datetime.now(dt.timezone.utc)
    for back in range(0, 48, 6):
        t = now - dt.timedelta(hours=back)
        cyc = f"{t:%Y%m%d}", f"{t.hour // 6 * 6:02d}"
        base = f"{NOMADS}/gfs.{cyc[0]}/{cyc[1]}/wave/station/bulls.t{cyc[1]}z"
        if all(exists(f"{base}/gfswave.{s}.spec") for s in STATIONS):
            return cyc, base
    raise RuntimeError("no complete GFS-Wave cycle found on NOMADS in the last 48 h")


def km(a, b):
    return math.hypot((a[0] - b[0]) * 110.6, (a[1] - b[1]) * 111.3 * math.cos(math.radians(a[0])))


def main(force=False):
    cyc, base = latest_cycle()
    tag = f"{cyc[0]}_{cyc[1]}z"
    if not force and OUT.exists() and json.loads(OUT.read_text()).get("cycle") == tag:
        print("already have", tag)
        return False
    RAW.mkdir(parents=True, exist_ok=True)
    data = {}
    for s in STATIONS:
        p = RAW / f"{s}.spec"
        urllib.request.urlretrieve(f"{base}/gfswave.{s}.spec", p)
        data[s] = read_spec(p)
        print("downloaded", s, p.stat().st_size, "bytes")

    tab = json.loads((ROOT / "bathy" / "refraction_ww3.json").read_text())
    freqs, _, _, _, _ = data[STATIONS[0]]
    assert np.allclose(freqs, tab["freqs"], rtol=6e-3)  # file prints 3 significant figures, "WW3 frequency grid changed: rerun tools/raytrace.py"
    df = band_widths(freqs)
    dth = 2 * math.pi / 36
    T = 1 / freqs
    times = sorted(set.intersection(*(set(d[2]) for d in data.values())))
    # WW3 writes nautical "going to" directions; +180 gives the "from" bearing the tables use (bins at 5,15,..355)
    order = {s: np.array([int(round(((math.degrees(a) + 180) % 360 - 5) / 10)) % 36 for a in data[s][1]]) for s in STATIONS}

    def on_tab_grid(s, i):
        E = np.zeros((36, len(freqs)))
        E[order[s]] = data[s][3][data[s][2].index(times[i])]
        return E

    out = {"source": "NOAA GFS-Wave (WAVEWATCH III) point spectra via NOMADS", "cycle": tag,
           "generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"),
           "stations": [{"name": s, "lat": data[s][4][0]["lat"], "lon": data[s][4][0]["lon"]} for s in STATIONS],
           "spots": {}}
    for spot, pt in SPOTS.items():
        w = {s: 1 / km(pt, (data[s][4][0]["lat"], data[s][4][0]["lon"])) ** 2 for s in STATIONS}
        wsum = sum(w.values())
        K2 = np.array(tab["spots"][spot]).T            # [dir, freq]
        reach = K2 > 0                                  # bins that can reach the spot at all
        cols = {k: [] for k in ("h0", "te", "hs", "tp", "dp")}
        for i in range(len(times)):
            E = sum(w[s] / wsum * on_tab_grid(s, i) for s in STATIONS) * df[None, :] * dth   # m^2 per bin
            Et = E * K2
            m0 = Et.sum()
            cols["h0"].append(4 * math.sqrt(m0))
            cols["te"].append(float((Et.sum(axis=0) * T ** 0.8).sum() / m0) ** 1.25 if m0 > 0 else None)
            Ein = np.where(reach, E, 0)                 # offshore swell headed for the spot
            cols["hs"].append(4 * math.sqrt(Ein.sum()))
            fi = int(np.argmax(Ein.sum(axis=0)))
            cols["tp"].append(float(T[fi]))
            cols["dp"].append(float(tab["dirs"][int(np.argmax(Ein[:, fi]))]))
        # hourly in local time (output is hourly to +120 h then 3-hourly): interpolate
        tu = np.array([dt.datetime.fromisoformat(t).replace(tzinfo=dt.timezone.utc).timestamp() for t in times])
        hours = np.arange(tu[0], tu[-1] + 1, 3600)
        rec = {"t": [(dt.datetime.fromtimestamp(x, dt.timezone.utc) + dt.timedelta(hours=LOCAL_UTC_OFFSET)).strftime("%Y-%m-%dT%H:00") for x in hours]}
        for k, v in cols.items():
            v = np.array([np.nan if x is None else x for x in v])
            if k == "dp":   # nearest, not interpolated, for directions
                vi = v[np.clip(np.searchsorted(tu, hours), 0, len(tu) - 1)]
            else:
                vi = np.interp(hours, tu, v)
            rec[k] = [None if np.isnan(x) else round(float(x), 3 if k in ("h0", "hs") else 1) for x in vi]
        out["spots"][spot] = rec
        j = 0
        print(f"{spot:11s} {rec['t'][j]}  H0' {rec['h0'][j]:.2f} m  Te {rec['te'][j]} s | offshore in-window Hs {rec['hs'][j]:.2f} m Tp {rec['tp'][j]} s from {rec['dp'][j]:.0f}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    print("wrote", OUT, OUT.stat().st_size, "bytes, cycle", tag)
    return True


if __name__ == "__main__":
    import sys
    main(force="--force" in sys.argv)
