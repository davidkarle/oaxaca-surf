"""Parametric tropical-cyclone swell for the Oaxaca spots from NHC eastern Pacific forecast advisories.

Global wave models (0.25 deg) under-resolve a hurricane's small, intense core, so swell from storms close
to Mexico is often forecast too small. For every active EP storm this reads NHC's forecast advisory (TCM),
interpolates the track hourly and, at each hour:

  1. Radius of max wind  Rmax = 46.4 exp(-0.0155 Vmax + 0.0169 |lat|) km      (Willoughby & Rahn 2004)
  2. Equivalent fetch    F/R' = a Vmax^2 + b Vmax Vfm + c Vfm^2 + d Vmax + e Vfm + f,
                         R' = 22.5e3 log10(Rmax) - 70.8e3                      (Young 1988, SI units)
     which captures the "trapped fetch" when the storm moves at about the waves' group speed.
  3. Peak waves          g Hs/V^2 = 0.0016 (gF/V^2)^0.5,  Vmax fp/g = 3.5 (gF/V^2)^-0.33   (JONSWAP)
  4. Emission            a moving storm radiates its biggest, longest swell ahead of its track and much less
                         behind: H x (0.75 + 0.25 cos A), Tp x (0.9 + 0.1 cos A), A = angle between storm
                         heading and the direction to the spot (Young 1988's asymmetric wave field, simplified).
     Decay to the spot   H = Hs sqrt(r0/r) (spreading over an arc), with r0 fitted per advisory so the decay
                         reaches 4 m at NHC's mean analysed 4 m-seas radius.
  5. Arrival             t + r/cg, cg = g Tp / 4 pi; direction = bearing from the spot to the storm.

The app runs the resulting deep-water swell through each spot's GEBCO refraction table.
Writes storms/latest.json.

Usage: python tools/fetch_storms.py
       python tools/fetch_storms.py --replay <advisory .shtml URL> [--shift-to-now]   (test on an archived storm)
"""
import datetime as dt, html, json, math, pathlib, re, sys, urllib.request
from fetch_spectra import SPOTS, LOCAL_UTC_OFFSET

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "storms" / "latest.json"
G, KT = 9.81, 0.514444
MONTHS = "JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split()


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "oaxaca-surf (personal forecast)"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")


def parse_tcm(text):
    t = html.unescape(re.sub(r"<[^>]+>", "", text))
    hm, mon, day, yr = re.search(r"(\d{4}) UTC \w{3} (\w{3}) (\d{1,2}) (\d{4})", t).groups()
    issued = dt.datetime(int(yr), MONTHS.index(mon) + 1, int(day), int(hm[:2]), tzinfo=dt.timezone.utc)

    def when(d, h):
        w = issued.replace(day=1) + dt.timedelta(days=int(d) - 1, hours=int(h) - issued.hour)
        return w if int(d) >= issued.day - 15 else (w + dt.timedelta(days=32)).replace(day=int(d))

    ll = lambda a, ns, b, ew: (float(a) * (1 if ns == "N" else -1), float(b) * (-1 if ew == "W" else 1))
    name = re.search(r"\n\s*(?:\w+[- ])*?(?:HURRICANE|TROPICAL STORM|TROPICAL DEPRESSION|POTENTIAL TROPICAL CYCLONE|"
                     r"POST-TROPICAL CYCLONE|REMNANTS OF)\s+(\w[\w-]*)\s+FORECAST/ADVISORY NUMBER\s+(\d+)", t)
    c = re.search(r"CENTER LOCATED NEAR\s+([\d.]+)([NS])\s+([\d.]+)([EW]) AT (\d\d)/(\d\d)(\d\d)Z", t)
    mv = re.search(r"(\d+) DEGREES AT\s+(\d+) KT", t)
    vmax = int(re.search(r"MAX SUSTAINED WINDS\s+(\d+) KT", t)[1])
    seas = re.search(r"4 M SEAS\.+\s*(\d+)NE\s+(\d+)SE\s+(\d+)SW\s+(\d+)NW", t)
    pts = [dict(t=when(c[5], c[6]), pos=ll(*c.groups()[:4]), vmax=vmax, inland=False)]
    for m in re.finditer(r"(?:FORECAST|OUTLOOK) VALID (\d\d)/(\d\d)(\d\d)Z\s+([\d.]+)([NS])\s+([\d.]+)([EW])([^\n]*)\s*\n\s*MAX WIND\s+(\d+) KT", t):
        pts.append(dict(t=when(m[1], m[2]), pos=ll(m[4], m[5], m[6], m[7]), vmax=int(m[9]), inland="INLAND" in m[8]))
    return dict(name=name[1].title() if name else "?", adv=int(name[2]) if name else None, issued=issued,
                motion=(int(mv[1]), int(mv[2])) if mv else (0, 0),
                seas4_nm=[int(x) for x in seas.groups()] if seas else None, track=pts)


def dist_bearing(a, b):
    """Great-circle distance (km) and initial bearing (deg) from a to b."""
    p1, p2, dl = math.radians(a[0]), math.radians(b[0]), math.radians(b[1] - a[1])
    d = 2 * 6371 * math.asin(math.sqrt(math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2))
    brg = math.degrees(math.atan2(math.sin(dl) * math.cos(p2), math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)))
    return d, brg % 360


def young(vmax, vfm, lat):
    """Max Hs (m), Tp (s), Rmax (km), fetch (km) for a moving tropical cyclone."""
    v = min(max(vmax, 17.0), 68.0)                               # model's calibrated range
    vf = min(max(vfm, 0.0), 12.0)
    rmax = 46.4 * math.exp(-0.0155 * v + 0.0169 * abs(lat))
    R = min(max(rmax, 15.0), 60.0) * 1e3
    Rp = 22.5e3 * math.log10(R) - 70.8e3
    F = Rp * (-2.175e-3 * v * v + 1.506e-2 * v * vf - 1.223e-1 * vf * vf + 2.190e-1 * v + 6.737e-1 * vf + 7.980e-1)
    gf = min(G * F / (v * v), 2.3e4)                              # fully developed limit
    return 0.0016 * math.sqrt(gf) * v * v / G, v / (3.5 * G) * gf ** 0.33, rmax, F / 1e3


def hourly_track(adv):
    pts = adv["track"]
    out = []
    for p, q in zip(pts, pts[1:]):
        n = int((q["t"] - p["t"]).total_seconds() // 3600)
        d, brg = dist_bearing(p["pos"], q["pos"])
        vfm = d * 1e3 / max(1, n * 3600)
        for k in range(n):
            f = k / n
            out.append(dict(t=p["t"] + dt.timedelta(hours=k), lat=p["pos"][0] + f * (q["pos"][0] - p["pos"][0]),
                            lon=p["pos"][1] + f * (q["pos"][1] - p["pos"][1]), vmax=(p["vmax"] + f * (q["vmax"] - p["vmax"])) * KT,
                            vfm=vfm, heading=brg, inland=q["inland"] if f >= 0.5 else p["inland"]))
    return out


def storm_swell(adv):
    track = hourly_track(adv)
    if not track:
        return None, {}
    t0 = track[0]
    hs0, tp0, rmax0, f0 = young(t0["vmax"], t0["vfm"], t0["lat"])
    # decay scale r0/Rmax from NHC's mean 4 m-seas radius (fallback 2 Rmax)
    k = 2.0
    if adv["seas4_nm"] and hs0 > 4 and any(adv["seas4_nm"]):
        r4 = sum(adv["seas4_nm"]) / 4 * 1.852
        k = min(6.0, max(1.0, r4 * (4 / hs0) ** 2 / rmax0))
    spots = {}
    for spot, pt in SPOTS.items():
        arr = {}
        for s in track:
            if s["inland"] or s["vmax"] < 17.0:                  # below tropical-storm strength: leave to the models
                continue
            hs, tp, rmax, _ = young(s["vmax"], s["vfm"], s["lat"])
            r, frm = dist_bearing(pt, (s["lat"], s["lon"]))
            ca = math.cos(math.radians((frm + 180) - s["heading"])) if s["vfm"] > 0.5 else 0.0
            hs, tp = hs * (0.75 + 0.25 * ca), tp * (0.9 + 0.1 * ca)
            r0 = k * rmax
            h = hs if r <= r0 else hs * math.sqrt(r0 / r)
            ta = s["t"] + dt.timedelta(seconds=r * 1e3 / (G * tp / (4 * math.pi)))
            key = ta.replace(minute=0, second=0, microsecond=0)
            if h > arr.get(key, {"h": 0})["h"]:
                arr[key] = dict(h=h, tp=tp, d=frm, r=r)
        spots[spot] = arr
    info = dict(name=adv["name"], adv=adv["adv"], issued=adv["issued"].isoformat(), lat=t0["lat"], lon=t0["lon"],
                vmax_kt=round(t0["vmax"] / KT), motion=adv["motion"], hs=round(hs0, 1), tp=round(tp0, 1),
                rmax_km=round(rmax0), fetch_km=round(f0),
                track=[dict(t=p["t"].isoformat(), lat=p["pos"][0], lon=p["pos"][1], vmax=p["vmax"], inland=p["inland"]) for p in adv["track"]])
    return info, spots


def build(advisories, shift=None):
    out = {"source": "NHC eastern Pacific forecast advisories + Young (1988) parametric model",
           "generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="minutes"), "storms": [], "spots": {}}
    merged = {s: {} for s in SPOTS}
    for adv in advisories:
        if shift:
            adv["issued"] += shift
            for p in adv["track"]:
                p["t"] += shift
        info, spots = storm_swell(adv)
        if not info:
            continue
        out["storms"].append(info)
        for s, arr in spots.items():
            for k, v in arr.items():
                if v["h"] > merged[s].get(k, {"h": 0})["h"]:
                    merged[s][k] = dict(v, storm=info["name"])
    for s, arr in merged.items():
        # arrivals bunch up while a storm approaches and spread out as it leaves: fill hourly gaps linearly
        keys = sorted(arr)
        for a, b in zip(keys, keys[1:]):
            n = int((b - a).total_seconds() // 3600)
            for k in range(1, n):
                f, va, vb = k / n, arr[a], arr[b]
                dd = (vb["d"] - va["d"] + 180) % 360 - 180
                arr[a + dt.timedelta(hours=k)] = dict(h=va["h"] + f * (vb["h"] - va["h"]), tp=va["tp"] + f * (vb["tp"] - va["tp"]),
                                                      d=(va["d"] + f * dd) % 360, r=va["r"] + f * (vb["r"] - va["r"]),
                                                      storm=va["storm"] if f < 0.5 else vb["storm"])
        keys = sorted(arr)
        loc = lambda k: (k + dt.timedelta(hours=LOCAL_UTC_OFFSET)).strftime("%Y-%m-%dT%H:00")
        out["spots"][s] = {"t": [loc(k) for k in keys], "h": [round(arr[k]["h"], 2) for k in keys],
                           "tp": [round(arr[k]["tp"], 1) for k in keys], "d": [round(arr[k]["d"]) for k in keys],
                           "km": [round(arr[k]["r"]) for k in keys], "storm": [arr[k]["storm"] for k in keys]}
    return out


def main(argv=()):
    argv = list(argv)
    if "--replay" in argv:
        adv = parse_tcm(get(argv[argv.index("--replay") + 1]))
        shift = None
        if "--shift-to-now" in argv:
            now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
            shift = now - adv["issued"]
        out = build([adv], shift)
        path = ROOT / "storms" / "replay.json"
    else:
        cur = json.loads(get("https://www.nhc.noaa.gov/CurrentStorms.json"))
        advs = []
        for s in cur.get("activeStorms", []):
            if s["id"].startswith("ep") and s.get("forecastAdvisory"):
                try:
                    advs.append(parse_tcm(get(s["forecastAdvisory"]["url"])))
                except Exception as e:
                    print("skip", s["id"], e)
        out = build(advs)
        path = OUT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, separators=(",", ":")))
    for st in out["storms"]:
        print(f"{st['name']} adv {st['adv']}: {st['lat']:.1f},{st['lon']:.1f} {st['vmax_kt']} kt, motion {st['motion']} -> "
              f"Hs {st['hs']} m Tp {st['tp']} s (Rmax {st['rmax_km']} km, fetch {st['fetch_km']} km)")
    for s, v in out["spots"].items():
        if v["h"]:
            i = max(range(len(v["h"])), key=v["h"].__getitem__)
            print(f"  {s:11s} first {v['t'][0]}  peak {v['h'][i]:.1f} m {v['tp'][i]} s from {v['d'][i]} deg "
                  f"({v['km'][i]} km away) at {v['t'][i]} local")
    print("wrote", path, "storms:", len(out["storms"]))
    return out


if __name__ == "__main__":
    main(sys.argv[1:])
