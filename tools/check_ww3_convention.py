"""Check WW3 spec direction convention/units against Open-Meteo's GFS-Wave fields at the same point."""
import json, sys, urllib.request
import numpy as np
from ww3_spec import read_spec, band_widths

path = sys.argv[1]
f, d, t, E, m = read_spec(path)
lat, lon = m[0]["lat"], m[0]["lon"]
url = (f"https://marine-api.open-meteo.com/v1/marine?latitude={lat}&longitude={lon}&models=ncep_gfswave025&timezone=GMT"
       "&hourly=wave_height,wave_direction,wave_period,swell_wave_height,swell_wave_direction,swell_wave_period&past_days=1&forecast_days=3")
om = json.load(urllib.request.urlopen(url))["hourly"]
df, dth = band_widths(f), 2 * np.pi / len(d)
deg = np.degrees(d)
conv = {"from=deg": deg % 360, "from=deg+180": (deg + 180) % 360,
        "from=270-deg": (270 - deg) % 360, "from=90-deg": (90 - deg) % 360}

def mean_dir(Ef, dirs_from):
    w = Ef.sum(axis=1)                                     # energy per direction
    a = np.radians(dirs_from)
    return np.degrees(np.arctan2((w * np.sin(a)).sum(), (w * np.cos(a)).sum())) % 360

for i in range(0, 24, 6):
    ti = t[i]
    if ti not in om["time"]:
        continue
    j = om["time"].index(ti)
    Ei = E[i] * df[None, :] * dth
    hs = 4 * np.sqrt(Ei.sum())
    sw = Ei[:, f < 0.11]                                   # swell band, T > ~9 s
    fp = f[np.argmax(Ei.sum(axis=0))]
    print(f"{ti}  spec Hs {hs:.2f} vs OM {om['wave_height'][j]:.2f} | Tp {1/fp:.1f}s vs OM Tm {om['wave_period'][j]} "
          f"swellT {om['swell_wave_period'][j]} | OM mean dir {om['wave_direction'][j]} swell dir {om['swell_wave_direction'][j]}")
    print("    mean dir (all):  " + "  ".join(f"{k} {mean_dir(Ei, v):.0f}" for k, v in conv.items()))
    print("    mean dir (swell):" + "  ".join(f"{k} {mean_dir(sw, v):.0f}" for k, v in conv.items()))
