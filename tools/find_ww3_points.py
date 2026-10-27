"""List GFS-Wave (WW3) point-output stations inside a lat/lon box.

Usage: python tools/find_ww3_points.py bathy/wave_gfs_buoys.txt [lat_min lat_max lon_min lon_max]
The station file is NOAA-EMC global-workflow parm/wave/wave_gfs.buoys.full.
"""
import re, sys

box = [float(x) for x in sys.argv[2:6]] or [13, 17.5, -100, -93]
pts = []
for line in open(sys.argv[1], encoding="utf-8-sig"):
    if line.startswith("$") or not line.strip():
        continue
    m = re.match(r"\s*(-?[\d.]+)\s+(-?[\d.]+)\s+'([^']+)'", line)
    if not m:
        continue
    lon, lat, name = float(m[1]), float(m[2]), m[3].strip()
    if lon > 180:
        lon -= 360
    if box[0] <= lat <= box[1] and box[2] <= lon <= box[3]:
        pts.append((lat, lon, name))
for p in sorted(pts):
    print(p)
print(len(pts), "stations")
