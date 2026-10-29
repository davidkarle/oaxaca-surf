# Oaxaca Surf

Breaking-wave forecast for Zicatela, La Punta, Puertecito, Plataforma and Barra de la Cruz.

- **Wave models:** ECMWF WAM, GFS-Wave, MF MFWAM and DWD GWAM via Open-Meteo, plus NOAA's full GFS-Wave 2D spectra.
- **Bathymetry:** refraction per spot from backward ray tracing over the GEBCO 2026 depth grid.
- **Breaking:** Komar & Gaughan breaking height.
- **Hurricanes:** parametric hurricane swell (Young 1988) from NHC eastern Pacific advisories.
- **Calibration:** per-spot correction from logged sessions: your own (stored in your browser) plus the shared ones in `reports/reports.json`. Reports can be spoken or typed in plain words ("surfed Zicatela this morning, 5 feet"), or started from a photo or video: its GPS picks the spot, its capture time the hour, and tapping the wave face and a surfer on it measures the size.

## Shared reports

The Share button on a logged report opens a GitHub issue with the report as JSON. Label the issue `approved` (create the label once) and `.github/workflows/add-report.yml` validates it with `tools/add_report.py`, adds it to `reports/reports.json`, closes the issue and redeploys the site. Only the repository owner's label triggers it.

The page explains the method in detail under "How the prediction works".

## Run locally

```
python tools/serve.py 8765
```

Open http://localhost:8765/oaxaca-surf.html. The server refreshes `spectra/latest.json` and `storms/latest.json` hourly. Opening the HTML file directly also works, but without the spectral and hurricane layers.

## Hosting

`.github/workflows/pages.yml` rebuilds the data hourly and deploys to GitHub Pages (Settings → Pages → Source: GitHub Actions).
GitHub pauses scheduled workflows after 60 days with no repository activity: re-enable it in the Actions tab if that happens.

## Tools

| Script | What it does |
|---|---|
| `tools/raytrace.py` | GEBCO ray tracing → `bathy/refraction*.json`, embedded into the HTML. Rerun if a spot moves. |
| `tools/fetch_spectra.py` | GFS-Wave point spectra from NOMADS → per-spot transformed spectra |
| `tools/fetch_storms.py` | NHC advisories → parametric hurricane swell; `--replay <advisory URL> --shift-to-now` to test on an archived storm |
| `tools/add_report.py` | Validates a shared report from an issue body and adds it to `reports/reports.json` |
| `tools/ww3_spec.py`, `find_ww3_points.py`, `check_ww3_convention.py` | WAVEWATCH III spectra reader and checks |

Data: Open-Meteo, NOAA NCEP/NOMADS, NOAA NHC, GEBCO Compilation Group (2026) GEBCO 2026 Grid.
