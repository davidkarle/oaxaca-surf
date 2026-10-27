"""Serve the app locally and keep spectra/latest.json (GFS-Wave) and storms/latest.json (NHC) fresh, hourly.

Usage: python tools/serve.py [port]    then open http://localhost:<port>/oaxaca-surf.html
"""
import functools, http.server, pathlib, sys, threading, time, traceback

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import fetch_spectra, fetch_storms

ROOT = pathlib.Path(__file__).resolve().parent.parent


def refresher():
    while True:
        for job in (fetch_spectra.main, fetch_storms.main):
            try:
                job()
            except Exception:
                traceback.print_exc()
        time.sleep(3600)


if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    threading.Thread(target=refresher, daemon=True).start()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT))
    print(f"serving http://localhost:{port}/oaxaca-surf.html", flush=True)
    http.server.ThreadingHTTPServer(("127.0.0.1", port), handler).serve_forever()
