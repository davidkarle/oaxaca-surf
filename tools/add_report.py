"""Add a session report shared from the app (a GitHub issue) to reports/reports.json.

The app's Share button opens an issue whose body holds the report as a ```json block. When the owner labels
the issue "approved", .github/workflows/add-report.yml runs this with the issue body in $ISSUE_BODY.
The body is untrusted: the report is rebuilt field by field from checked values, so nothing else gets in.

Usage: ISSUE_BODY="..." python tools/add_report.py      (or pass the body on stdin)
Exit 0 = added or already there (message on stdout), 1 = not a valid report.
"""
import datetime as dt, json, math, os, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FILE = ROOT / "reports" / "reports.json"
SPOTS = {"zicatela", "lapunta", "puertecito", "plataforma", "barra"}
KINDS = {"swell", "swell2", "swell3", "wind"}


class Bad(Exception):
    pass


def num(v, lo, hi, what, null=False):
    if v is None and null:
        return None
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not lo <= v <= hi:
        raise Bad(f"{what} must be a number from {lo} to {hi}")
    return v


def text(v, n, what):
    if not isinstance(v, str) or len(v) > n:
        raise Bad(f"{what} must be text of at most {n} characters")
    return v


def clean(r):
    if not isinstance(r, dict):
        raise Bad("the report must be a JSON object")
    if r.get("spot") not in SPOTS:
        raise Bad(f"unknown spot {r.get('spot')!r}")
    t = text(r.get("t"), 16, "t")
    if not re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:00", t):
        raise Bad("t must look like 2026-10-27T10:00")
    when = dt.datetime.fromisoformat(t)
    if not dt.datetime(2024, 1, 1) <= when <= dt.datetime.now() + dt.timedelta(days=1):
        raise Bad("t is not a plausible session time")
    parts = r.get("parts")
    if not isinstance(parts, list) or not 1 <= len(parts) <= 10:
        raise Bad("parts must list each model's swell partitions (reports without model data don't calibrate)")
    out_parts = []
    for model in parts:
        if not isinstance(model, list) or not 1 <= len(model) <= 6:
            raise Bad("each entry of parts must be a list of 1-6 partitions")
        out = []
        for p in model:
            if not isinstance(p, dict) or p.get("kind") not in KINDS:
                raise Bad("each partition needs h, t, d and a known kind")
            q = {"h": num(p.get("h"), 0, 30, "h"), "t": num(p.get("t"), 0, 30, "t"),
                 "d": num(p.get("d"), 0, 360, "d"), "kind": p["kind"]}
            if p.get("pre") is True:
                q["pre"] = True
            out.append(q)
        out_parts.append(out)
    models = r.get("models")
    if not isinstance(models, list) or len(models) != len(parts):
        raise Bad("models must name one model per entry of parts")
    return {"spot": r["spot"], "t": t, "obs": num(r.get("obs"), 0, 60, "obs"),
            "sl": num(r.get("sl"), 0, 60, "sl", null=True), "note": text(r.get("note", ""), 200, "note"),
            "models": [text(m, 40, "model name") for m in models], "parts": out_parts,
            "tp": num(r.get("tp"), 0, 30, "tp", null=True), "dir": num(r.get("dir"), 0, 360, "dir", null=True)}


def extract(body):
    m = re.search(r"```json\s*(.*?)```", body, re.S)
    if not m:
        raise Bad("no ```json block in the issue")
    try:
        return json.loads(m.group(1))
    except json.JSONDecodeError as e:
        raise Bad(f"the json block doesn't parse: {e}")


def main():
    body = os.environ.get("ISSUE_BODY") or sys.stdin.read()
    try:
        r = clean(extract(body))
    except Bad as e:
        print(f"Not added: {e}.")
        return 1
    reports = json.loads(FILE.read_text(encoding="utf-8")) if FILE.exists() else []
    if any(o["spot"] == r["spot"] and o["t"] == r["t"] and o["obs"] == r["obs"] for o in reports):
        print(f"Already in the shared reports ({r['spot']} {r['t']}, {r['obs']} ft).")
        return 0
    reports.append(r)
    reports.sort(key=lambda o: (o["t"], o["spot"]))
    # One report per line, like the hand-written file, so diffs stay readable.
    FILE.write_text("[\n" + ",\n".join(json.dumps(o, separators=(",", ":"), ensure_ascii=False) for o in reports)
                    + "\n]\n", encoding="utf-8")
    print(f"Added {r['spot']} {r['t']}, {r['obs']} ft. The site redeploys in a few minutes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
