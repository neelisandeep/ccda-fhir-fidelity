"""Convert every Synthea C-CDA document to FHIR R4 with python-fhir-converter (CCD template).
Usage: python -I convert.py <ccda_dir> <out_dir> <log_csv>"""
import csv, sys, time, traceback
from pathlib import Path
from fhir_converter.renderers import CcdaRenderer

src, out, log = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
renderer = CcdaRenderer()
rows = []
for f in sorted(src.glob("*.xml")):
    t0 = time.perf_counter()
    status, err = "ok", ""
    try:
        with f.open(encoding="utf-8") as fh:
            res = renderer.render_fhir_string("CCD", fh)
        (out / (f.stem + ".json")).write_text(res, encoding="utf-8")
    except Exception as e:  # record, do not stop
        status, err = "error", f"{type(e).__name__}: {str(e)[:200]}"
    rows.append({"file": f.name, "bytes": f.stat().st_size, "seconds": round(time.perf_counter() - t0, 3),
                 "status": status, "error": err})
with log.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
ok = sum(r["status"] == "ok" for r in rows)
print(f"converted {ok}/{len(rows)}")
