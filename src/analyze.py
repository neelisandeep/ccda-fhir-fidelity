"""Measure C-CDA -> FHIR conversion fidelity against Synthea's native FHIR export.

Usage: python -I src/analyze.py <study_root>
Reads   data/output/{fhir,ccda}, results/converted, results/convert_log.csv
Writes  results/summary.json, results/family_metrics.csv, results/validation.csv,
        results/obs_export_loss.csv, results/per_patient.csv
"""
from __future__ import annotations

import csv
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract import ccda_events, fhir_events, has_time, is_canonical  # noqa: E402

from fhir.resources.R4B import construct_fhir_element  # noqa: E402

FAMILIES = ["Condition", "Procedure", "Encounter", "Immunization", "Medication", "Allergy",
            "Observation", "DiagnosticReport"]
CLINICAL_TYPES = {"Condition", "Procedure", "Encounter", "Immunization", "MedicationRequest",
                  "MedicationStatement", "Medication", "AllergyIntolerance", "Observation",
                  "DiagnosticReport", "Patient"}
ACTCODE = "http://terminology.hl7.org/CodeSystem/v3-ActCode"


def code_of(e: dict, lenient: bool):
    """Strict: the code in its proper element. Lenient: also accept a code the converter moved
    to another element (Encounter.class, AllergyIntolerance.reaction.substance) and a source
    translation code that replaced the source's primary code."""
    if not lenient:
        return e["code"]
    if e.get("code_l"):
        return e["code_l"]
    if e["code"] is None and e.get("relocated"):
        return e["relocated"].get("code")
    return e["code"]


def annotate_recoding(src: list[dict], cv: list[dict]) -> None:
    """Mark converted events whose code is a source *translation* rather than the source's own
    code, and link each to the source event it came from (same family, day and value).
    'ambiguous' = two or more different source codes share that translation on that day,
    so the converted code can no longer tell them apart (e.g. left vs right eye)."""
    tsrc = defaultdict(list)
    for e in src:
        for t in e.get("alt") or []:
            tsrc[(e["family"], t, e["day"])].append(e)
    primary = Counter((e["family"], e["code"], e["day"]) for e in src)
    used = set()
    for e in cv:
        c = e["code"] if e["code"] is not None else (e.get("relocated") or {}).get("code")
        k = (e["family"], c, e["day"])
        if k not in tsrc or primary[k]:
            continue
        cands = [s for s in tsrc[k] if id(s) not in used]
        if not cands:
            continue
        s = next((s for s in cands if same_value(s["value"], e["value"])), cands[0])
        used.add(id(s))
        e["code_l"] = s["code"]
        e["recoded"] = True
        e["ambiguous"] = len({x["code"] for x in tsrc[k]}) > 1


def same_value(a, b) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-9)
    return str(a) == str(b)


def match(ref: list[dict], cand: list[dict], lenient: bool):
    """Multiset match on (family, code, day). Inside a key, equal values pair first."""
    gr, gc = defaultdict(list), defaultdict(list)
    for e in ref:
        gr[(e["family"], code_of(e, lenient), e["day"])].append(e)
    for e in cand:
        gc[(e["family"], code_of(e, lenient), e["day"])].append(e)
    pairs, un_ref, un_cand = [], [], []
    for k in set(gr) | set(gc):
        xr, xc = list(gr.get(k, [])), list(gc.get(k, []))
        if k[1] is None:  # an event without a code can never match
            un_ref += xr
            un_cand += xc
            continue
        rest_r = []
        for a in xr:
            j = next((i for i, b in enumerate(xc) if same_value(a["value"], b["value"])), None)
            if j is None:
                rest_r.append(a)
            else:
                pairs.append((a, xc.pop(j)))
        n = min(len(rest_r), len(xc))
        pairs += list(zip(rest_r[:n], xc[:n]))
        un_ref += rest_r[n:]
        un_cand += xc[n:]
    return pairs, un_ref, un_cand


def validate(bundle: dict) -> list[tuple[str, str]]:
    """(resourceType, first error location) for every resource failing R4B structure checks."""
    bad = []
    for ent in bundle.get("entry", []):
        r = ent["resource"]
        try:
            construct_fhir_element(r["resourceType"], r)
        except Exception as ex:  # pydantic ValidationError
            lines = [ln.strip() for ln in str(ex).splitlines()[1:3]]
            bad.append((r["resourceType"], " ".join(lines)))
    return bad


def us_core_presence(bundle: dict) -> Counter:
    """Counts of resources and of resources carrying selected US Core-required elements."""
    c = Counter()
    for ent in bundle.get("entry", []):
        r = ent["resource"]
        rt = r["resourceType"]
        if rt == "AllergyIntolerance":
            c["allergy_n"] += 1
            c["allergy_code"] += bool(r.get("code"))
        elif rt == "Condition":
            c["condition_n"] += 1
            c["condition_category"] += bool(r.get("category"))
            c["condition_clinical_status"] += bool(r.get("clinicalStatus"))
        elif rt == "Encounter":
            c["encounter_n"] += 1
            c["encounter_type"] += bool(r.get("type"))
            c["encounter_class_actcode"] += (r.get("class") or {}).get("system") == ACTCODE
        elif rt == "Observation":
            c["observation_n"] += 1
            c["observation_category"] += bool(r.get("category"))
        elif rt == "DiagnosticReport":
            c["report_n"] += 1
            c["report_category"] += bool(r.get("category"))
        elif rt in ("MedicationRequest", "MedicationStatement"):
            c["med_n"] += 1
            c["med_is_request"] += rt == "MedicationRequest"
            c["med_status"] += bool(r.get("status"))
    return c


def main(root: Path) -> None:
    gt_dir, src_dir = root / "data/output/fhir", root / "data/output/ccda"
    cv_dir, out = root / "results/converted", root / "results"
    names = sorted(p.stem for p in cv_dir.glob("*.json") if (gt_dir / (p.stem + ".json")).exists())

    n = defaultdict(Counter)            # family -> layer counts
    m = defaultdict(Counter)            # family -> matched counts per comparison
    fid = defaultdict(Counter)          # field fidelity counters per family
    obs_loss = Counter()                # GT observation category -> not exported
    obs_total = Counter()
    dr_gt_kind = Counter()
    resource_gt, resource_cv = Counter(), Counter()
    val = {"GT": Counter(), "CV": Counter()}
    val_n = {"GT": Counter(), "CV": Counter()}
    val_err = Counter()
    usc = {"GT": Counter(), "CV": Counter()}
    per_patient = []
    allergen_lost = Counter()
    recode_pairs = Counter()
    src_obs_by_section, src_obs_lost_by_section = Counter(), Counter()

    for name in names:
        gt_b = json.loads((gt_dir / f"{name}.json").read_text(encoding="utf-8"))
        cv_b = json.loads((cv_dir / f"{name}.json").read_text(encoding="utf-8"))
        gt, cv = fhir_events(gt_b), fhir_events(cv_b)
        src = ccda_events(src_dir / f"{name}.xml")
        annotate_recoding(src, cv)
        for e in cv:
            if e["family"] == "Allergy" and e["code"] is None and not (e.get("relocated") or {}).get("code"):
                allergen_lost["n"] += 1
        resource_gt.update(e["resource"]["resourceType"] for e in gt_b["entry"])
        resource_cv.update(e["resource"]["resourceType"] for e in cv_b["entry"])
        for lay, b in (("GT", gt_b), ("CV", cv_b)):
            val_n[lay].update(e["resource"]["resourceType"] for e in b["entry"])
            for rt, loc in validate(b):
                val[lay][rt] += 1
                if lay == "CV":
                    val_err[(rt, loc)] += 1
            usc[lay].update(us_core_presence(b))

        for lay, ev in (("gt", gt), ("src", src), ("cv", cv)):
            n_f = Counter(e["family"] for e in ev)
            for f in FAMILIES:
                n[f][lay] += n_f[f]

        p_exp, un_gt_exp, _ = match(gt, src, lenient=False)
        p_cs, _, un_cv_s = match(src, cv, lenient=False)
        p_cl, un_src_l, un_cv_l = match(src, cv, lenient=True)
        for e in src:
            if e["family"] == "Observation":
                src_obs_by_section[e.get("section")] += 1
        for e in un_src_l:
            if e["family"] == "Observation":
                src_obs_lost_by_section[e.get("section")] += 1
        p_e2e, _, _ = match(gt, cv, lenient=True)
        for name_, pairs in (("export", p_exp), ("conv_strict", p_cs), ("conv_lenient", p_cl),
                             ("e2e", p_e2e)):
            for a, _b in pairs:
                m[a["family"]][name_] += 1
        for e in un_cv_l:
            m[e["family"]]["cv_unmatched"] += 1

        for e in gt:
            if e["family"] == "Observation":
                obs_total[e.get("category") or "none"] += 1
        for e in un_gt_exp:
            if e["family"] == "Observation":
                obs_loss[e.get("category") or "none"] += 1
        for r in (x["resource"] for x in gt_b["entry"]):
            if r["resourceType"] == "DiagnosticReport":
                dr_gt_kind["note" if r.get("presentedForm") else "results"] += 1

        # field fidelity: conversion step (source vs converted)
        for s, c in p_cl:
            f = s["family"]
            fid[f]["pairs"] += 1
            if has_time(s["ts"]):
                fid[f]["src_has_time"] += 1
                fid[f]["time_kept"] += has_time(c["ts"])
            sys_raw = c["system_raw"] if c["code"] is not None else (c.get("relocated") or {}).get("system")
            fid[f]["canonical_system"] += is_canonical(sys_raw)
            fid[f]["relocated"] += c["code"] is None
            fid[f]["recoded"] += bool(c.get("recoded"))
            fid[f]["recoded_ambiguous"] += bool(c.get("ambiguous"))
            if c.get("recoded"):
                recode_pairs[(f, s["code"], c["code"] if c["code"] is not None
                              else (c.get("relocated") or {}).get("code"))] += 1
            if s["value"] is not None:
                fid[f]["src_value"] += 1
                fid[f]["value_kept"] += same_value(s["value"], c["value"])
                if s["unit"]:
                    fid[f]["src_unit"] += 1
                    fid[f]["unit_kept"] += (s["unit"] == c["unit"])
            if f == "Condition" and s.get("abatement"):
                fid[f]["src_abatement"] += 1
                fid[f]["abatement_kept"] += (s["abatement"] == c.get("abatement"))
        # status and encounter links: end to end (native vs converted)
        for g, c in p_e2e:
            f = g["family"]
            fid[f]["e2e_pairs"] += 1
            if g.get("status") is not None:
                fid[f]["gt_status"] += 1
                fid[f]["status_present"] += c.get("status") is not None
                fid[f]["status_equal"] += c.get("status") == g.get("status")
            if f != "Encounter":
                fid[f]["gt_enc_link"] += bool(g.get("has_enc"))
                fid[f]["cv_enc_link"] += bool(c.get("has_enc"))

        tot_src = len([e for e in src if e["family"] in FAMILIES])
        tot_conv = sum(1 for a, _ in p_cl)
        tot_gt = len(gt)
        per_patient.append({"patient": name, "gt_events": tot_gt, "src_events": tot_src,
                            "conv_matched": tot_conv,
                            "conv_recall": round(tot_conv / tot_src, 4) if tot_src else None,
                            "export_coverage": round(len(p_exp) / tot_gt, 4) if tot_gt else None})

    # ---- write outputs
    rows = []
    for f in FAMILIES:
        a, b = n[f], m[f]
        rows.append({
            "family": f, "native_fhir": a["gt"], "ccda": a["src"], "converted": a["cv"],
            "export_coverage": b["export"] / a["gt"] if a["gt"] else None,
            "conversion_recall_strict": b["conv_strict"] / a["src"] if a["src"] else None,
            "conversion_recall_lenient": b["conv_lenient"] / a["src"] if a["src"] else None,
            "conversion_precision_lenient": b["conv_lenient"] / a["cv"] if a["cv"] else None,
            "end_to_end_recall": b["e2e"] / a["gt"] if a["gt"] else None,
            "converted_unmatched": b["cv_unmatched"],
        })
    with (out / "family_metrics.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    with (out / "validation.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["layer", "resource_type", "resources", "invalid"])
        for lay in ("GT", "CV"):
            for rt in sorted(val_n[lay]):
                w.writerow([lay, rt, val_n[lay][rt], val[lay][rt]])
    with (out / "obs_export_loss.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["category", "native_observations", "not_in_ccda"])
        for cat, tot in obs_total.most_common():
            w.writerow([cat, tot, obs_loss[cat]])
    with (out / "per_patient.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=per_patient[0].keys())
        w.writeheader()
        w.writerows(per_patient)

    log = list(csv.DictReader((out / "convert_log.csv").open()))
    secs = sorted(float(r["seconds"]) for r in log)
    sizes = sorted(int(r["bytes"]) for r in log)
    conv_r = sorted(p["conv_recall"] for p in per_patient if p["conv_recall"] is not None)
    exp_r = sorted(p["export_coverage"] for p in per_patient if p["export_coverage"] is not None)

    def q(xs, p):
        return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))]

    summary = {
        "patients": len(names),
        "ccda_files": len(log),
        "conversion_errors": sum(r["status"] != "ok" for r in log),
        "seconds_total": round(sum(secs), 1),
        "seconds_median": q(secs, 0.5), "seconds_p95": q(secs, 0.95), "seconds_max": secs[-1],
        "ccda_kb_median": round(q(sizes, 0.5) / 1024, 1), "ccda_kb_max": round(sizes[-1] / 1024, 1),
        "events_native": sum(n[f]["gt"] for f in FAMILIES),
        "events_ccda": sum(n[f]["src"] for f in FAMILIES),
        "events_converted": sum(n[f]["cv"] for f in FAMILIES),
        "matched_export": sum(m[f]["export"] for f in FAMILIES),
        "matched_conv_strict": sum(m[f]["conv_strict"] for f in FAMILIES),
        "matched_conv_lenient": sum(m[f]["conv_lenient"] for f in FAMILIES),
        "matched_e2e": sum(m[f]["e2e"] for f in FAMILIES),
        "per_patient_conv_recall": {"min": conv_r[0], "median": q(conv_r, 0.5), "max": conv_r[-1]},
        "per_patient_export_coverage": {"min": exp_r[0], "median": q(exp_r, 0.5), "max": exp_r[-1]},
        "fidelity": {f: dict(fid[f]) for f in FAMILIES},
        "us_core_presence": {k: dict(v) for k, v in usc.items()},
        "validation_top_errors": [[rt, loc, cnt] for (rt, loc), cnt in val_err.most_common(15)],
        "resources_native": dict(resource_gt.most_common()),
        "resources_converted": dict(resource_cv.most_common()),
        "diagnostic_report_native_kind": dict(dr_gt_kind),
        "allergen_lost_entirely": allergen_lost["n"],
        "recodings": [[f, a, b, k] for (f, a, b), k in recode_pairs.most_common()],
        "ccda_observations_by_section": dict(src_obs_by_section),
        "ccda_observations_lost_by_section": dict(src_obs_lost_by_section),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, dict)}, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
