"""Extract comparable clinical events from three layers for one synthetic patient.

Layer GT : Synthea's native FHIR R4 bundle (reference standard)
Layer SRC: Synthea's C-CDA document for the same patient (what the converter receives)
Layer CV : FHIR R4 produced from the C-CDA by python-fhir-converter

Every event is a dict with: family, system, code, day, ts (raw timestamp), value, unit,
status, extra (family-specific fields). Families: Condition, Observation, Procedure,
Immunization, Medication, Allergy, Encounter, DiagnosticReport.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from lxml import etree

NS = {"h": "urn:hl7-org:v3", "xsi": "http://www.w3.org/2001/XMLSchema-instance"}

OID = {
    "2.16.840.1.113883.6.96": "SNOMED",
    "2.16.840.1.113883.6.1": "LOINC",
    "2.16.840.1.113883.6.88": "RXNORM",
    "2.16.840.1.113883.12.292": "CVX",
    "2.16.840.1.113883.6.12": "CPT",
    "2.16.840.1.113883.6.90": "ICD10CM",
    "2.16.840.1.113883.6.69": "NDC",
    "2.16.840.1.113883.6.238": "CDCREC",
    "2.16.840.1.113883.5.4": "ACTCODE",
}
CANONICAL = {
    "http://snomed.info/sct": "SNOMED",
    "http://loinc.org": "LOINC",
    "http://www.nlm.nih.gov/research/umls/rxnorm": "RXNORM",
    "http://hl7.org/fhir/sid/cvx": "CVX",
    "http://www.ama-assn.org/go/cpt": "CPT",
    "http://hl7.org/fhir/sid/icd-10-cm": "ICD10CM",
    "http://hl7.org/fhir/sid/ndc": "NDC",
    "http://terminology.hl7.org/CodeSystem/v3-ActCode": "ACTCODE",
}


def norm_system(system: str | None) -> str:
    """Map a FHIR system URI, urn:oid URI or bare OID to a short terminology name."""
    if not system:
        return "NONE"
    if system in CANONICAL:
        return CANONICAL[system]
    oid = system[8:] if system.startswith("urn:oid:") else system
    return OID.get(oid, system)


def is_canonical(system: str | None) -> bool:
    return bool(system) and system in CANONICAL


def day_of(ts: str | None) -> str | None:
    """YYYY-MM-DD from a FHIR dateTime or an HL7 v3 TS (YYYYMMDDhhmmss)."""
    if not ts:
        return None
    if re.fullmatch(r"\d{8}.*", ts):
        return f"{ts[0:4]}-{ts[4:6]}-{ts[6:8]}"
    return ts[:10]


def has_time(ts: str | None) -> bool:
    if not ts:
        return False
    return ("T" in ts) if "-" in ts else len(re.sub(r"\D", "", ts)) > 8


# ---------------------------------------------------------------- FHIR layers
def _first_coding(cc: dict | None) -> dict:
    if not cc:
        return {}
    codings = cc.get("coding") or []
    return codings[0] if codings else {}


def _ev(family, coding, ts, rid, value=None, unit=None, status=None, **extra):
    return {
        "family": family,
        "system_raw": coding.get("system"),
        "system": norm_system(coding.get("system")),
        "code": coding.get("code"),
        "ts": ts,
        "day": day_of(ts),
        "value": value,
        "unit": unit,
        "status": status,
        "rid": rid,
        **extra,
    }


def _obs_value(r: dict):
    if "valueQuantity" in r:
        q = r["valueQuantity"]
        return q.get("value"), q.get("unit") or q.get("code")
    if "valueCodeableConcept" in r:
        return _first_coding(r["valueCodeableConcept"]).get("code"), None
    for k in ("valueString", "valueInteger", "valueBoolean"):
        if k in r:
            return r[k], None
    return None, None


def fhir_events(bundle: dict) -> list[dict]:
    res = [e["resource"] for e in bundle.get("entry", [])]
    by_ref = {}
    for r in res:
        by_ref[f"{r['resourceType']}/{r.get('id')}"] = r
    for e in bundle.get("entry", []):
        if e.get("fullUrl"):
            by_ref[e["fullUrl"]] = e["resource"]

    out = []
    for r in res:
        rt, rid = r["resourceType"], r.get("id")
        enc = bool((r.get("encounter") or r.get("context") or {}).get("reference"))
        prof = bool((r.get("meta") or {}).get("profile"))
        if rt == "Condition":
            c = _first_coding(r.get("code"))
            status = _first_coding(r.get("clinicalStatus")).get("code")
            out.append(_ev("Condition", c, r.get("onsetDateTime") or r.get("recordedDate"), rid,
                           status=status, abatement=day_of(r.get("abatementDateTime")),
                           has_enc=enc, profile=prof))
        elif rt == "Observation":
            ts = r.get("effectiveDateTime") or (r.get("effectivePeriod") or {}).get("start")
            cat = _first_coding((r.get("category") or [{}])[0]).get("code")
            v, u = _obs_value(r)
            if v is not None:
                out.append(_ev("Observation", _first_coding(r.get("code")), ts, rid, v, u,
                               r.get("status"), category=cat, has_enc=enc, profile=prof))
            for comp in r.get("component") or []:
                cv, cu = _obs_value(comp)
                out.append(_ev("Observation", _first_coding(comp.get("code")), ts, rid, cv, cu,
                               r.get("status"), category=cat, has_enc=enc, profile=prof,
                               component=True))
        elif rt == "Procedure":
            ts = r.get("performedDateTime") or (r.get("performedPeriod") or {}).get("start")
            out.append(_ev("Procedure", _first_coding(r.get("code")), ts, rid,
                           status=r.get("status"), has_enc=enc, profile=prof))
        elif rt == "Immunization":
            out.append(_ev("Immunization", _first_coding(r.get("vaccineCode")),
                           r.get("occurrenceDateTime"), rid, status=r.get("status"),
                           has_enc=enc, profile=prof))
        elif rt in ("MedicationRequest", "MedicationStatement"):
            cc = r.get("medicationCodeableConcept")
            if cc is None:
                ref = (r.get("medicationReference") or {}).get("reference")
                cc = (by_ref.get(ref) or {}).get("code")
            if rt == "MedicationRequest":
                ts = r.get("authoredOn")
            else:
                ts = (r.get("effectivePeriod") or {}).get("start") or r.get("effectiveDateTime")
            out.append(_ev("Medication", _first_coding(cc), ts, rid, status=r.get("status"),
                           resource_type=rt, has_enc=enc, profile=prof))
        elif rt == "AllergyIntolerance":
            # The allergen belongs in AllergyIntolerance.code; record where it actually sits.
            reac = (r.get("reaction") or [{}])[0]
            out.append(_ev("Allergy", _first_coding(r.get("code")),
                           r.get("onsetDateTime") or r.get("recordedDate"), rid,
                           status=_first_coding(r.get("clinicalStatus")).get("code"),
                           relocated=_first_coding(reac.get("substance")),
                           has_enc=enc, profile=prof))
        elif rt == "Encounter":
            # The encounter type belongs in Encounter.type; record Encounter.class as well.
            t = _first_coding((r.get("type") or [{}])[0])
            cls = r.get("class") or {}
            out.append(_ev("Encounter", t, (r.get("period") or {}).get("start"), rid,
                           status=r.get("status"), class_system=cls.get("system"),
                           relocated={"system": cls.get("system"), "code": cls.get("code")},
                           profile=prof, has_enc=True))
        elif rt == "DiagnosticReport":
            ts = r.get("effectiveDateTime") or (r.get("effectivePeriod") or {}).get("start")
            out.append(_ev("DiagnosticReport", _first_coding(r.get("code")), ts, rid,
                           status=r.get("status"),
                           n_results=len(r.get("result") or []), has_enc=enc, profile=prof))
    return out


def load_fhir(path: Path) -> list[dict]:
    return fhir_events(json.loads(path.read_text(encoding="utf-8")))


# ---------------------------------------------------------------- C-CDA layer
SEC = {
    "allergies": "2.16.840.1.113883.10.20.22.2.6.1",
    "medications": "2.16.840.1.113883.10.20.22.2.1.1",
    "results": "2.16.840.1.113883.10.20.22.2.3.1",
    "problems": "2.16.840.1.113883.10.20.22.2.5.1",
    "procedures": "2.16.840.1.113883.10.20.22.2.7.1",
    "encounters": "2.16.840.1.113883.10.20.22.2.22",
    "vitals": "2.16.840.1.113883.10.20.22.2.4.1",
    "immunizations": "2.16.840.1.113883.10.20.22.2.2",
    "social": "2.16.840.1.113883.10.20.22.2.17",
    "functional": "2.16.840.1.113883.10.20.22.2.14",
    "plan": "2.16.840.1.113883.10.20.22.2.10",
}


def _x(el, path):
    return el.xpath(path, namespaces=NS)


def _ts(el) -> str | None:
    """effectiveTime value, else low/@value."""
    et = _x(el, "h:effectiveTime")
    if not et:
        return None
    et = et[0]
    if et.get("value"):
        return et.get("value")
    low = _x(et, "h:low/@value")
    return low[0] if low else None


def _high(el) -> str | None:
    v = _x(el, "h:effectiveTime/h:high/@value")
    return v[0] if v else None


def _code(el) -> dict:
    c = _x(el, "h:code")
    if not c:
        return {}
    return {"system": c[0].get("codeSystem"), "code": c[0].get("code"),
            "alt": [t.get("code") for t in _x(c[0], "h:translation") if t.get("code")]}


def _cd_value(el) -> dict:
    v = _x(el, "h:value")
    if not v:
        return {}
    return {"system": v[0].get("codeSystem"), "code": v[0].get("code")}


def _obs_val(el):
    v = _x(el, "h:value")
    if not v:
        return None, None
    v = v[0]
    xtype = v.get("{http://www.w3.org/2001/XMLSchema-instance}type")
    if xtype in ("PQ", "INT", "REAL"):
        try:
            return float(v.get("value")), v.get("unit")
        except (TypeError, ValueError):
            return v.get("value"), v.get("unit")
    if xtype in ("CD", "CE", "CO"):
        return v.get("code"), None
    return (v.text or v.get("value")), None


def _src(family, coding, ts, value=None, unit=None, **extra):
    return {
        "family": family, "system_raw": coding.get("system"),
        "system": norm_system(coding.get("system")), "code": coding.get("code"),
        "alt": coding.get("alt", []),
        "ts": ts, "day": day_of(ts), "value": value, "unit": unit, **extra,
    }


def ccda_events(path: Path) -> list[dict]:
    doc = etree.parse(str(path))
    out = []

    def section(key):
        s = _x(doc, f"//h:section[h:templateId/@root='{SEC[key]}']")
        return s[0] if s else None

    s = section("problems")
    if s is not None:
        for obs in _x(s, "h:entry/h:act/h:entryRelationship/h:observation"):
            out.append(_src("Condition", _cd_value(obs), _ts(obs), abatement=day_of(_high(obs))))
    s = section("procedures")
    if s is not None:
        for p in _x(s, "h:entry/h:procedure"):
            out.append(_src("Procedure", _code(p), _ts(p)))
    s = section("immunizations")
    if s is not None:
        for sa in _x(s, "h:entry/h:substanceAdministration"):
            c = _x(sa, "h:consumable/h:manufacturedProduct/h:manufacturedMaterial")
            out.append(_src("Immunization", _code(c[0]) if c else {}, _ts(sa)))
    s = section("medications")
    if s is not None:
        for sa in _x(s, "h:entry/h:substanceAdministration"):
            c = _x(sa, "h:consumable/h:manufacturedProduct/h:manufacturedMaterial")
            out.append(_src("Medication", _code(c[0]) if c else {}, _ts(sa),
                            mood=sa.get("moodCode")))
    s = section("allergies")
    if s is not None:
        for obs in _x(s, "h:entry/h:act/h:entryRelationship/h:observation"):
            pe = _x(obs, "h:participant/h:participantRole/h:playingEntity")
            out.append(_src("Allergy", _code(pe[0]) if pe else {}, _ts(obs) or _ts(obs.getparent().getparent())))
    s = section("encounters")
    if s is not None:
        for enc in _x(s, "h:entry/h:encounter"):
            out.append(_src("Encounter", _code(enc), _ts(enc)))
    for key in ("results", "vitals"):
        s = section(key)
        if s is None:
            continue
        for org in _x(s, "h:entry/h:organizer"):
            if key == "results":
                out.append(_src("DiagnosticReport", _code(org), _ts(org),
                                n_results=len(_x(org, "h:component/h:observation"))))
            for obs in _x(org, "h:component/h:observation"):
                v, u = _obs_val(obs)
                out.append(_src("Observation", _code(obs), _ts(obs) or _ts(org), v, u, section=key))
    for key in ("social", "functional", "plan"):
        s = section(key)
        if s is None:
            continue
        for obs in _x(s, "h:entry/h:observation"):
            v, u = _obs_val(obs)
            out.append(_src("Observation", _code(obs), _ts(obs), v, u, section=key))
    return out
