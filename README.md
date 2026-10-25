# C-CDA to FHIR conversion fidelity with paired synthetic records

Code for the preprint *How Much Clinical Content Survives C-CDA to FHIR Conversion? Evaluating an Open-Source Converter with Paired Synthetic Records*.

Synthea exports each synthetic patient both as a C-CDA document and as a native FHIR R4 bundle. This code converts the C-CDA with python-fhir-converter and compares three layers of the same record, event by event:

- **N** native FHIR (reference standard)
- **S** C-CDA source
- **C** converted FHIR

No real patient data is used.

## Requirements

- Java 17+ (for Synthea)
- Python 3.13 (the macOS system `python3` is 3.9 and will not work; use `python3.13`)

```bash
python3.13 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install --no-deps -r requirements-nodeps.txt
```

`python-fhir-converter==0.3.0` is kept out of `requirements.txt` because it hard-pins old versions of its dependencies (e.g. `pyjson5==1.6.5`, `lxml==5.3.0`) that do not build on Python 3.13. Installing it with `--no-deps` lets it run on the newer versions in `requirements.txt`; `pip check` will report the version mismatch, which is expected.

## Reproduce

```bash
# 1. Synthetic data (master-branch-latest synthea-with-dependencies.jar from the Synthea GitHub releases)
java -jar synthea-with-dependencies.jar -s 20261024 -cs 20261024 -p 300 \
  --exporter.baseDirectory=data/output \
  --exporter.ccda.export=true --exporter.fhir.export=true \
  --exporter.hospital.fhir.export=false --exporter.practitioner.fhir.export=false \
  --exporter.csv.export=false --exporter.years_of_history=10 Massachusetts

# 2. Convert every C-CDA document
python -I src/convert.py data/output/ccda results/converted results/convert_log.csv

# 3. Measure fidelity
python -I src/analyze.py .
```

## Files

| File | Purpose |
| --- | --- |
| `src/convert.py` | Converts each C-CDA with the CCD template and logs time and errors |
| `src/extract.py` | Extracts comparable events from native FHIR, C-CDA and converted FHIR |
| `src/analyze.py` | Matching, recall and precision, field fidelity, validation, US Core checks |
| `results/family_metrics.csv` | Recall and precision by clinical family |
| `results/validation.csv` | Structural validation failures by resource type |
| `results/obs_export_loss.csv` | Native observations absent from the C-CDA, by category |
| `results/per_patient.csv` | Per-patient recall and export coverage |
| `results/summary.json` | All headline figures, field fidelity and recoding counts |

## License

MIT
