#!/usr/bin/env python3
"""Create and publish the QESPM Zenodo deposit (production zenodo.org).

Usage:
    ZENODO_TOKEN=your_token python3 zenodo_deposit.py [archive.zip]

The script reads the token ONLY from the environment (never from the command
line).  It creates a software deposition, uploads the repository archive
produced by `git archive --format=zip -o /tmp/qespm_v1.0.0.zip HEAD`,
publishes it, and prints the final DOI.

After running, paste the printed DOI into the manuscript's Data Availability
section and the cover letter.
"""
import json
import os
import sys
import urllib.request
import urllib.parse

TOKEN = os.environ.get("ZENODO_TOKEN") or os.environ.get("ZENODO_ACCESS_TOKEN")
if not TOKEN:
    import getpass
    TOKEN = getpass.getpass("Zenodo token (input hidden, used once, never saved): ").strip()
if not TOKEN:
    sys.exit("No token provided; cannot create a Zenodo deposit.")

ARCHIVE = sys.argv[1] if len(sys.argv) > 1 else "/tmp/qespm_v1.0.0.zip"
BASE = "https://zenodo.org/api"

metadata = {
    "metadata": {
        "title": ("QESPM: Quantum Electrostatic Scanning Probe Microscopy "
                  "--- code, manuscript and synthetic benchmark dataset"),
        "upload_type": "software",
        "publication_type": "other",
        "description": (
            "Reference implementation, LaTeX sources and synthetic benchmark "
            "dataset for the manuscript 'Quantum Limits of Surface Metrology: "
            "Spatial Information Extraction by a Single Trapped Ion'. "
            "Contents: the forward/inverse modelling chain, figure-generation "
            "scripts, the unit-and-factor audit, and benchmark/qespm_benchmark_"
            "data.npz (reference surfaces, simulated ion-signal maps and "
            "reference reconstructions)."),
        "creators": [
            {"name": "Kucharski, Dawid",
             "affiliation": "Poznan University of Technology"},
        ],
        "license": "mit",
        "access_right": "open",
        "version": "1.0.0",
        "keywords": [
            "surface metrology", "trapped ion", "quantum sensing",
            "instrument transfer function", "inverse problems",
            "measurement uncertainty",
        ],
        "related_identifiers": [
            {"relation": "isSupplementTo",
             "identifier": "https://github.com/dawidkucharski/qespm",
             "resource_type": "software"},
        ],
    }
}


def request(method, url, token=None, data=None, headers=None, raw=None):
    h = {"Authorization": f"Bearer {TOKEN}"}
    if headers:
        h.update(headers)
    body = None
    if data is not None:
        body = json.dumps(data).encode()
        h["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, headers=h, method=method)
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode()) if r.headers.get_content_type() == "application/json" else r.read()


dep = request("POST", f"{BASE}/depositions", data=metadata)
dep_id = dep["id"]
print(f"deposition {dep_id} created (draft), DOI reserved: "
      f"{dep['metadata']['prereserve_doi']['doi']}")

# upload archive
url = f"{dep['links']['bucket']}/{os.path.basename(ARCHIVE)}"
with open(ARCHIVE, "rb") as f:
    raw = f.read()
request("PUT", url, headers={"Content-Type": "application/zip"}, raw=raw)
print("archive uploaded")

# publish
published = request("POST", f"{BASE}/depositions/{dep_id}/actions/publish")
doi = published["metadata"]["prereserve_doi"]["doi"]
print(f"PUBLISHED. Final DOI: {doi}")
print("Insert this DOI into manuscript/main.tex (Data Availability) and the cover letter.")
