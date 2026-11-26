#!/usr/bin/env python3
"""P1A.02: read-only source preflight. No files are downloaded to disk."""
import concurrent.futures
import json
import re
import subprocess
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit
import csv

root = Path(__file__).resolve().parents[2]
manifest = root / "manifests/P1A02_SOURCE_CANDIDATES.tsv"

def curl(url, head=False, limit=5000000):
    args = ["curl", "-sS", "-L", "--connect-timeout", "10",
            "--max-time", "35", "--retry", "0"]
    if head:
        args += ["-I", "-o", "/dev/null", "-w",
                 "%{http_code}\t%{content_type}\t%{url_effective}"]
    else:
        args += ["--fail", "--max-filesize", str(limit)]
    return subprocess.run(args + [url], capture_output=True, text=True)

def probe(row):
    p = curl(row["url"], head=True)
    result = p.stdout.strip()
    return (row["id"], row["role"], result or "NO_HTTP_RESULT",
            f"curl_exit={p.returncode}", p.stderr.strip()[:200])

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = set()
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            href = dict(attrs).get("href", "")
            filename = urlsplit(href).path.rsplit("/", 1)[-1]
            if filename.endswith(".json"):
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*\.json", filename):
                    raise ValueError(f"Unexpected JSON filename {filename!r}")
                self.links.add(filename)

print("CFA_P1A02_SOURCE_PREFLIGHT_V1_1", flush=True)
print("No bulk transfers; HEAD results are accessibility checks, not integrity certification.")
version_url = "https://reactome.org/ContentService/data/database/version"
p = curl(version_url, limit=10000)
print("REACTOME_VERSION", repr(p.stdout.strip()), "curl_exit", p.returncode)
if p.stderr:
    print("VERSION_ERROR", p.stderr.strip()[:300])
print("EXPECTED_REACTOME_VERSION 98")

with manifest.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
print("\nID\tROLE\tHTTP_STATUS\tCONTENT_TYPE\tEFFECTIVE_URL\tCURL_EXIT\tERROR")
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for result in pool.map(probe, rows):
        print("\t".join(result), flush=True)

base = "https://release.geneontology.org/2026-08-05/go-cams/json/"
p = curl(base + "index.html")
parser = Links()
if p.returncode == 0:
    parser.feed(p.stdout)
print("\nGO_CAM_INDEX", "curl_exit", p.returncode, "MODEL_LINKS", len(parser.links))
if p.stderr:
    print("INDEX_ERROR", p.stderr.strip()[:300])
if parser.links:
    names = sorted(parser.links)
    for name in sorted({names[0], names[len(names)//2], names[-1]}):
        p = curl(base + name)
        try:
            obj = json.loads(p.stdout) if p.returncode == 0 else None
            if not isinstance(obj, dict):
                raise ValueError("expected JSON object")
            print("GO_CAM_PILOT", name, "KEYS", sorted(obj))
            for k,v in obj.items():
                print(" ", k, type(v).__name__,
                      len(v) if isinstance(v,(list,dict)) else str(v)[:100])
        except (ValueError, TypeError) as e:
            print("GO_CAM_PILOT_ERROR", name, p.returncode, str(e)[:200], p.stderr[:200])
print("\nPREFLIGHT_COMPLETE; SOURCE_CERTIFICATION=PENDING")
