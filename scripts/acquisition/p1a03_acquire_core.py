#!/usr/bin/env python3
"""CFA P1A.03 core acquisition. Standard library + curl; no extraction."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from urllib.parse import urlsplit
import zipfile

ROOT = Path(__file__).resolve().parents[2]
CANDIDATES = ROOT / "manifests/P1A02_SOURCE_CANDIDATES.tsv"
RAW = ROOT / "data/raw"
VERSION_URL = "https://reactome.org/ContentService/data/database/version"
EXPECTED = "98"

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def curl_base():
    return ["curl", "--fail", "--silent", "--show-error", "--location",
            "--proto", "=https", "--proto-redir", "=https",
            "--retry", "3", "--connect-timeout", "30"]

def version():
    result = subprocess.check_output(curl_base() + [
        "--max-time", "120", VERSION_URL], text=True).strip()
    if result != EXPECTED:
        raise RuntimeError(f"Reactome version {result!r}; expected {EXPECTED}. STOP.")
    return result

def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(8*1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def check_payload(p, name):
    if p.stat().st_size == 0:
        raise ValueError("empty payload")
    with p.open("rb") as f:
        prefix = f.read(8192)
    html = prefix.lstrip().lower().startswith((b"<!doctype html", b"<html"))
    if html and not name.endswith(".html"):
        raise ValueError("unexpected HTML response")
    if name.endswith(".zip"):
        with zipfile.ZipFile(p) as z:
            if not z.namelist():
                raise ValueError("empty ZIP")
            bad = z.testzip()
            if bad is not None:
                raise ValueError(f"ZIP CRC failure: {bad}")
        return "ZIP_CRC_PASS"
    if name.endswith((".gz", ".tgz")):
        with gzip.open(p, "rb") as f:
            if not f.read(1):
                raise ValueError("empty gzip payload")
            while f.read(8*1024*1024):
                pass
        return "GZIP_CRC_PASS"
    if name.endswith(".json"):
        with p.open(encoding="utf-8") as f:
            json.load(f)
        return "JSON_PARSE_PASS"
    if name.endswith(".dump"):
        return "NONEMPTY_ONLY_NATIVE_FORMAT_UNVERIFIED"
    if name.endswith(".html"):
        if not html and b"<html" not in prefix.lower():
            raise ValueError("index does not resemble HTML")
        return "HTML_PREFIX_PASS"
    if b"\0" in prefix:
        raise ValueError("unexpected binary prefix in text source")
    prefix.decode("utf-8", errors="strict")
    return "TEXT_PREFIX_PASS"

def destination(row, raw=RAW):
    version_dir = {"Reactome": "reactome_v98",
                   "GO": "go_2026-08-05",
                   "HGNC": "hgnc_snapshot"}[row["resource"]]
    name = Path(urlsplit(row["url"]).path).name
    return raw / version_dir / row["role"] / name

def atomic_json(path, obj):
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("x", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)

def acquire(row):
    p = destination(row)
    meta = p.with_name(p.name + ".receipt.json")
    tmp = p.with_name(p.name + ".partial")
    meta_tmp = meta.with_name(meta.name + ".partial")
    reactome = row["resource"] == "Reactome"
    if p.exists() or meta.exists() or tmp.exists() or meta_tmp.exists():
        if not p.is_file() or not meta.is_file() or tmp.exists() or meta_tmp.exists():
            raise RuntimeError(f"Unpaired/partial artifact: {p}; inspect before retry")
        m = json.loads(meta.read_text())
        for k in ("id", "role", "resource", "release_candidate", "url"):
            if m.get(k) != row[k]:
                raise RuntimeError(f"Receipt identity mismatch {k}: {p}")
        if m["bytes"] != p.stat().st_size or m["sha256"] != sha256(p):
            raise RuntimeError(f"Existing file integrity mismatch: {p}")
        if reactome and (m.get("reactome_version_before") != EXPECTED
                         or m.get("reactome_version_after") != EXPECTED):
            raise RuntimeError(f"Existing Reactome version receipt mismatch: {p}")
        print("REUSE", row["id"], m["bytes"], m["sha256"], flush=True)
        return m
    p.parent.mkdir(parents=True, exist_ok=True)
    before = version() if reactome else None
    started = utc()
    print("GET", row["id"], row["url"], flush=True)
    try:
        info = subprocess.check_output(curl_base() + [
            "--max-time", "3600", "--output", str(tmp), "--write-out",
            "%{http_code}\t%{url_effective}\t%{content_type}",
            row["url"]], text=True).split("\t", 2)
        if len(info) != 3 or info[0] != "200":
            raise RuntimeError(f"Unexpected HTTP response: {info}")
        after = version() if reactome else None
        validation = check_payload(tmp, p.name)
        m = {**row, "path": str(p.relative_to(ROOT)),
             "bytes": tmp.stat().st_size, "sha256": sha256(tmp),
             "http_status": info[0], "effective_url": info[1],
             "content_type": info[2], "started_utc": started,
             "acquired_utc": utc(), "reactome_version_before": before,
             "reactome_version_after": after, "payload_check": validation,
             "license_status": "PENDING_REVIEW",
             "status": "ACQUIRED_NOT_CERTIFIED"}
        # Fail closed: a crash between file and receipt leaves an unpaired file.
        os.replace(tmp, p)
        atomic_json(meta, m)
        print("ACQUIRED", row["id"], m["bytes"], m["sha256"], validation, flush=True)
        return m
    except Exception:
        print("STOP", row["id"], "retain partial artifacts for inspection", flush=True)
        raise

def index_summary(path):
    obj = json.loads(path.read_text())
    print("TAXON_INDEX_ROOT", type(obj).__name__, "length", len(obj) if isinstance(obj,(list,dict)) else "NA")
    if isinstance(obj, dict):
        print("TAXON_INDEX_KEYS", list(obj)[:30])
        for key, val in list(obj.items())[:4]:
            print("TAXON_INDEX_ENTRY", key, type(val).__name__, json.dumps(val)[:1500])
        for key, val in obj.items():
            if str(key) in ("NCBITaxon:9606","9606","Homo sapiens"):
                print("HUMAN_KEY", key, type(val).__name__, "length",
                      len(val) if isinstance(val,(dict,list)) else "NA", json.dumps(val)[:2500])
    elif isinstance(obj,list):
        for val in obj[:3]:
            print("TAXON_INDEX_ENTRY", json.dumps(val)[:1500])

def self_test():
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"payload"
        p.write_text("test\n")
        assert check_payload(p,"test.txt")=="TEXT_PREFIX_PASS"
        assert sha256(p)==hashlib.sha256(b"test\n").hexdigest()
        p.write_text("<html>error</html>")
        try:
            check_payload(p,"test.txt")
        except ValueError:
            pass
        else:
            raise AssertionError("HTML rejection failed")
        with zipfile.ZipFile(p,"w") as z:
            z.writestr("fixture.txt","test")
        assert check_payload(p,"test.zip")=="ZIP_CRC_PASS"
        with gzip.open(p,"wb") as f:
            f.write(b"test")
        assert check_payload(p,"test.gz")=="GZIP_CRC_PASS"
        p.write_text('{"test": true}')
        assert check_payload(p,"test.json")=="JSON_PARSE_PASS"
        atomic_json(Path(d)/"receipt.json", {"test":True})
        assert json.loads((Path(d)/"receipt.json").read_text())=={"test":True}
    print("SELF_TEST_PASS (local payload/checksum fixtures only)")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--execute",action="store_true")
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return
    with CANDIDATES.open() as f:
        rows=[r for r in csv.DictReader(f,delimiter="\t") if r["id"]!="R14"]
    if len(rows)!=19 or len({r["id"] for r in rows})!=19:
        raise RuntimeError("unexpected source manifest; review")
    print("CFA_P1A03_CORE_ACQUISITION_V1", flush=True)
    print("MODE", "EXECUTE" if args.execute else "PLAN_ONLY", flush=True)
    print("FILES", len(rows), "FILESYSTEM_FREE_BYTES", shutil.disk_usage(ROOT).free)
    for r in rows:
        print("PLAN",r["id"],r["role"],str(destination(r).relative_to(ROOT)))
    if not args.execute:
        print("No changes. Confirm storage allocation and site download policy before --execute.")
        return
    version()
    receipts=[]
    for row in rows:
        receipts.append(acquire(row))
    version()
    # Local snapshot, ignored by Git. Immutable per-file receipts are authoritative.
    snapshot=RAW/"P1A03_CORE_RECEIPTS.json"
    if snapshot.exists():
        if json.loads(snapshot.read_text()) != receipts:
            raise RuntimeError("Existing aggregate receipt mismatch")
    else:
        atomic_json(snapshot,receipts)
    index_summary(destination(next(r for r in rows if r["id"]=="G02")))
    print("CORE_ACQUIRED",len(receipts),"TOTAL_BYTES",sum(m["bytes"] for m in receipts))
    print("REACTOME_VERSION_FINAL",EXPECTED)
    print("GO_CAM_INDIVIDUAL_MODELS_PENDING")
    print("SOURCE_CERTIFICATION_PENDING")

if __name__=="__main__":
    main()
