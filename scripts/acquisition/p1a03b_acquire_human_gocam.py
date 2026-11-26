#!/usr/bin/env python3
"""P1A.03B: acquire human-indexed GO-CAMs, retaining all contexts and evidence."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
from urllib.parse import urlsplit

from p1a03_acquire_core import ROOT, RAW, sha256, utc, curl_base, atomic_json

BASE = "https://release.geneontology.org/2026-08-05/go-cams/json/"
HUMAN = "NCBITaxon:9606"
DIRECTORY = RAW / "go_2026-08-05" / "human_indexed_models"
ID = re.compile(r"[0-9a-f]+")
EXPECTED = ROOT / "manifests/P1A03_CORE_REPORTED_SHA256.tsv"

class Index(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids=set()
    def handle_starttag(self, tag, attrs):
        if tag!="a":
            return
        name=Path(urlsplit(dict(attrs).get("href","")).path).name
        if name.endswith(".json") and ID.fullmatch(name[:-5]):
            self.ids.add(name[:-5])

def load_indices():
    with EXPECTED.open() as f:
        rows={r["id"]:r for r in csv.DictReader(f,delimiter="\t")}
    paths={"G01": RAW/"go_2026-08-05/primary_index/index.html",
           "G02": RAW/"go_2026-08-05/scope_audit/taxon_index.json"}
    for key,p in paths.items():
        if p.stat().st_size!=int(rows[key]["bytes"]) or sha256(p)!=rows[key]["sha256"]:
            raise RuntimeError(f"Index snapshot mismatch {key}")
    idx=Index()
    idx.feed(paths["G01"].read_text(encoding="utf-8"))
    tax=json.loads(paths["G02"].read_text())
    if not isinstance(tax,dict) or HUMAN not in tax:
        raise ValueError("Unexpected taxon index")
    memberships=defaultdict(set)
    for k,v in tax.items():
        if not isinstance(v,list) or any(not isinstance(i,str) or not ID.fullmatch(i) for i in v):
            raise ValueError(f"Unexpected index entries {k}")
        for i in v:
            memberships[i].add(k)
    raw=tax[HUMAN]
    ids=sorted(set(raw))
    missing=set(ids)-idx.ids
    print("HTML_UNIQUE_IDS",len(idx.ids),"TAXON_KEYS",len(tax),
          "TAXON_UNION",len(memberships))
    print("HUMAN_INDEX_ENTRIES",len(raw),"HUMAN_UNIQUE_IDS",len(ids),
          "HUMAN_DUPLICATE_ENTRIES",len(raw)-len(ids),
          "HUMAN_NOT_IN_HTML",len(missing))
    print("HTML_NOT_IN_TAXON_INDEX",len(idx.ids-set(memberships)),
          "TAXON_INDEX_NOT_IN_HTML",len(set(memberships)-idx.ids))
    print("HUMAN_INDEX_MULTITAXON",sum(len(memberships[i])>1 for i in ids))
    print("HUMAN_INDEX_ADDITIONAL_TAXA",dict(Counter(k for i in ids for k in memberships[i] if k!=HUMAN)))
    if len(raw)!=1520 or len(idx.ids)!=2078 or len(tax)!=67 or missing or len(raw)!=len(ids):
        raise RuntimeError("Unexpected human index cardinality or mismatch; review before download")
    return ids,memberships

def inspect(p, ident):
    obj=json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(obj,dict) or obj.get("id")!="gomodel:"+ident:
        raise ValueError(f"Model identity mismatch {ident}")
    def count(k):
        return len(obj[k]) if isinstance(obj.get(k),list) else None
    return {"model_id":ident,"model_taxon":obj.get("taxon"),
            "model_status":obj.get("status"),"date_modified":obj.get("date_modified"),
            "top_level_keys":sorted(obj),
            **{k+"_n":count(k) for k in ("activities","molecules","objects","provenances")},
            "taxon_exact_human":obj.get("taxon")==HUMAN}

def acquire(ident, memberships):
    p=DIRECTORY/(ident+".json")
    rec=DIRECTORY/(ident+".receipt.json")
    url=BASE+ident+".json"
    if p.exists() or rec.exists():
        if not p.is_file() or not rec.is_file():
            raise RuntimeError(f"Unpaired file/receipt {ident}; inspect")
        m=json.loads(rec.read_text())
        if (m.get("url")!=url or m.get("sha256")!=sha256(p)
                or m.get("bytes")!=p.stat().st_size
                or m.get("indexed_taxa")!=sorted(memberships)
                or m.get("model")!=inspect(p,ident)):
            raise RuntimeError(f"Receipt verification failed {ident}")
        return m,False
    # A failed GET leaves no final file; a crash during publication is detected above.
    with tempfile.TemporaryDirectory(prefix=".transfer-",dir=DIRECTORY) as tmpdir:
        temp=Path(tmpdir)/"model.json"
        info=subprocess.check_output(curl_base()+[
            "--max-time","180","--max-filesize","50000000",
            "--output",str(temp),"--write-out","%{http_code}\t%{url_effective}",url],
            text=True).split("\t",1)
        if len(info)!=2 or info[0]!="200":
            raise RuntimeError(f"Unexpected HTTP response {ident}: {info}")
        summary=inspect(temp,ident)
        m={"url":url,"effective_url":info[1],"http_status":200,
           "release_archive":"2026-08-05","bytes":temp.stat().st_size,
           "sha256":sha256(temp),"acquired_utc":utc(),
           "indexed_taxa":sorted(memberships),"model":summary,
           "status":"ACQUIRED_NOT_CERTIFIED"}
        os.replace(temp,p)
        atomic_json(rec,m)
    time.sleep(0.25)
    return m,True

def summarize(receipts):
    models=[r["model"] for r in receipts]
    print("MODEL_TAXA",dict(Counter(json.dumps(m["model_taxon"],sort_keys=True) for m in models)))
    print("MODEL_STATUS",dict(Counter(json.dumps(m["model_status"],sort_keys=True) for m in models)))
    print("MODEL_NONEXACT_HUMAN_TAXON",sum(not m["taxon_exact_human"] for m in models))
    print("MULTITAXON_INDEX_MODELS",sum(len(r["indexed_taxa"])>1 for r in receipts))
    print("TOP_LEVEL_KEY_COUNTS",dict(Counter(k for m in models for k in m["top_level_keys"])))
    for k in ("activities","molecules","objects","provenances"):
        vals=[m[k+"_n"] for m in models]
        print("MODEL_FIELD",k,"total_entries",sum(v for v in vals if v is not None),
              "missing_or_nonlist",sum(v is None for v in vals),
              "empty",sum(v==0 for v in vals))
    print("TOTAL_MODEL_BYTES",sum(r["bytes"] for r in receipts))

def self_test():
    idx=Index()
    idx.feed('<a href="abc123.json">model</a><a href="index.html">index</a>')
    assert idx.ids=={"abc123"}
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"test.json"
        p.write_text(json.dumps({"id":"gomodel:abc123","taxon":HUMAN,"activities":[{}]}))
        m=inspect(p,"abc123")
        assert m["taxon_exact_human"] and m["activities_n"]==1 and m["objects_n"] is None
        try:
            inspect(p,"abc124")
        except ValueError:
            pass
        else:
            raise AssertionError("ID mismatch not detected")
    print("SELF_TEST_PASS (index and model identity fixtures only)")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--execute",action="store_true")
    ap.add_argument("--pilot",action="store_true",help="first ten sorted IDs only")
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return
    ids,memberships=load_indices()
    selected=ids[:10] if args.pilot else ids
    print("CFA_P1A03B_HUMAN_INDEXED_GO_CAM_V1")
    print("MODE","EXECUTE" if args.execute else "PLAN_ONLY","SELECTED",len(selected))
    if not args.execute:
        return
    DIRECTORY.mkdir(parents=True,exist_ok=True)
    # Prevent simultaneous instances; crash leaves lock requiring inspection.
    lock=DIRECTORY/".acquisition.lock"
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    os.close(fd)
    receipts=[]
    try:
        expected={i+".json" for i in ids}|{i+".receipt.json" for i in ids}
        extras={p.name for p in DIRECTORY.glob("*.json")}-expected-{"MODEL_INVENTORY.json"}
        if extras:
            raise RuntimeError(f"Unexpected local model files: {sorted(extras)[:10]}")
        # Batch transfers to limit queued requests after a failure.
        with ThreadPoolExecutor(max_workers=3) as pool:
            for start in range(0,len(selected),12):
                batch=selected[start:start+12]
                results=list(pool.map(lambda i: acquire(i,memberships[i]),batch))
                receipts.extend(m for m,_ in results)
                print("PROGRESS",len(receipts),"/",len(selected),
                      "new_in_batch",sum(new for _,new in results),flush=True)
        summarize(receipts)
        if not args.pilot:
            path=DIRECTORY/"MODEL_INVENTORY.json"
            if path.exists():
                if json.loads(path.read_text())!=receipts:
                    raise RuntimeError("Existing inventory mismatch")
            else:
                atomic_json(path,receipts)
            print("GO_CAM_ACQUIRED",len(receipts),"SOURCE_CERTIFICATION_PENDING")
        else:
            print("PILOT_ACQUIRED",len(receipts),"FULL_COLLECTION_PENDING")
    finally:
        lock.unlink()

if __name__=="__main__":
    main()
