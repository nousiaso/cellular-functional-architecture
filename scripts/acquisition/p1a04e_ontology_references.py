#!/usr/bin/env python3
"""Acquire immutable RO/ECO reference snapshots; audit observed GO-CAM terms."""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import os
import sqlite3
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/"data/raw"
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
OUT=RAW/"ontology_interpretation"
SOURCES={
 "RO":("oborel/obo-relations","091f4761ec84657cc336d623e02006d6c351722c","ro.obo","d43af196c85721f9ee7a0547c1a4d8bdb5fe8830","2026-09-04"),
 "ECO":("evidenceontology/evidenceontology","6bbc0730de4f3d04a5a1347946cd231e40fde4bd","eco.obo","13e647b9851f76caba2c06a0dc27703873d766bb","2026-07-10")
}
START=time.monotonic()
def emit(*args):
    print(*args,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)
def parse_obo(data):
    header=defaultdict(list)
    entries={}
    current=header
    def save(stanza):
        if stanza is header or "id" not in stanza:
            return
        ident=stanza["id"][0]
        need(ident not in entries,"duplicate OBO stanza ID "+ident)
        entries[ident]=stanza
    for line in data.decode("utf-8-sig").splitlines():
        line=line.strip()
        if not line or line.startswith("!"):
            continue
        if line.startswith("[") and line.endswith("]"):
            save(current)
            current=defaultdict(list)
            current["_stanza"].append(line)
        elif ":" in line:
            key,val=line.split(":",1)
            current[key].append(val.strip())
    save(current)
    return header,entries
def get_source(name,spec,execute):
    repo,commit,filename,blob,date=spec
    url=f"https://raw.githubusercontent.com/{repo}/{commit}/{filename}"
    p=OUT/name/date/filename
    rec=p.with_name(p.name+".receipt.json")
    emit("SOURCE",name,date,commit,blob,url)
    if not execute:
        return None
    p.parent.mkdir(parents=True,exist_ok=True)
    def inspect(data):
        h=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
        need(h==blob,"Git blob hash mismatch "+name)
        head,entries=parse_obo(data)
        need(head.get("data-version")==["releases/"+date],"ontology version mismatch "+name)
        need(entries,"no OBO stanzas "+name)
        receipt={"name":name,"release":date,"commit":commit,"git_blob_sha1":blob,
                 "url":url,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}
        return receipt,entries
    if p.exists() or rec.exists():
        need(p.is_file() and rec.is_file(),"unpaired ontology files "+name)
        receipt,entries=inspect(p.read_bytes())
        need(json.loads(rec.read_text())==receipt,"receipt mismatch "+name)
        emit("REUSE",name,receipt)
        return entries
    with tempfile.TemporaryDirectory(dir=p.parent,prefix=".transfer-") as tmp:
        tmp=Path(tmp)
        q=tmp/filename
        subprocess.run(["curl","--fail","--silent","--show-error","--location",
                        "--proto","=https","--proto-redir","=https","--retry","3",
                        "--connect-timeout","30","--max-time","180","--max-filesize","30000000",
                        "--output",str(q),url],check=True)
        receipt,entries=inspect(q.read_bytes())
        r=tmp/"receipt.json"
        r.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
        # Exclusive publication, no overwrite; unpaired artifact fails closed on rerun.
        os.link(q,p)
        os.link(r,rec)
        emit("ACQUIRED",name,receipt)
        return entries
def audit_terms():
    need(DB.stat().st_size==1451589632,"unexpected staged database size")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    meta=dict(db.execute("SELECT key,value FROM metadata"))
    need(meta.get("schema_version")=="P1A04C_v1","staged schema mismatch")
    ro=Counter()
    for field,pred,n in db.execute("SELECT field,predicate,COUNT(*) FROM go_association GROUP BY field,predicate"):
        ro[(field,pred)]=n
    need(sum(ro.values())==15491 and len({p for f,p in ro})==15,"unexpected staged association inventory")
    eco=Counter()
    go=Counter()
    def walk(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k=="term" and isinstance(v,str):
                    if v.startswith("ECO:"): eco[v]+=1
                    if v.startswith("GO:"): go[v]+=1
                walk(v)
        elif isinstance(x,list):
            for v in x: walk(v)
    n=0
    for (payload,) in db.execute("SELECT payload_json FROM go_model"):
        walk(json.loads(payload))
        n+=1
    need(n==1520,"unexpected model count")
    db.close()
    emit("OBSERVED_RO_PREDICATES",sorted(ro.items()))
    emit("OBSERVED_ECO_TERM_OCCURRENCES",sorted(eco.items()))
    emit("OBSERVED_GO_TERM_IDS",len(go))
    return ro,eco,go
def resolve(label,observed,entries,details):
    alt=defaultdict(list)
    for ident,entry in entries.items():
        for aid in entry.get("alt_id",[]):
            alt[aid.split()[0]].append(ident)
    counts=Counter()
    for ident,n in sorted(observed.items()):
        if ident in entries:
            target=ident
            status="exact_obsolete" if entries[target].get("is_obsolete")==["true"] else "exact_active"
        elif len(alt[ident])==1:
            target=alt[ident][0]
            status="alternate_to_obsolete" if entries[target].get("is_obsolete")==["true"] else "alternate_to_active"
        else:
            target=None
            status="ambiguous_alt" if alt[ident] else "unresolved"
        counts[status]+=1
        if details or status!="exact_active":
            e=entries.get(target,{})
            emit("TERM",label,ident,"OCCURRENCES",n,"STATUS",status,"TARGET",target,
                 "STANZA",e.get("_stanza"),"NAME",e.get("name"),"DEFINITION",e.get("def"),
                 "IS_A",e.get("is_a"),"RELATIONSHIP",e.get("relationship"),
                 "REPLACED_BY",e.get("replaced_by"),"CONSIDER",e.get("consider"))
    emit("TERM_RESOLUTION_SUMMARY",label,"DISTINCT",len(observed),dict(counts))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--execute",action="store_true")
    a=ap.parse_args()
    emit("CFA_P1A04E_ONTOLOGY_REFERENCES_V1", "EXECUTE" if a.execute else "PLAN_ONLY")
    snapshots={k:get_source(k,v,a.execute) for k,v in SOURCES.items()}
    if not a.execute:
        return
    ro,eco,go=audit_terms()
    rc=Counter()
    for (field,pred),n in ro.items():
        rc[pred]+=n
    resolve("RO",rc,snapshots["RO"],True)
    resolve("ECO",eco,snapshots["ECO"],True)
    p=RAW/"go_2026-08-05/ontology_identity/go.obo"
    raw=p.read_bytes()
    need(hashlib.sha256(raw).hexdigest()=="d3593751d885ca160b2ab7baf6c7eccd88ca3c4599f79436c674bad661095ff0","GO ontology checksum mismatch")
    h,entries=parse_obo(raw)
    emit("GO_DATA_VERSION",h.get("data-version"))
    resolve("GO",go,entries,False)
    emit("P1A04E_ONTOLOGY_REFERENCE_AUDIT_COMPLETE; REVIEW_TERM_SEMANTICS_BEFORE_GRAPH")
if __name__=="__main__":
    main()
