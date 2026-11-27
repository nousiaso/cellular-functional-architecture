#!/usr/bin/env python3
"""Inventory chemical identifiers and optionally acquire/audit one ChEBI FULL snapshot."""
import argparse
from collections import Counter, defaultdict
import datetime as dt
import gzip
import json
from pathlib import Path
import re
import sqlite3
import subprocess
import tempfile
from p1a04k_scope_feasibility_audit import ROOT,STATE,rows,sha,need,check_manifest

DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
URL="https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"
OUT=ROOT/"data/raw/ontology_interpretation/ChEBI/monthly_snapshot"
def emit(*x):print(*x,flush=True)
def chebi_id(s):
    m=re.fullmatch(r"(?:CHEBI:)?([0-9]+)",s,re.I)
    return "CHEBI:"+m.group(1) if m else None
def inventory():
    check_manifest(STATE,"P1A04H_v1")
    need(DB.stat().st_size==1451589632,"staging size")
    bp=set();bad=set();raw=defaultdict(set);n=0
    for r in rows(STATE/"biopax_states.tsv"):
        if r["class"]!="SmallMolecule":continue
        n+=1
        for ref,xs in json.loads(r["reference_xrefs_json"]):
            for ns,acc in xs:
                raw[ns].add(acc)
                if ns=="ChEBI":
                    x=chebi_id(acc)
                    if x:bp.add(x)
                    else:bad.add(acc)
    need(n==5451,"chemical state inventory")
    emit("BP_CHEMICAL_REFERENCE_TOKENS",[(k,len(v),sorted(v)[:5]) for k,v in sorted(raw.items())])
    emit("BP_CHEBI_NONSTANDARD_TOKENS",sorted(bad))
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","staging schema")
    go=set(); profiles=Counter(); key_profiles=Counter(); examples=defaultdict(list)
    def walk(x,path=()):
        if isinstance(x,dict):
            for k,v in x.items():
                if k=="term" and isinstance(v,str):yield ".".join(path+(k,)),v
                yield from walk(v,path+(k,))
        elif isinstance(x,list):
            for v in x:yield from walk(v,path+("[]",))
    for field in ("molecules","activities"):
        for mid,sid,payload in db.execute("SELECT model_id,source_id,payload_json FROM go_record WHERE field=?",(field,)):
            obj=json.loads(payload)
            key_profiles[(field,tuple(sorted(obj)))]+=1
            for loc,t in walk(obj):
                ns=t.split(":",1)[0]
                profiles[(field,loc,ns)]+=1
                if ns=="CHEBI":
                    x=chebi_id(t)
                    need(x is not None,"unrecognized GO ChEBI token "+t)
                    go.add(x)
                if len(examples[(field,loc,ns)])<2:examples[(field,loc,ns)].append((mid,sid,t))
    db.close()
    emit("GO_RECORD_KEY_PROFILES",sorted(key_profiles.items()))
    emit("GO_TERM_NAMESPACE_BY_FIELD_PATH",sorted(profiles.items()))
    emit("GO_TERM_EXAMPLES",sorted(examples.items()))
    emit("CHEBI_DISTINCT_IDS","BP",len(bp),"GO",len(go),"shared",len(bp&go),"union",len(bp|go))
    return bp,go
def inspect(path,wanted):
    header=[]; found=defaultdict(list); aliases=defaultdict(list); stanza=None; fields=defaultdict(list); n=0
    def finish():
        nonlocal n
        if stanza!="[Term]" or not fields.get("id"):return
        ident=fields["id"][0]
        n+=1
        hits=set(fields.get("alt_id",[]))&wanted
        if ident in wanted or hits:
            item={k:fields[k] for k in ("id","name","is_obsolete","alt_id","replaced_by","consider","xref","property_value","is_a") if k in fields}
            if ident in wanted:found[ident].append(item)
            for x in hits:aliases[x].append(item)
    with gzip.open(path,"rt",encoding="utf-8") as f:
        for line in f:
            s=line.rstrip("\r\n")
            if s.startswith("[") and s.endswith("]"):
                finish();stanza=s;fields=defaultdict(list);continue
            if stanza is None:header.append(s);continue
            if stanza=="[Term]" and ": " in s:
                k,v=s.split(": ",1)
                if k in ("id","name","is_obsolete","alt_id","replaced_by","consider","xref","property_value","is_a"):
                    fields[k].append(v)
        finish()
    version=[s for s in header if s.startswith("data-version:")]
    emit("CHEBI_HEADER",header[:30])
    need(len(version)==1 and n>0,"ChEBI version or term schema not established")
    emit("CHEBI_TERM_COUNT",n)
    return version[0],found,aliases
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--execute",action="store_true",help="download and verify the official monthly FULL OBO snapshot")
    a=ap.parse_args()
    bp,go=inventory()
    if not a.execute:
        emit("PLAN_ONLY",URL,"NO_DOWNLOAD")
        return
    wanted=bp|go
    if OUT.exists():
        p=OUT/"chebi.obo.gz";r=json.loads((OUT/"receipt.json").read_text())
        need(r["url"]==URL and p.stat().st_size==r["bytes"] and sha(p)==r["sha256"],"ChEBI receipt mismatch")
        version,found,aliases=inspect(p,wanted)
        need(version==r["data_version"],"ChEBI version mismatch")
    else:
        OUT.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="chebi_snapshot_",dir=OUT.parent) as td:
            temp=Path(td);p=temp/"chebi.obo.gz"
            result=subprocess.run(["curl","--fail","--location","--silent","--show-error","--retry","3","--connect-timeout","30","--max-time","900","--proto","=https","--proto-redir","=https","--output",str(p),"--write-out","%{url_effective}",URL],check=True,capture_output=True,text=True)
            version,found,aliases=inspect(p,wanted)
            r={"url":URL,"effective_url":result.stdout,"data_version":version,
               "bytes":p.stat().st_size,"sha256":sha(p),"acquired_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
               "license_reference":"https://www.ebi.ac.uk/chebi/about/","scope":"monthly FULL ontology; release observed from payload"}
            (temp/"receipt.json").write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
            need(not OUT.exists(),"snapshot destination appeared")
            temp.rename(OUT)
    emit("CHEBI_SNAPSHOT_RECEIPT",json.dumps(r,sort_keys=True))
    def classify(x):
        f=found.get(x,[]);al=aliases.get(x,[])
        if len(f)==1 and not al:return "exact_obsolete" if "true" in f[0].get("is_obsolete",[]) else "exact_active"
        if not f and len(al)==1:return "alternate_id_candidate"
        if not f and not al:return "not_found"
        return "ambiguous"
    for name,ids in (("BP",bp),("GO",go),("UNION",wanted)):
        c=Counter(classify(x) for x in ids)
        emit("CHEBI_ID_COVERAGE",name,dict(c))
        emit("CHEBI_REVIEW_IDS",name,[(x,classify(x),[(v.get("id"),v.get("name"),v.get("replaced_by"),v.get("consider")) for v in found.get(x,[])+aliases.get(x,[])]) for x in sorted(ids) if classify(x)!="exact_active"][:40])
    emit("P1A04P_CHEMICAL_REFERENCE_AUDIT_COMPLETE; NO_MERGES")
if __name__=="__main__":main()
