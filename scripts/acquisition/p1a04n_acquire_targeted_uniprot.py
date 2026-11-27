#!/usr/bin/env python3
"""Acquire release-attested UniProt identity records for exactly the unresolved GO enabler tokens."""
import argparse
import base64
from collections import Counter, defaultdict
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"diagnostics"))
from p1a04k_scope_feasibility_audit import ROOT, STATE, rows, check_manifest, need
from p1a04m_uniprot_identity_pilot import parent, BASE

RELEASE="2026_04"
RELEASE_DATE="04-November-2026"
OUT=ROOT/"data/raw/uniprot_2026_04/targeted_go_enablers"
TOKEN=re.compile(r"[A-Z0-9]+(?:-(?:[0-9]+|PRO_[0-9]+))?")
def emit(*a): print(*a,flush=True)
def digest(b): return hashlib.sha256(b).hexdigest()
def atomic(path,obj):
    tmp=path.with_name(path.name+".partial")
    with tmp.open("x",encoding="utf-8") as f:
        json.dump(obj,f,sort_keys=True,ensure_ascii=False,indent=2)
        f.write("\n");f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
def tokens():
    check_manifest(STATE,"P1A04H_v1")
    bt=defaultdict(set); unknown=set()
    for r in rows(STATE/"biopax_states.tsv"):
        if r["class"]!="Protein" or not json.loads(r["reference_ids_json"]):continue
        for _,xs in json.loads(r["reference_xrefs_json"]):
            for ns,acc in xs:
                if ns in ("UniProt","UniProt Isoform"):
                    bt[acc].update(json.loads(r["taxon_union_json"]))
                    if r["unresolved_taxon"]=="1":unknown.add(acc)
    def terms(obj):
        if isinstance(obj,dict):
            t=obj.get("term")
            if isinstance(t,str) and t.startswith("UniProtKB:"):yield t.split(":",1)[1]
            yield from terms(obj.get("has_part",[]))
        elif isinstance(obj,list):
            for x in obj:yield from terms(x)
    go=set()
    for r in rows(STATE/"gocam_enablers.tsv"):go.update(terms(json.loads(r["enabler_json"])))
    result=sorted(x for x in go if not bt.get(x) or x in unknown)
    need(len(result)==520,"unresolved inventory differs from P1A04L")
    need(all(TOKEN.fullmatch(x) for x in result),"unexpected token syntax")
    return result
def identity(token,obj):
    acc=parent(token); primary=obj.get("primaryAccession")
    org=obj.get("organism") or {}
    need(isinstance(primary,str) and isinstance(org.get("taxonId"),int),"identity schema")
    status="primary" if acc==primary else "secondary" if acc in obj.get("secondaryAccessions",[]) else "unresolved_returned_identity"
    isoforms={x for c in obj.get("comments",[]) if c.get("commentType")=="ALTERNATIVE PRODUCTS"
              for i in c.get("isoforms",[]) for x in i.get("isoformIds",[])}
    features=[f for f in obj.get("features",[]) if f.get("featureId") and token.endswith("-"+f["featureId"])]
    if token==acc: sub="parent_token"
    elif "-PRO_" in token:sub="feature_exact" if len(features)==1 else "feature_unresolved"
    else:sub="isoform_exact" if token in isoforms else "isoform_unresolved"
    hgnc=sorted({r["id"] for r in obj.get("uniProtKBCrossReferences",[])
                 if r.get("database")=="HGNC" and isinstance(r.get("id"),str)})
    return {"token":token,"queried_parent":acc,"primary":primary,"identity":status,
            "substate":sub,"taxon_id":org["taxonId"] if status!="unresolved_returned_identity" else None,
            "scientific_name":org.get("scientificName"),"entry_type":obj.get("entryType"),
            "gene_names":[g.get("geneName",{}).get("value") for g in obj.get("genes",[])],
            "HGNC_xrefs":hgnc,"matching_features":[{"id":f["featureId"],"type":f.get("type"),"location":f.get("location")} for f in features]}
def fetch(acc):
    url=BASE+acc+".json"
    for attempt in range(4):
        req=Request(url,headers={"Accept":"application/json","User-Agent":"CFA-identity-acquisition/1.0"})
        try:
            with urlopen(req,timeout=45) as r:
                status=r.status; body=r.read(); headers=r.headers; final=r.geturl()
        except HTTPError as e:
            if e.code in (404,410):
                status=e.code;body=e.read();headers=e.headers;final=e.geturl()
            elif e.code in (429,500,502,503,504) and attempt<3:
                delay=e.headers.get("Retry-After","")
                time.sleep(min(120,int(delay)) if delay.isdigit() else 2**(attempt+1))
                continue
            else:raise
        except (URLError,TimeoutError):
            if attempt==3:raise
            time.sleep(2**(attempt+1));continue
        need(urlsplit(final).scheme=="https" and urlsplit(final).hostname=="rest.uniprot.org","unexpected endpoint")
        release=headers.get("X-UniProt-Release"); date=headers.get("X-UniProt-Release-Date")
        if status==200:
            need((release,date)==(RELEASE,RELEASE_DATE),"release drift or missing header")
            identity(acc,json.loads(body))
        else:
            # Negative lookups with no release headers remain unresolved, not certified absences.
            need(release in (None,RELEASE) and date in (None,RELEASE_DATE),"release drift on error response")
        return {"accession":acc,"url":url,"effective_url":final,"status":status,
                "release":release,"release_date":date,
                "acquired_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
                "bytes":len(body),"sha256":digest(body),"body_base64":base64.b64encode(body).decode("ascii")}
    raise RuntimeError("retry exhausted")
def verify(acc,rec):
    need(rec["accession"]==acc and rec["url"]==BASE+acc+".json","receipt identity")
    u=urlsplit(rec["effective_url"]);need(u.scheme=="https" and u.hostname=="rest.uniprot.org","receipt endpoint")
    body=base64.b64decode(rec["body_base64"],validate=True)
    need(len(body)==rec["bytes"] and digest(body)==rec["sha256"],"receipt payload hash")
    if rec["status"]==200:
        need((rec["release"],rec["release_date"])==(RELEASE,RELEASE_DATE),"receipt release")
        obj=json.loads(body);identity(acc,obj);return obj
    need(rec["status"] in (404,410),"unexpected stored status")
    need(rec["release"] in (None,RELEASE) and rec["release_date"] in (None,RELEASE_DATE),"error receipt release")
    return None
def self_test():
    obj={"primaryAccession":"P12345","secondaryAccessions":["Q12345"],"organism":{"taxonId":9606},
         "comments":[{"commentType":"ALTERNATIVE PRODUCTS","isoforms":[{"isoformIds":["P12345-2"]}]}],
         "features":[{"featureId":"PRO_0000000001","type":"Chain","location":{"start":{"value":10}}}]}
    need(identity("P12345-2",obj)["substate"]=="isoform_exact","isoform fixture")
    need(identity("P12345-3",obj)["substate"]=="isoform_unresolved","missing isoform fixture")
    need(identity("P12345-PRO_0000000001",obj)["substate"]=="feature_exact","feature fixture")
    need(identity("Q12345",obj)["identity"]=="secondary","secondary fixture")
    need(identity("Q99999",obj)["taxon_id"] is None,"unconfirmed identity fixture")
    emit("SELF_TEST_PASS")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--self-test",action="store_true");ap.add_argument("--execute",action="store_true")
    a=ap.parse_args()
    if a.self_test:self_test();return
    ts=tokens(); parents=sorted({parent(t) for t in ts})
    plan={"schema_version":"P1A04N_v1","release":RELEASE,"release_date":RELEASE_DATE,"tokens":ts,"parents":parents}
    emit("ACQUISITION_PLAN","tokens",len(ts),"distinct_parent_requests",len(parents),"release",RELEASE)
    if not a.execute:emit("PLAN_ONLY");return
    OUT.mkdir(parents=True,exist_ok=True)
    pp=OUT/"PLAN.json"
    if pp.exists():need(json.loads(pp.read_text())==plan,"stored plan mismatch")
    else:atomic(pp,plan)
    records=OUT/"records";records.mkdir(exist_ok=True)
    need(not list(OUT.rglob("*.partial")),"partial artifact requires inspection")
    need(set(x.name for x in records.iterdir())<={x+".json" for x in parents},"unexpected record file")
    by_parent={}; hashes={}; summary=Counter()
    for i,acc in enumerate(parents,1):
        path=records/(acc+".json")
        if path.exists():rec=json.loads(path.read_text());summary["reused"]+=1
        else:
            rec=fetch(acc);verify(acc,rec);atomic(path,rec);summary["fetched"]+=1
            time.sleep(0.25)
        by_parent[acc]=verify(acc,rec)
        hashes[acc]={"status":rec["status"],"bytes":rec["bytes"],"payload_sha256":rec["sha256"],
                     "receipt_sha256":digest(path.read_bytes())}
        if i%25==0 or i==len(parents):emit("PARENT_PROGRESS",i,len(parents),dict(summary))
    result=[]; counts=Counter(); taxa=Counter()
    for token in ts:
        acc=parent(token);obj=by_parent[acc]
        row=identity(token,obj) if obj is not None else {"token":token,"queried_parent":acc,"identity":"lookup_unresolved","substate":"unresolved","taxon_id":None}
        result.append(row);counts[(row["identity"],row["substate"])]+=1
        t=row["taxon_id"];taxa["human" if t==9606 else "unresolved" if t is None else "nonhuman"]+=1
    resolved={"schema_version":"P1A04N_v1","release":RELEASE,"token_records":result,"parent_receipts":hashes}
    rp=OUT/"IDENTITY_INVENTORY.json"
    if rp.exists():need(json.loads(rp.read_text())==resolved,"stored identity inventory mismatch")
    else:atomic(rp,resolved)
    emit("IDENTITY_STATUS",sorted(counts.items()))
    emit("TAXONOMY_STATUS",dict(taxa))
    emit("UNRESOLVED_EXACT_TOKENS",[r["token"] for r in result if r["identity"] not in ("primary","secondary") or r["substate"] not in ("parent_token","isoform_exact","feature_exact")])
    emit("IDENTITY_INVENTORY_SHA256",digest(rp.read_bytes()))
    emit("P1A04N_TARGETED_UNIPROT_ACQUISITION_COMPLETE; NO_CANONICAL_MERGES")
if __name__=="__main__":main()
