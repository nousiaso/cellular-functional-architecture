#!/usr/bin/env python3
"""Read-only publication and assertion provenance reconciliation."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import time

ROOT=Path(__file__).resolve().parents[2]
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
TYPED=ROOT/"results/P1A04F_TYPED_SOURCE_LAYER_V1"
START=time.monotonic()
def emit(*x): print(*x,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def norm_pmid(x):
    if not isinstance(x,str): return None
    m=re.fullmatch(r"PMID:\s*([0-9]+)",x.strip())
    return m.group(1) if m else None
def main():
    with (TYPED/"biopax_operations.tsv").open() as f:
        events={r["rdf_id"]:r["stable_id"] for r in csv.DictReader(f,delimiter="\t")}
    need(len(events)==16187,"event count")
    universe=set(events.values())
    recs=json.loads((ROOT/"data/raw/P1A03_CORE_RECEIPTS.json").read_text())
    matches=[r for r in recs if r.get("id")=="R06"]
    need(len(matches)==1,"R06 receipt cardinality")
    r=matches[0]; p=(ROOT/r["path"]).resolve()
    need(p.is_relative_to((ROOT/"data/raw").resolve()),"R06 path")
    data=p.read_bytes()
    need(len(data)==r["bytes"]==1136864 and hashlib.sha256(data).hexdigest()==r["sha256"]=="3b2169fc94cbe51b3eea84e28c11c37512ebd834c7bca97f18719f0027a5abf8","R06 integrity")
    raw=set(); numeric=defaultdict(set); atypical=[]
    for row in csv.reader(data.decode("utf-8").splitlines(),delimiter="\t"):
        if not row or row[0] not in universe: continue
        raw.add(row[0])
        if len(row)==2 and re.fullmatch(r"[0-9]+",row[1]):
            numeric[row[0]].add(row[1])
        else: atypical.append(row)
    emit("R06_EVENT_COUNTS","raw",len(raw),"numeric",len(numeric),"numeric_missing",sorted(raw-numeric.keys()))
    emit("R06_ATYPICAL_ROWS",len(atypical),atypical[:30])
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","schema")
    # Join direct xref properties only. Do not infer citations from pathway ancestry.
    bp=defaultdict(set); pub_other=Counter(); count=0
    sql="""SELECT e.rdf_id,x.rdf_id,p.predicate_uri,p.literal_text
    FROM bp_property l JOIN bp_record e ON e.ordinal=l.record_ordinal
    JOIN bp_record x ON x.rdf_id=l.target_local_id
    JOIN bp_property p ON p.record_ordinal=x.ordinal
    WHERE l.predicate_uri=? AND x.class_uri=?
    AND p.predicate_uri IN (?,?)"""
    prefix="{http://www.biopax.org/release/biopax-level3.owl#}"
    by_link=defaultdict(lambda:defaultdict(set))
    for event,xref,pred,value in db.execute(sql,(prefix+"xref",prefix+"PublicationXref",prefix+"db",prefix+"id")):
        if event in events: by_link[(event,xref)][pred.rsplit("}",1)[-1]].add(value)
    for (event,xref),fields in by_link.items():
        namespaces=fields["db"]; ids=fields["id"]
        if any(isinstance(x,str) and x.lower()=="pubmed" for x in namespaces):
            bp[events[event]].update(x for x in ids if isinstance(x,str) and re.fullmatch(r"[0-9]+",x))
        else: pub_other[tuple(sorted(str(x) for x in namespaces))]+=1
    bp={k:v for k,v in bp.items() if v}
    emit("DIRECT_BP_PUBLICATION_XREF_LINKS",len(by_link),"OTHER_NAMESPACE",dict(pub_other))
    emit("BP_R06_PMID_EVENT_COVERAGE","direct_bp",len(bp),"numeric_export",len(numeric),
         "both",len(bp.keys()&numeric.keys()),"direct_only",len(bp.keys()-numeric.keys()),
         "export_only",len(numeric.keys()-bp.keys()),"neither",len(universe-(bp.keys()|numeric.keys())))
    tags=Counter()
    for e in bp.keys()&numeric.keys():
        a,b=bp[e],numeric[e]
        tags["equal" if a==b else "direct_subset" if a<b else "export_subset" if b<a else "partial" if a&b else "disjoint"]+=1
    emit("BP_R06_PMID_SET_COMPARISON",dict(tags))
    # Activity payload evidence includes nested enabler and causal assertions.
    # Separately inspect association evidence directly attached to each staged edge.
    def source(mid): return "Reactome_namespace" if mid.startswith("R-HSA-") else "other_namespace"
    def walk(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k=="reference" and isinstance(v,str):yield v
                yield from walk(v)
        elif isinstance(x,list):
            for v in x:yield from walk(v)
    def summarize(refs):
        return (bool(refs),any(norm_pmid(x) for x in refs),
                any(x.strip().startswith("Reactome:") for x in refs))
    ac=Counter(); edges=Counter(); eco=Counter(); pubs=defaultdict(set)
    for mid,payload in db.execute("SELECT model_id,payload_json FROM go_record WHERE field='activities'"):
        refs=set(walk(json.loads(payload))); ac[(source(mid),*summarize(refs))]+=1
        pubs[source(mid)].update(x for s in refs if (x:=norm_pmid(s)))
    for mid,field,pred,payload in db.execute("SELECT model_id,field,predicate,payload_json FROM go_association"):
        obj=json.loads(payload); ev=obj.get("evidence")
        refs=set(walk(ev)); edges[(source(mid),field,*summarize(refs))]+=1
        if isinstance(ev,list):
            for item in ev:
                if isinstance(item,dict):eco[(source(mid),field,str(item.get("term")))]+=1
    emit("GO_ACTIVITY_REFERENCE_PROFILES",sorted(ac.items()))
    emit("GO_ASSOCIATION_DIRECT_EVIDENCE_PROFILES",sorted(edges.items()))
    emit("GO_ASSOCIATION_ECO_OCCURRENCES",sorted(eco.items()))
    all_bp=set().union(*bp.values(),*numeric.values())
    for k,v in sorted(pubs.items()):
        emit("GO_ACTIVITY_PMID_UNIVERSE",k,len(v),"overlap_bp_or_export",len(v&all_bp))
    db.close()
    emit("P1A04J_PROVENANCE_AUDIT_COMPLETE; NO_MERGES_OR_EVIDENCE_WEIGHTS")
if __name__=="__main__": main()
