#!/usr/bin/env python3
"""Final bounded source receipt and scope preflight. No DB table scans or new downloads."""
from collections import Counter,defaultdict
from functools import lru_cache
import csv
import json
import re
import sqlite3
from p1a04s_remaining_identity_and_format_preflight import core,ROOT,STATE,DB,rows,sha,need
from p1a04k_scope_feasibility_audit import TYPED,check_manifest
from p1a04q_state_correspondence_audit import CHEBI_HASH,UNI_HASH
from p1a04s_remaining_identity_and_format_preflight import PART_HASH
def emit(*a):print(*a,flush=True)
STABLE=re.compile(r"R-[A-Z]+-[0-9]+")
def pairs(path):
    with path.open(newline="",encoding="utf-8") as f:
        for i,r in enumerate(csv.reader(f,delimiter="\t"),1):
            if not r:continue
            need(len(r)==2,"unexpected two-column schema "+path.name+" row "+str(i))
            yield r
def main():
    check_manifest(TYPED,"P1A04F_v1")
    check_manifest(STATE,"P1A04H_v1")
    raw=ROOT/"data/raw"
    pins=list(rows(ROOT/"manifests/P1A03_CORE_REPORTED_SHA256.tsv"))
    receipts=json.loads((raw/"P1A03_CORE_RECEIPTS.json").read_text())
    need(len(pins)==len(receipts)==19,"core source count")
    byid={r["id"]:r for r in receipts}
    need(set(byid)=={r["id"] for r in pins},"core receipt IDs")
    for x in pins:
        p=core(x["id"]);r=byid[x["id"]]
        side=p.with_name(p.name+".receipt.json")
        need(json.loads(side.read_text())==r,"core sidecar mismatch "+x["id"])
        emit("CORE_SOURCE_VERIFIED",x["id"],r["bytes"],r["sha256"],r["release_candidate"],r["role"])
    # Recheck pinned inventories; payload-level checks already passed in A04A/O/S.
    for rel,h in (
        ("go_2026-08-05/human_indexed_models/MODEL_INVENTORY.json","03418cc71f0195d49f89f3bcbad09ca1e6416867c58afa0d0fe48304a4034d09"),
        ("ontology_interpretation/RO/2026-09-04/ro.obo","f2145eaecd649b395f5a908906f429fd6ac6dca2bc6816045b5c8a80da9a6b10"),
        ("ontology_interpretation/ECO/2026-07-10/eco.obo","c8b52a2da2b0f92224b50a1d8773fd7e3bca0977d920d4e31ef9fde960c84fd4"),
        ("ontology_interpretation/ChEBI/monthly_snapshot/chebi.obo.gz",CHEBI_HASH),
        ("uniprot_2026_04/targeted_go_enablers/IDENTITY_INVENTORY.json",UNI_HASH),
        ("uniprot_2026_04/targeted_go_participants/IDENTITY_INVENTORY.json",PART_HASH)):
        p=raw/rel
        need(sha(p)==h,"pinned reference mismatch "+rel)
        emit("REFERENCE_OR_INVENTORY_VERIFIED",rel,p.stat().st_size,h)
    # Database was rehashed in T v2. Avoid redundant table scans here.
    need(DB.stat().st_size==1451589632,"staging size")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    for sql in ("SELECT ordinal,rdf_id,class_uri FROM bp_record","SELECT record_ordinal,predicate_uri,target_local_id,literal_text FROM bp_property"):
        emit("SQLITE_QUERY_PLAN",sql,db.execute("EXPLAIN QUERY PLAN "+sql).fetchall())
    db.close()
    # Disease hierarchy is an annotation-only scope diagnostic, not a graph input.
    disease={a:b for a,b in pairs(core("R11")) if STABLE.fullmatch(a)}
    need(disease,"empty disease pathway inventory")
    parents=defaultdict(set);edges=set()
    for a,b in pairs(core("R09")):
        need(STABLE.fullmatch(a) and STABLE.fullmatch(b),"unexpected pathway relation row "+str((a,b)))
        if a.startswith("R-HSA-") and b.startswith("R-HSA-"):
            parents[b].add(a);edges.add((a,b))
    emit("HUMAN_PATHWAY_RELATIONS",len(edges),"DISEASE_PATHWAYS",len(disease))
    visiting=set()
    @lru_cache(None)
    def ancestry(p):
        need(p not in visiting,"pathway hierarchy cycle "+p)
        visiting.add(p)
        try:
            found=set()
            for q in parents[p]:
                found.add(q);found.update(ancestry(q))
            return frozenset(found)
        finally:visiting.remove(p)
    allp={a for a,b in edges}|{b for a,b in edges}
    for p in sorted(allp):ancestry(p)
    ops={r["stable_id"]:r for r in rows(TYPED/"biopax_operations.tsv")}
    need(len(ops)==16187,"event count")
    ep=defaultdict(set)
    for r in rows(core("R04")):
        if r["reaction_id"] in ops:ep[r["reaction_id"]].add(r["pathway_id"])
    need(len(ep)==15771,"R04 event count")
    def cls(ps,ds):
        if not ps:return "no_R04_mapping"
        n=ps&ds
        return "disease_only" if n==ps else "mixed" if n else "no_listed_disease_context"
    used={p for ps in ep.values() for p in ps}
    inherited={p for p in used if p in disease or ancestry(p)&disease.keys()}
    direct=set(disease)
    emit("R04_PATHWAY_CONTEXT_COUNTS","distinct",len(used),"direct_disease",len(used&direct),"disease_or_descendant",len(inherited),"not_in_human_relation_graph",len(used-allp))
    comparison=Counter();examples=defaultdict(list)
    for e in ops:
        ps=ep[e];a,b=cls(ps,direct),cls(ps,inherited)
        comparison[a,b]+=1
        if a!=b and len(examples[a,b])<5:examples[a,b].append((e,sorted(ps)))
    emit("DISEASE_CONTEXT_DIRECT_VS_ANCESTRY",sorted(comparison.items()))
    emit("DISEASE_CONTEXT_CHANGED_EXAMPLES",sorted(examples.items()))
    emit("DEFERRED_NATIVE_DUMP","R02","bytes verified; native semantic verification not performed; not a primary input")
    emit("DEFERRED_FULL_ARCHIVE","R14","not acquired")
    emit("P1A04U_SOURCE_SCOPE_PREFLIGHT_COMPLETE; NO_FILTERS_OR_MERGES")
if __name__=="__main__":main()
