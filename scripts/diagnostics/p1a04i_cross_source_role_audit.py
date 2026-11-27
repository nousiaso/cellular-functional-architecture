#!/usr/bin/env python3
"""P1A.04I: bounded, read-only source overlap and role-specific identity audit."""
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import time

ROOT=Path(__file__).resolve().parents[2]
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
TYPED=ROOT/"results/P1A04F_TYPED_SOURCE_LAYER_V1"
STATE=ROOT/"results/P1A04H_STATE_CATALOGUE_V1"
START=time.monotonic()
def emit(*x): print(*x,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def rows(p):
    with p.open(newline="",encoding="utf-8") as f:
        yield from csv.DictReader(f,delimiter="\t")
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def check_manifest(directory,version):
    m=json.loads((directory/"manifest.json").read_text())
    need(m["schema_version"]==version,"manifest schema mismatch")
    for name,rec in m["files"].items():
        need(Path(name).name==name,"invalid manifest file name")
        p=directory/name
        need(p.stat().st_size==rec["bytes"] and sha(p)==rec["sha256"],"derived file mismatch "+name)
    emit("DERIVED_MANIFEST_REHASH_PASS",version,len(m["files"]))
def compare(a,b):
    if not a and not b:return "both_empty"
    if not a:return "go_only"
    if not b:return "biopax_only"
    if a==b:return "equal_nonempty"
    if a<b:return "go_subset"
    if b<a:return "biopax_subset"
    if a&b:return "partial_overlap"
    return "disjoint"
def pmid(s):
    if not isinstance(s,str):return None
    m=re.fullmatch(r"PMID:\s*([0-9]+)",s.strip())
    return m.group(1) if m else None
def refs(obj):
    if isinstance(obj,dict):
        for k,v in obj.items():
            if k=="reference" and isinstance(v,str):yield v
            yield from refs(v)
    elif isinstance(obj,list):
        for v in obj:yield from refs(v)
def term_ids(obj):
    # Exact UniProt accession tokens only; never strip isoform suffixes.
    if isinstance(obj,dict):
        for k,v in obj.items():
            if k=="term" and isinstance(v,str) and v.startswith("UniProtKB:"):
                yield v.split(":",1)[1]
            # Restrict recursion to enabler composition, not evidence or part_of.
            if k=="has_part":yield from term_ids(v)
    elif isinstance(obj,list):
        for v in obj:yield from term_ids(v)
def main():
    check_manifest(TYPED,"P1A04F_v1")
    check_manifest(STATE,"P1A04H_v1")
    need(DB.stat().st_size==1451589632,"staged DB size")
    states={r["rdf_id"]:r for r in rows(STATE/"biopax_states.tsv")}
    need(len(states)==61418,"state count")
    children=defaultdict(list)
    for r in rows(STATE/"biopax_composition.tsv"):
        children[r["owner_id"]].append(r["child_id"])
    @lru_cache(None)
    def proteins(e):
        if e not in states:return frozenset()
        acc=set()
        for _,xs in json.loads(states[e]["reference_xrefs_json"]):
            for ns,id_ in xs:
                if ns in ("UniProt","UniProt Isoform"):acc.add(id_)
        for child in children[e]:acc.update(proteins(child))
        return frozenset(acc)
    ops={r["stable_id"]:r["rdf_id"] for r in rows(TYPED/"biopax_operations.tsv")}
    need(len(ops)==16187,"event count")
    cat=set(); controls=defaultdict(set); controlled=defaultdict(set)
    for r in rows(TYPED/"biopax_control_records.tsv"):
        if r["class"]=="Catalysis":cat.add(r["rdf_id"])
    for r in rows(TYPED/"biopax_incidence.tsv"):
        if r["layer"]=="controller_of_record":controls[r["target_id"]].add(r["source_id"])
        elif r["layer"]=="controlled_by_record":controlled[r["source_id"]].add(r["target_id"])
    catalysts=defaultdict(set)
    for c in cat:
        for event in controlled[c]:catalysts[event].update(controls[c])
    def tax_summary(e):
        if e not in states:return "outside_state_catalogue"
        r=states[e]; t=set(json.loads(r["taxon_union_json"]))
        u=r["unresolved_taxon"]=="1"
        v="human_only_known" if t=={"9606"} else "mixed_human_nonhuman" if "9606" in t else "nonhuman_known" if t else "no_known_taxon"
        return v+("_unresolved" if u else "")
    # No human-catalyst assertion for an OR set solely from its taxon union.
    tc=Counter()
    for event in ops.values():
        ts={tax_summary(e) for e in catalysts[event]}
        tc[tuple(sorted(ts))]+=1
    emit("EVENT_CATALYST_TAXON_PROFILES",sorted(tc.items()))
    parts=defaultdict(lambda:defaultdict(set))
    for r in rows(STATE/"biopax_event_participants.tsv"):
        parts[r["event_rdf_id"]][r["role"]].add(r["entity_rdf_id"])
    # Explicit role usage of fingerprint groups; no group merges.
    use=defaultdict(set)
    for event,roles in parts.items():
        for role,es in roles.items():
            for e in es:use[e].add((event,role))
    for event,es in catalysts.items():
        for e in es:use[e].add((event,"catalyst"))
    groups=defaultdict(list)
    for r in rows(STATE/"state_fingerprint_candidates.tsv"):groups[r["fingerprint"]].append(r["rdf_id"])
    need(len(groups)==440 and sum(map(len,groups.values()))==947,"fingerprint group inventory")
    fc=Counter()
    for es in groups.values():
        f={x for e in es for x in json.loads(states[e]["flags_json"])}
        fc["groups_with_any_uncertainty" if f else "groups_without_recorded_flags"]+=1
        fc["groups_with_different_event_role_usage" if len({frozenset(use[e]) for e in es})>1 else "groups_with_identical_event_role_usage"]+=1
        fc["groups_with_any_event_role_usage" if any(use[e] for e in es) else "groups_without_event_role_usage"]+=1
        fc["groups_with_different_display_names" if len({states[e]["name_json"] for e in es})>1 else "groups_with_identical_display_names"]+=1
    emit("FINGERPRINT_GROUP_REVIEW",dict(fc))
    emit("STOICHIOMETRY_REVIEW_ROWS",list(rows(STATE/"reaction_stoichiometry_review.tsv")))
    # Reaction-PMID export is a provenance comparison, not the complete evidence set.
    r06=defaultdict(set)
    receipt_rows=json.loads((ROOT/"data/raw/P1A03_CORE_RECEIPTS.json").read_text())
    matches=[r for r in receipt_rows if r.get("id")=="R06"]
    need(len(matches)==1,"R06 receipt cardinality")
    receipt=matches[0]
    path=(ROOT/receipt["path"]).resolve()
    need(path.is_relative_to((ROOT/"data/raw").resolve()),"R06 path outside raw directory")
    need(path.is_file(),"R06 receipt path missing")
    need(path.stat().st_size==receipt["bytes"]==1136864,"R06 size mismatch")
    need(sha(path)==receipt["sha256"]=="3b2169fc94cbe51b3eea84e28c11c37512ebd834c7bca97f18719f0027a5abf8","R06 checksum mismatch")
    emit("R06_SOURCE_INTEGRITY_PASS",str(path.relative_to(ROOT)))
    with path.open(newline="",encoding="utf-8") as f:
        for row in csv.reader(f,delimiter="\t"):
            if len(row)==2 and row[0] in ops and row[1].isdigit():r06[row[0]].add(row[1])
    emit("REACTION_PMID_EXPORT_EVENTS",len(r06))
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","schema mismatch")
    models={m:(t,json.loads(ts)) for m,t,ts in db.execute("SELECT model_id,primary_taxon,indexed_taxa_json FROM go_model")}
    counts=Counter(); example=defaultdict(list); seen=set(); all_refs=Counter()
    for mid,sid,payload in db.execute("SELECT model_id,source_id,payload_json FROM go_record WHERE field='activities'"):
        obj=json.loads(payload); eb=obj.get("enabled_by") or {}
        rs=set(refs(obj)); ps={x for s in rs if (x:=pmid(s))}
        for s in rs:
            all_refs["PMID_exact_or_whitespace"]+=bool(pmid(s))
            all_refs["Reactome_reference"]+=s.strip().startswith("Reactome:")
            all_refs["other_reference"]+=not pmid(s) and not s.strip().startswith("Reactome:")
        raw=(sid or "").removeprefix("gomodel:")
        cited={s.strip().split(":",1)[1] for s in rs if s.strip().startswith("Reactome:")} & ops.keys()
        if raw not in ops:
            counts["nonbare_or_unmatched"]+=1
            counts["nonbare_or_unmatched_with_event_citation"]+=bool(cited)
            continue
        counts["bare_match_occurrences"]+=1; seen.add(raw)
        event=ops[raw]
        g=set(term_ids(eb))
        b={p for e in catalysts[event] for p in proteins(e)}
        tag=compare(g,b)
        counts["enabler_vs_catalyst:"+tag]+=1
        # These are flattened accession inventories for audit only, not AND/OR equivalence.
        if len(example[tag])<3:
            example[tag].append({"model":mid,"activity":sid,"enabler_type":eb.get("type"),
                 "go_accessions":sorted(g)[:20],"bp_catalyst_accessions":sorted(b)[:20],
                 "bp_catalyst_taxon_profiles":sorted({tax_summary(e) for e in catalysts[event]}),
                 "model_taxon":models[mid][0],"indexed_taxa":models[mid][1]})
        counts["bare_match_with_same_event_citation"]+=raw in cited
        counts["bare_match_with_direct_PMID"]+=bool(ps)
        counts["bare_match_with_export_PMID_overlap"]+=bool(ps & r06[raw])
    need(counts["bare_match_occurrences"]==682 and counts["nonbare_or_unmatched"]==6703,"GO activity crosswalk count changed")
    emit("GO_BARE_MATCH_DISTINCT_EVENTS",len(seen))
    emit("GO_EVENT_OVERLAP_AUDIT",dict(counts))
    emit("GO_ACTIVITY_REFERENCE_INVENTORY",dict(all_refs))
    for k,v in sorted(example.items()):emit("GO_ENABLER_COMPARISON_EXAMPLES",k,v)
    db.close()
    emit("P1A04I_CROSS_SOURCE_ROLE_AUDIT_COMPLETE; NO_MERGES_OR_SCOPE_EXCLUSIONS")
if __name__=="__main__":main()
