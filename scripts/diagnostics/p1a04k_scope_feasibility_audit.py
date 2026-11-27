#!/usr/bin/env python3
"""Read-only human-scope feasibility audit; not a biological inclusion filter."""
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time

ROOT=Path(__file__).resolve().parents[2]
TYPED=ROOT/"results/P1A04F_TYPED_SOURCE_LAYER_V1"
STATE=ROOT/"results/P1A04H_STATE_CATALOGUE_V1"
START=time.monotonic()
# H: a human-compatible realization exists; F: every recorded realization
# requires a nonhuman component; U: unresolved. Chemicals are taxon-neutral.
# These are compositional feasibility categories, not activity or expression.
def all_of(xs):
    xs=list(xs)
    return "F" if "F" in xs else "U" if not xs or "U" in xs else "H"
def any_of(xs):
    xs=list(xs)
    return "H" if "H" in xs else "U" if not xs or "U" in xs else "F"
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
    need(m["schema_version"]==version,"manifest schema")
    for name,r in m["files"].items():
        need(Path(name).name==name,"unsafe manifest name")
        p=directory/name
        need(p.stat().st_size==r["bytes"] and sha(p)==r["sha256"],"manifest mismatch "+name)
    emit("DERIVED_MANIFEST_REHASH_PASS",version,len(m["files"]))
def self_test():
    need(all_of(["H","F"])=="F","AND foreign")
    need(any_of(["H","F"])=="H","OR human alternative")
    need(all_of(["H","U"])=="U" and any_of(["F","U"])=="U","unknown")
    need(all_of([])==any_of([])=="U","empty")
    emit("SELF_TEST_PASS")
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--self-test",action="store_true")
    if ap.parse_args().self_test: self_test();return
    check_manifest(TYPED,"P1A04F_v1")
    check_manifest(STATE,"P1A04H_v1")
    states={r["rdf_id"]:r for r in rows(STATE/"biopax_states.tsv")}
    need(len(states)==61418,"state count")
    components=defaultdict(list);members=defaultdict(list)
    for r in rows(STATE/"biopax_composition.tsv"):
        target=components if r["relation"]=="component" else members if r["relation"]=="memberPhysicalEntity" else None
        need(target is not None,"unknown composition relation")
        target[r["owner_id"]].append(r["child_id"])
    need(sum(map(len,components.values()))==40083 and sum(map(len,members.values()))==36113,"composition inventory")
    def leaf(r):
        t=set(json.loads(r["taxon_union_json"]))
        u=r["unresolved_taxon"]=="1"
        if r["class"]=="SmallMolecule":return "H"
        if u or not t or ("9606" in t and len(t)>1):return "U"
        return "H" if t=={"9606"} else "F"
    @lru_cache(None)
    def feasibility(e):
        r=states[e]
        c,m=components[e],members[e]
        if c and m:return "U" # unreviewed hybrid semantics
        if c:return all_of(feasibility(x) for x in c)
        if m:return any_of(feasibility(x) for x in m)
        return leaf(r)
    sc=Counter(); comparison=Counter()
    for e,r in states.items():
        f=feasibility(e);sc[(r["class"],f)]+=1
        t=set(json.loads(r["taxon_union_json"]))
        if "9606" in t and len(t)>1:
            comparison[(r["class"],f)]+=1
    emit("ENTITY_COMPOSITION_FEASIBILITY",sorted(sc.items()))
    emit("MIXED_TAXON_UNION_FEASIBILITY",sorted(comparison.items()))
    ops={r["rdf_id"]:r for r in rows(TYPED/"biopax_operations.tsv")}
    need(len(ops)==16187,"event count")
    control={r["rdf_id"]:r for r in rows(TYPED/"biopax_control_records.tsv")}
    controller=defaultdict(set);controlled=defaultdict(set)
    for r in rows(TYPED/"biopax_incidence.tsv"):
        if r["layer"]=="controller_of_record":controller[r["target_id"]].add(r["source_id"])
        elif r["layer"]=="controlled_by_record":controlled[r["source_id"]].add(r["target_id"])
    cats=defaultdict(list); reg=defaultdict(list); targetclasses=Counter()
    for c,r in control.items():
        es=controller[c]
        # Multiple controllers in a single record: do not infer AND/OR.
        fs={feasibility(e) if e in states else "U" for e in es}
        f=next(iter(fs)) if len(fs)==1 and len(es)==1 else "U"
        for target in controlled[c]:
            targetclasses[(r["class"],"event" if target in ops else "control" if target in control else "other")]+=1
            if target in ops:
                (cats if r["class"]=="Catalysis" else reg)[target].append(f)
    emit("CONTROL_TARGET_CLASSES",sorted(targetclasses.items()))
    parts=defaultdict(lambda:defaultdict(set))
    for r in rows(STATE/"biopax_event_participants.tsv"):
        parts[r["event_rdf_id"]][r["role"]].add(r["entity_rdf_id"])
    def profile(es):
        if not es:return "none_recorded"
        fs={feasibility(e) if e in states else "U" for e in es}
        return "+".join(sorted(fs))
    def ctrl_profile(xs):
        if not xs:return "none_recorded"
        return "+".join(sorted(set(xs)))
    ec=Counter(); examples=defaultdict(list); detailed=Counter()
    for e,r in ops.items():
        roles=parts[e]
        inp=profile(roles["left"]|roles["template"])
        out=profile(roles["right"]|roles["product"])
        ca=ctrl_profile(cats[e]); re=ctrl_profile(reg[e])
        detailed[(r["class"],ca,re,inp,out)]+=1
        ec[(r["class"],ca)]+=1
        k=(ca,inp,out)
        if len(examples[k])<2:
            examples[k].append((r["stable_id"],r["name"]))
    emit("EVENT_CATALYTIC_FEASIBILITY_PROFILES",sorted(ec.items()))
    emit("EVENT_ROLE_FEASIBILITY_PROFILES",sorted(detailed.items()))
    for k,v in sorted(examples.items()):
        if "F" in "+".join(k) or "U" in "+".join(k):emit("EVENT_SCOPE_REVIEW_EXAMPLE",k,v)
    # Exact UniProt identity observations from reference-bearing BP protein records.
    # No HGNC inference of species, no isoform stripping, no model-taxon imputation.
    acc_tax=defaultdict(set);acc_unknown=set()
    for e,r in states.items():
        if r["class"]!="Protein" or not json.loads(r["reference_ids_json"]):continue
        for _,xs in json.loads(r["reference_xrefs_json"]):
            for ns,acc in xs:
                if ns not in ("UniProt","UniProt Isoform"):continue
                acc_tax[acc].update(json.loads(r["taxon_union_json"]))
                if r["unresolved_taxon"]=="1":acc_unknown.add(acc)
    def acc_status(acc):
        t=acc_tax.get(acc,set())
        if not t or acc in acc_unknown:return "unresolved"
        if t=={"9606"}:return "observed_human_only"
        if "9606" in t:return "conflicting_or_mixed_taxa"
        return "observed_nonhuman_only"
    def terms(obj):
        if isinstance(obj,dict):
            t=obj.get("term")
            if isinstance(t,str) and t.startswith("UniProtKB:"):yield t.split(":",1)[1]
            yield from terms(obj.get("has_part",[]))
        elif isinstance(obj,list):
            for v in obj:yield from terms(v)
    gc=Counter();ga=Counter();unresolved=set()
    for r in rows(STATE/"gocam_enablers.tsv"):
        en=json.loads(r["enabler_json"])
        acc=set(terms(en));statuses={acc_status(x) for x in acc}
        source="Reactome_namespace" if r["model_id"].startswith("R-HSA-") else "other_namespace"
        gc[(source,r["enabler_type"],"+".join(sorted(statuses)) if statuses else "no_UniProt_accession")]+=1
        for x in acc:
            s=acc_status(x);ga[(source,s)]+=1
            if s=="unresolved":unresolved.add(x)
    emit("GO_ENABLER_ACCESSION_TAXON_PROFILES",sorted(gc.items()))
    emit("GO_ENABLER_ACCESSION_OCCURRENCES",sorted(ga.items()))
    emit("GO_UNRESOLVED_DISTINCT_UNIPROT",len(unresolved),"EXAMPLES",sorted(unresolved)[:30])
    emit("P1A04K_SCOPE_FEASIBILITY_AUDIT_COMPLETE; NO_SCOPE_EXCLUSIONS")
if __name__=="__main__":main()
