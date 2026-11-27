#!/usr/bin/env python3
"""Read-only targeted scope exceptions, disease-context coverage and identity intake."""
from collections import Counter, defaultdict
from functools import lru_cache
import csv
import json
from pathlib import Path
import re
import time
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from p1a04k_scope_feasibility_audit import (ROOT,TYPED,STATE,rows,sha,need,check_manifest,all_of,any_of)

START=time.monotonic()
def emit(*x): print(*x,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def core_sources():
    receipts=json.loads((ROOT/"data/raw/P1A03_CORE_RECEIPTS.json").read_text())
    pinned={r["id"]:r for r in rows(ROOT/"manifests/P1A03_CORE_REPORTED_SHA256.tsv")}
    result={}
    for ident in ("R04","R11","H01"):
        matches=[r for r in receipts if r.get("id")==ident]
        need(len(matches)==1,"receipt cardinality "+ident)
        r=matches[0]; p=(ROOT/r["path"]).resolve()
        need(p.is_relative_to((ROOT/"data/raw").resolve()),"source outside raw directory")
        need(p.stat().st_size==r["bytes"]==int(pinned[ident]["bytes"]),"size "+ident)
        need(sha(p)==r["sha256"]==pinned[ident]["sha256"],"hash "+ident)
        result[ident]=p
        emit("CORE_SOURCE_REHASH_PASS",ident,p.name)
    return result
def main():
    sources=core_sources()
    check_manifest(TYPED,"P1A04F_v1")
    check_manifest(STATE,"P1A04H_v1")
    states={r["rdf_id"]:r for r in rows(STATE/"biopax_states.tsv")}
    need(len(states)==61418,"state count")
    comp=defaultdict(list); members=defaultdict(list)
    for r in rows(STATE/"biopax_composition.tsv"):
        dest=comp if r["relation"]=="component" else members if r["relation"]=="memberPhysicalEntity" else None
        need(dest is not None,"unknown composition relation")
        dest[r["owner_id"]].append(r["child_id"])
    def taxa(e):return set(json.loads(states[e]["taxon_union_json"]))
    @lru_cache(None)
    def feasibility(e):
        r=states[e];c,m=comp[e],members[e]
        if c and m:return "U"
        if c:return all_of(feasibility(x) for x in c)
        if m:return any_of(feasibility(x) for x in m)
        t=taxa(e)
        if r["class"]=="SmallMolecule":return "H"
        if r["unresolved_taxon"]=="1" or not t or ("9606" in t and len(t)>1):return "U"
        return "H" if t=={"9606"} else "F"
    @lru_cache(None)
    def leaves(e):
        children=comp[e]+members[e]
        return frozenset().union(*(leaves(c) for c in children)) if children else frozenset([e])
    def profile(e):
        p=set()
        for x in leaves(e):
            r=states[x];t=taxa(x)
            if r["class"]=="SmallMolecule":p.add("chemical")
            elif r["class"]=="Protein":
                if "9606" in t:p.add("human_protein")
                if t-{"9606"}:p.add("nonhuman_protein")
                if not t or r["unresolved_taxon"]=="1":p.add("unknown_protein")
            elif r["class"] in ("Dna","Rna"):
                if "9606" in t:p.add("human_nucleic_acid")
                if t-{"9606"}:p.add("nonhuman_nucleic_acid")
                if not t or r["unresolved_taxon"]=="1":p.add("unknown_nucleic_acid")
            else:p.add("other_leaf")
        return tuple(sorted(p))
    ops={r["rdf_id"]:r for r in rows(TYPED/"biopax_operations.tsv")}
    need(len(ops)==16187,"event count")
    ids={r["stable_id"]:e for e,r in ops.items()}
    control={r["rdf_id"]:r for r in rows(TYPED/"biopax_control_records.tsv")}
    controllers=defaultdict(set);targets=defaultdict(set)
    for r in rows(TYPED/"biopax_incidence.tsv"):
        if r["layer"]=="controller_of_record":controllers[r["target_id"]].add(r["source_id"])
        elif r["layer"]=="controlled_by_record":targets[r["source_id"]].add(r["target_id"])
    cat=defaultdict(set)
    for c,r in control.items():
        if r["class"]=="Catalysis":
            for e in targets[c]:cat[e].update(controllers[c])
    foreign=Counter(); examples=defaultdict(list); foreign_events=set()
    for e,es in cat.items():
        for x in es:
            if x not in states or feasibility(x)!="F":continue
            foreign_events.add(e)
            k=profile(x);foreign[k]+=1
            if len(examples[k])<3:examples[k].append((ops[e]["stable_id"],ops[e]["name"],x))
    emit("FOREIGN_REQUIRED_CATALYST_EVENTS",len(foreign_events))
    emit("FOREIGN_REQUIRED_CATALYST_LEAF_PROFILES",sorted(foreign.items()))
    for k,v in sorted(examples.items()):emit("FOREIGN_PROFILE_EXAMPLES",k,v)
    # Explicit component vs member paths; bounded examples retain all leaf IDs.
    for stable in ("R-HSA-8983680","R-HSA-8985104","R-HSA-8943080","R-HSA-167191","R-HSA-9661820","R-HSA-5210935"):
        if stable not in ids:emit("TARGET_EVENT_MISSING",stable);continue
        e=ids[stable]
        emit("TARGET_EVENT",stable,ops[e]["name"])
        for x in sorted(cat[e]):
            emit("TARGET_CATALYST",stable,x,states[x]["name_json"],feasibility(x),profile(x))
            visited=set()
            def trace(n):
                if n in visited:return
                visited.add(n);r=states[n]
                emit("COMPONENT_TRACE",stable,n,r["class"],r["name_json"],r["taxon_union_json"],r["reference_xrefs_json"],"components",comp[n],"members",members[n])
                for child in comp[n]+members[n]:trace(child)
            trace(x)
    # R04 only records events present in that export. Disease membership is an
    # annotation flag, not proof of disease-only function or failed reaction.
    disease={}
    with sources["R11"].open(newline="",encoding="utf-8") as f:
        for row in csv.reader(f,delimiter="\t"):
            if len(row)==2:disease[row[0]]=row[1]
    ep=defaultdict(set)
    for r in rows(sources["R04"]):
        if r["reaction_id"] in ids:ep[r["reaction_id"]].add(r["pathway_id"])
    need(len(ep)==15771,"R04 event coverage changed")
    def context(s):
        ps=ep.get(s,set())
        if not ps:return "no_R04_mapping"
        ds=ps&disease.keys()
        return "listed_disease_only" if ds==ps else "mixed_listed_and_unlisted" if ds else "no_listed_disease_pathway"
    emit("R04_DISEASE_CONTEXT_ALL_EVENTS",dict(Counter(context(s) for s in ids)))
    emit("R04_DISEASE_CONTEXT_FOREIGN_CATALYST",dict(Counter(context(ops[e]["stable_id"]) for e in foreign_events)))
    # Missing product is not equivalent to negation. Report explicit name cues
    # for human review, without automatically assigning failed-event status.
    out=set()
    for r in rows(STATE/"biopax_event_participants.tsv"):
        if r["role"] in ("right","product"):out.add(r["event_rdf_id"])
    cues=re.compile(r"\b(?:does not|do not|fails? to|unable to|defective|inactive|cannot|can not)\b",re.I)
    nc=Counter();ex=defaultdict(list)
    for e,r in ops.items():
        k=(r["class"],e in out,bool(cues.search(r["name"])),context(r["stable_id"]))
        nc[k]+=1
        if (not k[1] or k[2]) and len(ex[k])<2:ex[k].append((r["stable_id"],r["name"]))
    emit("EVENT_OUTPUT_AND_NAME_CUE_CONTEXT",sorted(nc.items()))
    for k,v in sorted(ex.items()):emit("OUTPUT_NAME_CUE_EXAMPLES",k,v)
    # GO unresolved accession inventory based on the exact prior audit rule.
    bt=defaultdict(set);unknown=set()
    for e,r in states.items():
        if r["class"]!="Protein" or not json.loads(r["reference_ids_json"]):continue
        for _,xs in json.loads(r["reference_xrefs_json"]):
            for ns,acc in xs:
                if ns in ("UniProt","UniProt Isoform"):
                    bt[acc].update(taxa(e))
                    if r["unresolved_taxon"]=="1":unknown.add(acc)
    def terms(obj):
        if isinstance(obj,dict):
            t=obj.get("term")
            if isinstance(t,str) and t.startswith("UniProtKB:"):yield t.split(":",1)[1]
            yield from terms(obj.get("has_part",[]))
        elif isinstance(obj,list):
            for v in obj:yield from terms(v)
    go=set()
    for r in rows(STATE/"gocam_enablers.tsv"):go.update(terms(json.loads(r["enabler_json"])))
    unresolved={x for x in go if not bt.get(x) or x in unknown}
    need(len(unresolved)==520,"GO unresolved accession inventory changed")
    hgnc=defaultdict(set)
    for r in rows(sources["H01"]):
        if r["status"]=="Approved":
            for x in r["uniprot_ids"].split("|"):
                if x.strip():hgnc[x.strip()].add(r["hgnc_id"])
    stats=Counter();needs_lookup=[]
    for x in sorted(unresolved):
        exact=hgnc.get(x,set())
        # Base candidates for isoforms or processed chains are not equivalent IDs.
        base=re.sub(r"-(?:[0-9]+|PRO_[0-9]+)$","",x)
        candidate=hgnc.get(base,set()) if base!=x else set()
        k="exact_unique" if len(exact)==1 else "exact_ambiguous" if exact else "base_candidate_unique" if len(candidate)==1 else "base_candidate_ambiguous" if candidate else "no_HGNC_candidate"
        stats[k]+=1
        if k!="exact_unique":needs_lookup.append(x)
    emit("GO_UNRESOLVED_HGNC_CANDIDATE_CLASSES",dict(stats))
    emit("GO_UNRESOLVED_ALL_EXACT_ACCESSIONS",sorted(unresolved))
    emit("GO_UNRESOLVED_NONEXACT_HGNC_ACCESSIONS",needs_lookup)
    emit("P1A04L_SCOPE_EXCEPTION_AUDIT_COMPLETE; NO_EXCLUSIONS")
if __name__=="__main__":main()
