#!/usr/bin/env python3
"""Read-only certification of acquired identity snapshot and crosswalk coverage."""
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"acquisition"))
from p1a04n_acquire_targeted_uniprot import (
    ROOT,STATE,OUT,RELEASE,RELEASE_DATE,tokens,parent,verify,identity,
    digest,need,rows,emit)
from p1a04k_scope_feasibility_audit import sha

EXPECTED="e9a48281b9b96e481902f09a9346a702aca19355623d81d2d97deecbd6d2b05d"
def main():
    ts=tokens()
    parents=sorted({parent(t) for t in ts})
    plan=json.loads((OUT/"PLAN.json").read_text())
    need(plan=={"schema_version":"P1A04N_v1","release":RELEASE,
                "release_date":RELEASE_DATE,"tokens":ts,"parents":parents},"plan drift")
    path=OUT/"IDENTITY_INVENTORY.json"
    need(sha(path)==EXPECTED,"inventory checksum")
    inv=json.loads(path.read_text())
    need(inv["schema_version"]=="P1A04N_v1_1" and inv["release"]==RELEASE,"inventory version")
    need(set(inv["parent_receipts"])==set(parents) and len(parents)==490,"parent inventory")
    need(not list(OUT.rglob("*.partial")),"partial artifacts")
    need({x.name for x in (OUT/"records").iterdir()}=={x+".json" for x in parents},"receipt directory inventory")
    objs={}; nbytes=0
    for acc in parents:
        p=OUT/"records"/(acc+".json")
        b=p.read_bytes();r=json.loads(b)
        objs[acc]=verify(acc,r)
        expected={"status":r["status"],"bytes":r["bytes"],"payload_sha256":r["sha256"],"receipt_sha256":digest(b)}
        need(inv["parent_receipts"][acc]==expected,"receipt mismatch "+acc)
        need(r["status"]==200 and objs[acc] is not None,"unexpected unsuccessful lookup "+acc)
        nbytes+=r["bytes"]
    regenerated=[identity(t,objs[parent(t)]) for t in ts]
    need(regenerated==inv["token_records"],"identity reconstruction mismatch")
    by_token={r["token"]:r for r in regenerated}
    need(len(by_token)==520,"token inventory")
    counts=Counter((r["identity"],r["substate"]) for r in regenerated)
    need(counts=={("primary","parent_token"):395,("primary","feature_exact"):73,
                  ("primary","isoform_exact"):51,("inactive_entry","unresolved"):1},"identity counts")
    emit("IDENTITY_INVENTORY_REHASH_AND_RECONSTRUCTION_PASS",len(parents),len(ts),"PAYLOAD_BYTES",nbytes)
    emit("TOKEN_TAXONOMY_BY_SUBSTATE",sorted(Counter((r["taxon_id"],r["substate"]) for r in regenerated).items(),key=str))
    # HGNC exact accession and UniProt HGNC xref are two separate observations.
    receipt=json.loads((ROOT/"data/raw/P1A03_CORE_RECEIPTS.json").read_text())
    matches=[r for r in receipt if r.get("id")=="H01"]
    need(len(matches)==1,"HGNC receipt")
    r=matches[0];p=(ROOT/r["path"]).resolve()
    need(p.is_relative_to((ROOT/"data/raw").resolve()),"HGNC path")
    pinned=[x for x in rows(ROOT/"manifests/P1A03_CORE_REPORTED_SHA256.tsv") if x["id"]=="H01"]
    need(len(pinned)==1 and sha(p)==r["sha256"]==pinned[0]["sha256"] and p.stat().st_size==r["bytes"]==int(pinned[0]["bytes"]),"HGNC integrity")
    hgnc=defaultdict(set);approved={}
    for r in rows(p):
        if r["status"]!="Approved":continue
        approved[r["hgnc_id"]]=r["symbol"]
        for acc in r["uniprot_ids"].split("|"):
            if acc.strip():hgnc[acc.strip()].add(r["hgnc_id"])
    hc=Counter();ex=defaultdict(list)
    for r in regenerated:
        if r["taxon_id"]!=9606:continue
        a=hgnc[r["primary"]];b=set(r["HGNC_xrefs"])
        tag="equal_unique" if a==b and len(a)==1 else "equal_multiple" if a==b and len(a)>1 else "both_empty" if not a and not b else "HGNC_only" if a and not b else "UniProt_only" if b and not a else "overlap_not_equal" if a&b else "disjoint"
        hc[tag]+=1
        if tag!="equal_unique" and len(ex[tag])<10:ex[tag].append((r["token"],sorted(a),sorted(b),r["gene_names"]))
        if b-approved.keys():hc["tokens_with_unapproved_or_absent_UniProt_HGNC_xref"]+=1
    emit("TARGETED_HUMAN_PARENT_HGNC_COMPARISON",dict(hc))
    emit("HGNC_REVIEW_EXAMPLES",dict(ex))
    # Preserve provenance of taxonomy: previously observed in BP versus targeted UniProt.
    bt=defaultdict(set);unknown=set();chem=Counter();ref_ids=defaultdict(set)
    for r in rows(STATE/"biopax_states.tsv"):
        cl=r["class"]
        for ref,xs in json.loads(r["reference_xrefs_json"]):
            if cl=="SmallMolecule":
                for ns,acc in xs:ref_ids[ns].add((ref,acc))
            if cl=="Protein":
                for ns,acc in xs:
                    if ns in ("UniProt","UniProt Isoform"):
                        bt[acc].update(json.loads(r["taxon_union_json"]))
                        if r["unresolved_taxon"]=="1":unknown.add(acc)
        if cl=="SmallMolecule":
            namespaces={ns for _,xs in json.loads(r["reference_xrefs_json"]) for ns,_ in xs}
            chem[tuple(sorted(namespaces))]+=1
    def terms(obj):
        if isinstance(obj,dict):
            t=obj.get("term")
            if isinstance(t,str):yield t
            yield from terms(obj.get("has_part",[]))
        elif isinstance(obj,list):
            for x in obj:yield from terms(x)
    def status(acc):
        if acc in by_token:
            t=by_token[acc]["taxon_id"]
            return "targeted_human" if t==9606 else "targeted_unresolved" if t is None else "targeted_nonhuman"
        t=bt.get(acc,set())
        return "BP_unresolved" if not t or acc in unknown else "BP_human" if t=={"9606"} else "BP_mixed" if "9606" in t else "BP_nonhuman"
    ac=Counter();occ=Counter();other=Counter();other_ids=defaultdict(set)
    for r in rows(STATE/"gocam_enablers.tsv"):
        obj=json.loads(r["enabler_json"]); tt=set(terms(obj))
        acc={x.split(":",1)[1] for x in tt if x.startswith("UniProtKB:")}
        profile=tuple(sorted({status(x) for x in acc}))
        ac[(r["enabler_type"],profile)]+=1
        for x in acc:occ[status(x)]+=1
        for x in tt:
            if not x.startswith("UniProtKB:"):
                ns=x.split(":",1)[0];other[ns]+=1;other_ids[ns].add(x)
        if "P03208" in acc:emit("INACTIVE_TOKEN_ACTIVITY",r["model_id"],r["activity_id"],r["enabler_type"],r["enabler_term"])
    emit("GO_ENABLER_UPDATED_TAXONOMY_PROFILES",sorted(ac.items()))
    emit("GO_ENABLER_ACCESSION_OCCURRENCE_COVERAGE",dict(occ))
    emit("GO_NON_UNIPROT_ENABLER_TERM_OCCURRENCES",dict(other))
    emit("GO_NON_UNIPROT_ENABLER_DISTINCT_TERMS",[(k,len(v),sorted(v)[:10]) for k,v in sorted(other_ids.items())])
    emit("BP_SMALL_MOLECULE_REFERENCE_NAMESPACE_PROFILES",sorted(chem.items()))
    emit("BP_SMALL_MOLECULE_DISTINCT_REFERENCE_ACCESSION_PAIRS",sorted((k,len(v)) for k,v in ref_ids.items()))
    emit("P1A04O_IDENTITY_COVERAGE_AUDIT_COMPLETE; NO_MERGES_OR_EXCLUSIONS")
if __name__=="__main__":main()
