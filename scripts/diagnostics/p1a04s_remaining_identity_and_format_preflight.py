#!/usr/bin/env python3
"""Bounded remaining identity review and SBML format reconnaissance; read only."""
from collections import Counter,defaultdict
import json
from pathlib import Path
import re
import sqlite3
import sys
import tarfile
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"acquisition"))
import p1a04r_remaining_participant_identities as r
import p1a04n_acquire_targeted_uniprot as u
from p1a04p_chemical_reference_audit import ROOT,STATE,DB,rows,sha,need
from p1a04e_ontology_references import parse_obo

PART_HASH="dabdab6f05b04f84a3173b147996564069587d7af2e9c6a8a07d6e1f900e49fb"
def emit(*x):print(*x,flush=True)
def core(id):
    receipts=json.loads((ROOT/"data/raw/P1A03_CORE_RECEIPTS.json").read_text())
    rr=[x for x in receipts if x["id"]==id]
    pins=[x for x in rows(ROOT/"manifests/P1A03_CORE_REPORTED_SHA256.tsv") if x["id"]==id]
    need(len(rr)==len(pins)==1,"core receipt "+id)
    rec=rr[0];p=(ROOT/rec["path"]).resolve()
    need(p.is_relative_to((ROOT/"data/raw").resolve()),"source path")
    need(p.stat().st_size==rec["bytes"]==int(pins[0]["bytes"]) and sha(p)==rec["sha256"]==pins[0]["sha256"],"source hash "+id)
    return p
def verify_participants():
    inv,gaps,_=r.inventory()
    p=r.OUT/"IDENTITY_INVENTORY.json"
    need(sha(p)==PART_HASH,"participant inventory hash")
    out=json.loads(p.read_text())
    need(out["schema_version"]=="P1A04R_v1" and out["release"]==u.RELEASE,"participant schema")
    ts=sorted(gaps);parents=sorted({u.parent(x) for x in ts})
    need(set(out["parent_receipts"])==set(parents),"parent inventory")
    need(json.loads((r.OUT/"PLAN.json").read_text())=={"schema_version":"P1A04R_v1","release":u.RELEASE,"source_inventory_sha256":r.UNI_HASH,"tokens":ts,"parents":parents},"plan mismatch")
    need(not list(r.OUT.rglob("*.partial")),"partial records")
    need({x.name for x in (r.OUT/"records").iterdir()}=={x+".json" for x in parents if x not in inv["parent_receipts"]},"participant receipt files")
    objs={};origins=Counter()
    for acc in parents:
        x=out["parent_receipts"][acc]
        origin="enabler_snapshot" if acc in inv["parent_receipts"] else "participant_snapshot"
        fp=(u.OUT if origin=="enabler_snapshot" else r.OUT)/"records"/(acc+".json")
        rec=json.loads(fp.read_text())
        need(x=={"origin":origin,"receipt_sha256":sha(fp),"payload_sha256":rec["sha256"],"status":rec["status"]},"receipt mismatch "+acc)
        if origin=="enabler_snapshot":need(sha(fp)==inv["parent_receipts"][acc]["receipt_sha256"],"enabler receipt mismatch")
        objs[acc]=u.verify(acc,rec)
        need(objs[acc] is not None,"lookup unresolved "+acc)
        origins[origin]+=1
    regen=[u.identity(t,objs[u.parent(t)]) for t in ts]
    need(regen==out["token_records"],"participant reconstruction mismatch")
    need(Counter((x["identity"],x["substate"]) for x in regen)=={("primary","feature_exact"):10,("primary","isoform_exact"):13,("primary","parent_token"):22},"participant identity counts")
    emit("PARTICIPANT_RECONSTRUCTION_PASS",len(ts),dict(origins))
    return regen
def hgnc_audit(regen):
    h=defaultdict(set)
    for x in rows(core("H01")):
        if x["status"]=="Approved":
            for acc in x["uniprot_ids"].split("|"):
                if acc.strip():h[acc.strip()].add(x["hgnc_id"])
    counts=Counter()
    for x in regen:
        if x["taxon_id"]!=9606:continue
        a=h[x["primary"]];b=set(x["HGNC_xrefs"])
        tag="equal_unique" if a==b and len(a)==1 else "review"
        counts[tag]+=1
        if tag=="review":emit("PARTICIPANT_HGNC_REVIEW",x["token"],sorted(a),sorted(b))
    emit("PARTICIPANT_HGNC_COMPARISON",dict(counts))
    # Audit remaining BP human protein reference mapping without flattening
    # isoforms or treating an ambiguous parent candidate as a certified state.
    refs=defaultdict(set);tax=defaultdict(set);unresolved=set()
    for x in rows(STATE/"biopax_states.tsv"):
        if x["class"]!="Protein":continue
        for ref,xs in json.loads(x["reference_xrefs_json"]):
            refs[ref].update(tuple(y) for y in xs)
            tax[ref].update(json.loads(x["taxon_union_json"]))
            if x["unresolved_taxon"]=="1":unresolved.add(ref)
    c=Counter();examples=defaultdict(list)
    for ref,xs in refs.items():
        if tax[ref]!={"9606"} or ref in unresolved:continue
        acc={v for ns,v in xs if ns in ("UniProt","UniProt Isoform")}
        exact=set().union(*(h[x] for x in acc))
        base=set().union(*(h[u.parent(x)] for x in acc))
        k="exact_unique" if len(exact)==1 else "exact_ambiguous" if exact else "parent_candidate_unique" if len(base)==1 else "parent_candidate_ambiguous" if base else "no_candidate"
        c[k]+=1
        if k!="exact_unique" and len(examples[k])<5:examples[k].append((ref,sorted(acc),sorted(exact),sorted(base)))
    emit("BP_HUMAN_REFERENCE_HGNC_CANDIDATES",dict(c))
    emit("BP_HUMAN_REFERENCE_HGNC_EXAMPLES",dict(examples))
def nonprotein():
    header,go=parse_obo(core("G03").read_bytes())
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    names=defaultdict(set);counts=Counter();complex_profiles=Counter();complex_ids=defaultdict(list)
    for mid,sid,payload in db.execute("SELECT model_id,source_id,payload_json FROM go_record WHERE field='activities'"):
        obj=json.loads(payload);en=obj["enabled_by"];t=en["term"];ns=t.split(":",1)[0]
        if ns=="UniProtKB":continue
        counts[(ns,"enabler")]+=1;names[ns].add(t)
        if ns in ("GO","ComplexPortal"):
            parts=en.get("has_part",[])
            k=(ns,"with_recorded_parts" if parts else "no_recorded_parts")
            complex_profiles[k]+=1
            complex_ids[t].append((mid,sid,len(parts)))
    for mid,sid,payload in db.execute("SELECT model_id,source_id,payload_json FROM go_record WHERE field='molecules'"):
        obj=json.loads(payload);t=obj["term"];ns=t.split(":",1)[0]
        if ns in ("UniProtKB","CHEBI"):continue
        counts[(ns,"molecule")]+=1;names[ns].add(t)
    db.close()
    emit("NON_UNIPROT_NON_CHEMICAL_AND_CHEBI_ENABLER_OCCURRENCES",sorted(counts.items()))
    emit("REMAINING_NON_UNIPROT_TERMS",[(k,len(v),sorted(v)) for k,v in sorted(names.items()) if k!="GO"])
    emit("GO_COMPLEX_ENABLER_COMPOSITION",sorted(complex_profiles.items()))
    emit("GO_TERM_IDENTITY",[(t,go.get(t,{}).get("name"),go.get(t,{}).get("namespace"),go.get(t,{}).get("is_obsolete")) for t in sorted(names["GO"])])
    emit("COMPLEXPORTAL_ENABLER_CONTEXT",[(t,v) for t,v in sorted(complex_ids.items()) if t.startswith("ComplexPortal:")])
def sbml():
    p=core("R03")
    def local(x):return x.rsplit("}",1)[-1]
    # No extraction to disk; actual archive member and XML schema inspection.
    with tarfile.open(p,"r:gz") as tf:
        files=[m for m in tf.getmembers() if m.isfile()]
        emit("SBML_ARCHIVE_REGULAR_FILES",len(files))
        emit("SBML_ARCHIVE_SUFFIXES",dict(Counter(Path(m.name).suffix for m in files)))
        emit("SBML_ARCHIVE_EXAMPLES",[(m.name,m.size) for m in files[:8]])
        xml=sorted([m for m in files if Path(m.name).suffix.lower() in (".xml",".sbml")],key=lambda m:m.size)
        if not xml:
            emit("SBML_NO_RECOGNIZED_XML_SUFFIX; REVIEW_REQUIRED");return
        sample=[xml[i] for i in sorted({0,len(xml)//2,len(xml)-1})]
        for m in sample:
            with tf.extractfile(m) as f:root=ET.parse(f).getroot()
            counts=Counter(local(x.tag) for x in root.iter());attrs=defaultdict(set);examples=defaultdict(list)
            for x in root.iter():
                k=local(x.tag)
                attrs[k].update(local(a) for a in x.attrib)
                if k in ("model","species","reaction","speciesReference","modifierSpeciesReference","compartment","li") and len(examples[k])<5:
                    examples[k].append(dict(x.attrib))
            emit("SBML_SAMPLE",m.name,m.size,root.tag,dict(root.attrib))
            emit("SBML_TAG_COUNTS",sorted(counts.items()))
            emit("SBML_TAG_ATTRIBUTES",[(k,sorted(v)) for k,v in sorted(attrs.items())])
            emit("SBML_RECORD_EXAMPLES",dict(examples))
    emit("SBML_RECONNAISSANCE_COMPLETE; NOT_CROSS_FORMAT_CERTIFICATION")
def main():
    regen=verify_participants();hgnc_audit(regen);nonprotein();sbml()
    emit("P1A04S_REMAINING_IDENTITY_AND_FORMAT_PREFLIGHT_COMPLETE; NO_MERGES")
if __name__=="__main__":main()
