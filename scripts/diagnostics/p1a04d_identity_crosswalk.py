#!/usr/bin/env python3
"""Read-only identity and event crosswalk audit; no canonicalization or graph."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import sqlite3
import time

ROOT=Path(__file__).resolve().parents[2]
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
RAW=ROOT/"data/raw"
BP="{http://www.biopax.org/release/biopax-level3.owl#}"
STABLE=re.compile(r"R-[A-Z]+-[0-9]+")
START=time.monotonic()
def emit(*args):
    print(*args, "ELAPSED_S",round(time.monotonic()-START,2),flush=True)

def main():
    if DB.stat().st_size!=1451589632:
        raise RuntimeError("Unexpected staged database size; inspect before proceeding")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    emit("CFA_P1A04D_IDENTITY_CROSSWALK_V1")
    meta=dict(db.execute("SELECT key,value FROM metadata"))
    emit("STAGING_METADATA",meta)
    if meta.get("schema_version")!="P1A04C_v1":
        raise RuntimeError("Unexpected staging schema")
    # Sequential scans rather than millions of indexed lookups on network storage.
    cls={}
    ord_to_id={}
    for ordinal,rid,uri in db.execute("SELECT ordinal,rdf_id,class_uri FROM bp_record"):
        cls[ordinal]=uri.rsplit("}",1)[-1]
        ord_to_id[ordinal]=rid
    id_to_ord={v:k for k,v in ord_to_id.items() if v is not None}
    emit("RECORD_INVENTORY",len(cls))
    keep={"xref","db","id","idVersion","organism","entityReference","cellularLocation",
          "feature","notFeature","memberPhysicalEntity","component","componentStoichiometry",
          "participantStoichiometry","physicalEntity","stoichiometricCoefficient",
          "left","right","controller","controlled","controlType","conversionDirection",
          "catalysisDirection","template","product","displayName","modificationType",
          "featureLocation","sequencePosition","positionStatus","sequenceIntervalBegin",
          "sequenceIntervalEnd","term"}
    props=defaultdict(lambda:defaultdict(list))
    for i,(ordinal,pred,target,text) in enumerate(db.execute(
        "SELECT record_ordinal,predicate_uri,target_local_id,literal_text FROM bp_property"),1):
        k=pred.rsplit("}",1)[-1]
        if k in keep:
            props[ordinal][k].append(target if target is not None else text)
        if i%500000==0:
            emit("PROPERTIES_SCANNED",i)
    emit("PROPERTY_SCAN_COMPLETE",i)
    def vals(o,k):
        return props.get(o,{}).get(k,[])
    def xrefs(o):
        out=[]
        for target in vals(o,"xref"):
            t=id_to_ord.get(target)
            if t is None:
                continue
            for namespace in vals(t,"db"):
                for accession in vals(t,"id"):
                    if namespace and accession:
                        out.append((cls[t],namespace,accession))
        return out
    def stable(o):
        return {acc for typ,ns,acc in xrefs(o) if typ=="UnificationXref"
                and ns=="Reactome" and STABLE.fullmatch(acc)}
    event_classes={"BiochemicalReaction","Degradation","TemplateReaction"}
    event_index=defaultdict(set)
    for c in sorted(event_classes):
        rows=[o for o,k in cls.items() if k==c]
        cardinality=Counter()
        for o in rows:
            ids=stable(o)
            cardinality[len(ids)]+=1
            for ident in ids:
                event_index[ident].add(o)
        emit("EVENT_STABLE_CARDINALITY",c,len(rows),dict(cardinality))
    emit("EVENT_STABLE_UNIVERSE",len(event_index),
         "ONE_RECORD",sum(len(v)==1 for v in event_index.values()),
         "MULTIPLE_RECORDS",sum(len(v)>1 for v in event_index.values()))
    r=RAW/"reactome_v98"
    tables={}
    with (r/"extraction_check/reactome_reaction_exporter.txt").open() as f:
        tables["R04"]={row["reaction_id"] for row in csv.DictReader(f,delimiter="\t")
                       if row["reaction_id"].startswith("R-HSA-")}
    for label,path,col in [
        ("R05","extraction_check/ProteinRoleReaction.txt",2),
        ("R06","provenance/ReactionPMIDS.txt",0),
        ("R07","extraction_check/UniProt2ReactomeReactions.txt",1)]:
        with (r/path).open() as f:
            tables[label]={row[col] for row in csv.reader(f,delimiter="\t")
                           if row[col].startswith("R-HSA-")}
    human={s for s in event_index if s.startswith("R-HSA-")}
    for label,ids in tables.items():
        emit("EVENT_EXPORT_CROSSWALK",label,"BOTH",len(human&ids),
             "BIOPAX_ONLY",len(human-ids),"EXPORT_ONLY",len(ids-human),
             "EXPORT_ONLY_EXAMPLES",sorted(ids-human)[:10])
    # State and molecular identity completeness; no default compartment or organism.
    mapping=defaultdict(set)
    with (RAW/"hgnc_snapshot/identity/hgnc_complete_set.txt").open() as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["status"]=="Approved":
                for acc in row["uniprot_ids"].split("|"):
                    if acc.strip():
                        mapping[acc.strip()].add(row["hgnc_id"])
    for c in ("ProteinReference","DnaReference","RnaReference","SmallMoleculeReference"):
        taxa=Counter()
        namespaces=Counter()
        uniprot=set()
        for o,k in cls.items():
            if k!=c:
                continue
            tax={acc for t in vals(o,"organism") if t in id_to_ord
                 for _,ns,acc in xrefs(id_to_ord[t]) if ns=="NCBI Taxonomy"}
            taxa[tuple(sorted(tax))]+=1
            for _,ns,acc in xrefs(o):
                namespaces[ns]+=1
                if ns in ("UniProt","UniProt Isoform"):
                    uniprot.add(acc)
        emit("REFERENCE_TAXON_SETS",c,taxa.most_common(30))
        emit("REFERENCE_XREF_NAMESPACES",c,namespaces.most_common(25))
        if uniprot:
            exact=Counter()
            for acc in uniprot:
                hits=mapping.get(acc,set())
                base_hits=mapping.get(acc.rsplit("-",1)[0],set()) if re.search(r"-[0-9]+$",acc) else set()
                exact["exact_unique" if len(hits)==1 else "exact_ambiguous" if hits else
                      "isoform_base_candidate_unique" if len(base_hits)==1 else
                      "isoform_base_candidate_ambiguous" if base_hits else "unmapped"]+=1
            emit("REFERENCE_HGNC_CANDIDATES",c,len(uniprot),dict(exact))
    for c in ("Protein","Complex","PhysicalEntity","SmallMolecule","Dna","Rna"):
        coverage=Counter()
        patterns=Counter()
        for o,k in cls.items():
            if k!=c:
                continue
            coverage["records"]+=1
            for p in ("entityReference","cellularLocation","feature","notFeature","component","memberPhysicalEntity"):
                coverage[p]+=bool(vals(o,p))
            patterns[(bool(vals(o,"component")),bool(vals(o,"memberPhysicalEntity")))]+=1
        emit("ENTITY_STATE_COVERAGE",c,dict(coverage),"COMPONENT_MEMBER_PATTERNS",patterns.most_common())
    # A named GO-CAM event ID is a candidate crosswalk, not a merged event.
    gc=Counter()
    examples=defaultdict(list)
    for model,ordinal,sid,payload in db.execute(
        "SELECT model_id,ordinal,source_id,payload_json FROM go_record WHERE field='activities'"):
        obj=json.loads(payload)
        raw=(sid or "").removeprefix("gomodel:")
        if STABLE.fullmatch(raw):
            n=len(event_index.get(raw,()))
            gc["bare_stable_occurrences"]+=1
            gc["bare_stable_unique" if n==1 else "bare_stable_ambiguous" if n else "bare_stable_unmatched"]+=1
            if n!=1 and len(examples["bare_review"])<10:
                examples["bare_review"].append((model,sid,n))
        else:
            gc["nonbare_activity_id"]+=1
        # References can identify pathways as well as reactions; never use them as equivalence.
        refs=set()
        def walk(x):
            if isinstance(x,dict):
                for k,v in x.items():
                    if k=="reference" and isinstance(v,str) and v.strip().startswith("Reactome:"):
                        refs.add(v.strip().split(":",1)[1])
                    walk(v)
            elif isinstance(x,list):
                for v in x:
                    walk(v)
        walk(obj)
        if refs:
            gc["activities_with_reactome_references"]+=1
            gc["activities_with_any_event_reference"]+=bool(refs & event_index.keys())
            if len(examples["reference_context"])<8:
                examples["reference_context"].append((model,sid,sorted(refs)))
    emit("GO_EVENT_CROSSWALK",dict(gc))
    emit("GO_EVENT_CROSSWALK_EXAMPLES",dict(examples))
    db.close()
    emit("P1A04D_READ_ONLY_AUDIT_COMPLETE; NO_CANONICAL_MERGES")

if __name__=="__main__":
    main()
