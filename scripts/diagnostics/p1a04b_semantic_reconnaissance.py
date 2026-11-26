#!/usr/bin/env python3
"""P1A.04B read-only semantic reconnaissance; no graph, no extraction, no network."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET
import zipfile

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/"data/raw"
R=RAW/"reactome_v98"
G=RAW/"go_2026-08-05"
H=RAW/"hgnc_snapshot/identity/hgnc_complete_set.txt"
RDF="{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
XML="{http://www.w3.org/XML/1998/namespace}"
HSA=re.compile(r"R-HSA-[0-9]+")
def local(tag):
    return tag.rsplit("}",1)[-1]

def emit(label,*args):
    print(label,*args,flush=True)

def brief(x,n=600):
    return json.dumps(x,ensure_ascii=False,sort_keys=True)[:n]

def hgnc():
    mapping=defaultdict(set)
    status=Counter()
    with H.open(encoding="utf-8",newline="") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            status[row["status"]]+=1
            if row["status"]=="Approved":
                for acc in (row["uniprot_ids"] or "").split("|"):
                    if acc.strip():
                        mapping[acc.strip()].add(row["hgnc_id"])
    emit("HGNC_STATUS",dict(status))
    emit("HGNC_UNIPROT_KEYS",len(mapping),"AMBIGUOUS_KEYS",sum(len(v)>1 for v in mapping.values()))
    return mapping

def mapping_stats(label,accessions,mapping):
    accessions=set(accessions)
    emit("UNIPROT_EXACT_HGNC",label,"DISTINCT",len(accessions),
         "UNIQUE",sum(len(mapping.get(x,()))==1 for x in accessions),
         "AMBIGUOUS",sum(len(mapping.get(x,()))>1 for x in accessions),
         "UNMAPPED",sum(x not in mapping for x in accessions))
    emit("UNMAPPED_EXAMPLES",label,sorted(accessions-set(mapping))[:15])

def tables(mapping):
    role_pairs=set()
    events=set()
    proteins=set()
    roles=Counter()
    with (R/"extraction_check/ProteinRoleReaction.txt").open() as f:
        for acc,role,event in csv.reader(f,delimiter="\t"):
            if HSA.fullmatch(event):
                roles[role]+=1
                role_pairs.add((acc,role,event))
                events.add(event)
                proteins.add(acc)
    emit("R05_HSA_ROWS_BY_ROLE",dict(roles),"DISTINCT_TRIPLES",len(role_pairs),"EVENTS",len(events))
    mapping_stats("R05_HSA",proteins,mapping)
    uniprot_events=set()
    accession_set=set()
    evidence=Counter()
    species=Counter()
    with (R/"extraction_check/UniProt2ReactomeReactions.txt").open() as f:
        for acc,event,url,name,code,organism in csv.reader(f,delimiter="\t"):
            if HSA.fullmatch(event):
                uniprot_events.add(event)
                accession_set.add(acc)
                evidence[code]+=1
                species[organism]+=1
    emit("R07_HSA_EVENTS",len(uniprot_events),"EVIDENCE_ROWS",dict(evidence),"SPECIES_ROWS",dict(species))
    mapping_stats("R07_HSA",accession_set,mapping)
    emit("R05_R07_EVENT_OVERLAP",len(events & uniprot_events),"R05_ONLY",len(events-uniprot_events),"R07_ONLY",len(uniprot_events-events))
    r04=set()
    r04_roles=Counter()
    with (R/"extraction_check/reactome_reaction_exporter.txt").open() as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if HSA.fullmatch(row["reaction_id"]):
                r04.add(row["reaction_id"])
                r04_roles[row["role_in_reaction"]]+=1
    emit("R04_HSA_EVENTS",len(r04),"ROLE_STRINGS",r04_roles.most_common(25))
    refs=defaultdict(set)
    with (R/"provenance/ReactionPMIDS.txt").open() as f:
        for event,pmid in csv.reader(f,delimiter="\t"):
            if HSA.fullmatch(event):
                refs[event].add(pmid)
    emit("R06_HSA_EVENTS",len(refs),"EVENT_PMID_PAIRS",sum(map(len,refs.values())),
         "DISTINCT_PMIDS",len(set().union(*refs.values())) if refs else 0)
    return {"R04":r04,"R05":events,"R07":uniprot_events,"R06":set(refs)}

def gocam(mapping,table_events):
    directory=G/"human_indexed_models"
    inventory=json.loads((directory/"MODEL_INVENTORY.json").read_text())
    global_a=defaultdict(set)
    global_m=defaultdict(set)
    records=[]
    for receipt in inventory:
        mid=receipt["model"]["model_id"]
        obj=json.loads((directory/(mid+".json")).read_text())
        for field,index in (("activities",global_a),("molecules",global_m)):
            for entry in obj[field]:
                index[entry["id"]].add(mid)
        records.append((mid,obj))
    emit("GO_ACTIVITY_ID_UNIVERSE",len(global_a),"SHARED_ACROSS_MODELS",sum(len(x)>1 for x in global_a.values()))
    emit("GO_MOLECULE_ID_UNIVERSE",len(global_m),"SHARED_ACROSS_MODELS",sum(len(x)>1 for x in global_m.values()))
    all_reactome_ids=set()
    for ns in ("Reactome_HSA","hexadecimal"):
        counts=Counter()
        predicates=Counter()
        enabler=Counter()
        evidence_codes=Counter()
        reference_types=Counter()
        refs=set()
        providers=Counter()
        endpoints=Counter()
        accession=set()
        activity_ids=set()
        examples={}
        def scan(x,field):
            if isinstance(x,dict):
                if "evidence" in x:
                    ev=x["evidence"]
                    counts["evidence_lists"]+=1
                    if isinstance(ev,list):
                        counts["empty_evidence_lists"]+=not ev
                        for e in ev:
                            if not isinstance(e,dict):
                                counts["unexpected_evidence_item"]+=1
                                continue
                            evidence_codes[(field,str(e.get("term")))]+=1
                            ref=e.get("reference")
                            if isinstance(ref,str):
                                refs.add(ref)
                                reference_types[ref.split(":",1)[0]]+=1
                            else:
                                counts["evidence_reference_missing_or_nonstring"]+=1
                    else:
                        counts["unexpected_evidence_shape"]+=1
                for k,v in x.items():
                    if k=="provided_by" and isinstance(v,list):
                        providers.update(str(s) for s in v)
                    if k=="term" and isinstance(v,str) and v.startswith("UniProtKB:"):
                        accession.add(v.split(":",1)[1])
                    if k not in ("evidence","provenances"):
                        scan(v,field+"."+k)
                # Count provider occurrences at nested evidence and assertion levels too.
                for k in ("evidence","provenances"):
                    if k in x:
                        scan(x[k],field+"."+k)
            elif isinstance(x,list):
                for item in x:
                    scan(item,field+"[]")
        for mid,obj in records:
            if ("Reactome_HSA" if HSA.fullmatch(mid) else "hexadecimal")!=ns:
                continue
            counts["models"]+=1
            aa=[a["id"] for a in obj["activities"]]
            mm=[m["id"] for m in obj["molecules"]]
            aset,mset=set(aa),set(mm)
            counts["duplicate_activity_ids_within_model"]+=len(aa)-len(aset)
            counts["duplicate_molecule_ids_within_model"]+=len(mm)-len(mset)
            for a in obj["activities"]:
                counts["activities"]+=1
                activity_ids.add(a["id"])
                raw=a["id"].removeprefix("gomodel:")
                if HSA.fullmatch(raw):
                    all_reactome_ids.add(raw)
                eb=a.get("enabled_by",{})
                enabler[str(eb.get("type")) if isinstance(eb,dict) else type(eb).__name__]+=1
                if isinstance(eb,dict):
                    counts["activities_with_has_part"]+=bool(eb.get("has_part"))
                for field,key,local_ids,global_ids in (
                    ("causal_associations","downstream_activity",aset,global_a),
                    ("molecular_associations","molecule",mset,global_m)):
                    vals=a.get(field,[])
                    if not isinstance(vals,list):
                        counts["nonlist_"+field]+=1
                        continue
                    for rel in vals:
                        if not isinstance(rel,dict):
                            counts["nondict_"+field]+=1
                            continue
                        pred=str(rel.get("predicate"))
                        predicates[(field,pred)]+=1
                        target=rel.get(key)
                        state=("local" if isinstance(target,str) and target in local_ids else
                               "other_model_only" if isinstance(target,str) and target in global_ids else
                               "unresolved_or_invalid")
                        endpoints[(field,state)]+=1
                        examples.setdefault((field,pred),[mid,a["id"],brief(rel)])
                        if state!="local":
                            examples.setdefault((field,state),[mid,a["id"],brief(rel)])
                for field in ("enabled_by","molecular_function","occurs_in","part_of","causal_associations","molecular_associations"):
                    if field in a:
                        scan(a[field],field)
            scan(obj["molecules"],"molecules")
            scan(obj["provenances"],"model_provenances")
        emit("GO_COUNTS",ns,dict(counts),"DISTINCT_ACTIVITY_IDS",len(activity_ids))
        emit("GO_ENABLER_TYPES",ns,dict(enabler))
        emit("GO_PREDICATES",ns,predicates.most_common())
        emit("GO_ENDPOINT_RESOLUTION",ns,endpoints.most_common())
        emit("GO_EVIDENCE_CODES",ns,evidence_codes.most_common(35))
        emit("GO_REFERENCE_TYPES",ns,dict(reference_types),"DISTINCT_REFERENCES",len(refs))
        emit("GO_DISTINCT_PMIDS",ns,sum(x.startswith("PMID:") for x in refs))
        emit("GO_PROVIDER_OCCURRENCES",ns,providers.most_common(20))
        mapping_stats("GO_"+ns+"_ACTIVITIES_AND_MOLECULES",accession,mapping)
        for key,value in sorted(examples.items()):
            emit("GO_RELATION_EXAMPLE",ns,key,brief(value,850))
    emit("GO_EXACT_REACTOME_ACTIVITY_IDS",len(all_reactome_ids))
    for source,ids in table_events.items():
        emit("GO_REACTOME_ACTIVITY_ID_TABLE_OVERLAP",source,len(all_reactome_ids & ids),len(all_reactome_ids-ids))

def biopax(table_events,mapping):
    path=R/"primary/biopax.zip"
    classes=Counter()
    properties=Counter()
    ref_properties=Counter()
    literal_values=defaultdict(Counter)
    samples={}
    nodes={}
    refs=Counter()
    # Keep top-level records only; no element tree retained between records.
    # ID spellings are normalized only for obvious same-document # references.
    # Other unresolved forms are explicitly reported, not guessed.
    stable_by_class=defaultdict(set)
    xref_db=Counter()
    uniprot=set()
    with zipfile.ZipFile(path) as z, z.open("Homo_sapiens.owl") as f:
        stack=[]
        base=""
        for event,e in ET.iterparse(f,events=("start","end")):
            if event=="start":
                stack.append(e)
                if len(stack)==1:
                    base=e.attrib.get(XML+"base","")
                    emit("BIOPAX_ROOT_ATTR",e.attrib)
                continue
            if len(stack)==2:
                cls=local(e.tag)
                classes[cls]+=1
                ident=e.get(RDF+"ID") or e.get(RDF+"about")
                def norm(s):
                    if s.startswith("#"):
                        return s[1:]
                    if base and s.startswith(base.rstrip("#")+"#"):
                        return s[len(base.rstrip("#"))+1:]
                    return s
                if ident is not None:
                    key=norm(ident)
                    nodes.setdefault(key,Counter())[cls]+=1
                props=defaultdict(list)
                for child in e:
                    prop=local(child.tag)
                    properties[(cls,prop)]+=1
                    target=child.get(RDF+"resource")
                    val=(child.text or "").strip()
                    if target is not None:
                        ref_properties[(cls,prop)]+=1
                        refs[(cls,prop,norm(target))]+=1
                    if prop in ("conversionDirection","catalysisDirection","controlType"):
                        literal_values[(cls,prop)][val]+=1
                    if prop in ("db","id","displayName","left","right","controller","controlled","component","memberPhysicalEntity","entityReference","cellularLocation","feature","notFeature","physicalEntity","stoichiometricCoefficient","nextStep","stepProcess","xref","organism","controlType","conversionDirection","catalysisDirection"):
                        if len(props[prop])<5:
                            props[prop].append(target if target is not None else val)
                if cls in ("UnificationXref","RelationshipXref","PublicationXref"):
                    for db in props.get("db",[]):
                        xref_db[db]+=1
                        for val in props.get("id",[]):
                            if HSA.fullmatch(val):
                                stable_by_class[cls].add(val)
                            if db.lower() in ("uniprot","uniprotkb","uniprot knowledgebase"):
                                uniprot.add(val)
                if cls not in samples:
                    samples[cls]=(ident,dict(props))
                stack[-2].remove(e)
                e.clear()
            stack.pop()
    emit("BIOPAX_TOP_LEVEL_CLASSES",classes.most_common())
    emit("BIOPAX_NODE_IDS",len(nodes),"REPEATED_ID_RECORDS",sum(sum(v.values())-1 for v in nodes.values()))
    unresolved=Counter()
    resolved=Counter()
    examples={}
    for (cls,prop,target),n in refs.items():
        key=(cls,prop)
        if target in nodes:
            resolved[key]+=n
        else:
            unresolved[key]+=n
            examples.setdefault(key,target)
    for key in sorted(ref_properties):
        emit("BIOPAX_REFERENCE",key,"N",ref_properties[key],"RESOLVED",resolved[key],
             "UNRESOLVED",unresolved[key],"EXAMPLE_UNRESOLVED",examples.get(key,""))
    for key,val in sorted(literal_values.items()):
        emit("BIOPAX_LITERAL",key,dict(val))
    emit("BIOPAX_XREF_DATABASES",xref_db.most_common(40))
    for cls,ids in sorted(stable_by_class.items()):
        emit("BIOPAX_XREF_HSA",cls,len(ids),
             "TABLE_OVERLAP",{k:len(ids & v) for k,v in table_events.items()})
    mapping_stats("BIOPAX_XREF_EXACT_UNIPROT_DATABASE_LABELS",uniprot,mapping)
    for cls,(ident,props) in sorted(samples.items()):
        emit("BIOPAX_RECORD_EXAMPLE",cls,ident,brief(props,1200))
    emit("BIOPAX_TOP_LEVEL_PROPERTY_COUNTS",properties.most_common(80))
    emit("NATIVE_TOOL_AVAILABILITY",{k:bool(shutil.which(k)) for k in ("neo4j","neo4j-admin","java","zstd")})

def main():
    emit("CFA_P1A04B_SEMANTIC_RECONNAISSANCE_V1")
    mapping=hgnc()
    tableset=tables(mapping)
    gocam(mapping,tableset)
    biopax(tableset,mapping)
    emit("P1A04B_RECONNAISSANCE_COMPLETE; SEMANTIC_CERTIFICATION_PENDING")

if __name__=="__main__":
    main()
