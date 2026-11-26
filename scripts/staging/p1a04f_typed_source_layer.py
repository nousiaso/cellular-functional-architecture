#!/usr/bin/env python3
"""Build source-scoped typed relationship tables; no canonical merging or inference."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
OUT=ROOT/"results/P1A04F_TYPED_SOURCE_LAYER_V1"
RAW=ROOT/"data/raw"
START=time.monotonic()
# label, direction relative to activity, qualitative effect (not a numerical weight)
RO={
"RO:0002233":("material_input","reverse",""),
"RO:0002234":("material_output","forward",""),
"RO:0012001":("small_molecule_activation","reverse","+"),
"RO:0012002":("small_molecule_inhibition","reverse","-"),
"RO:0002304":("upstream_positive_effect","forward","+"),
"RO:0002305":("upstream_negative_effect","forward","-"),
"RO:0002407":("indirect_positive_regulation","forward","+"),
"RO:0002409":("indirect_negative_regulation","forward","-"),
"RO:0002411":("causal_upstream_unspecified","forward",""),
"RO:0002413":("provides_input","forward",""),
"RO:0002578":("direct_regulation_unspecified","forward",""),
"RO:0002629":("direct_positive_regulation","forward","+"),
"RO:0002630":("direct_negative_regulation","forward","-"),
"RO:0012009":("constitutive_upstream","forward","+"),
"RO:0012010":("removes_input","forward","-")}
MOL={"RO:0002233","RO:0002234","RO:0012001","RO:0012002"}
EXPECTED=dict(zip(
["RO:0002233","RO:0002234","RO:0012001","RO:0012002","RO:0002304","RO:0002305",
 "RO:0002407","RO:0002409","RO:0002411","RO:0002413","RO:0002578","RO:0002629",
 "RO:0002630","RO:0012009","RO:0012010"],
[6641,2990,352,85,85,19,232,91,191,703,1,3166,844,11,80]))
EVENT={"BiochemicalReaction","Degradation","TemplateReaction"}
CONTROL={"Catalysis","Control","TemplateReactionRegulation"}
ENTITY={"Protein","Complex","PhysicalEntity","SmallMolecule","Dna","Rna"}
OBSOLETE={"GO:0006568","GO:0006904","GO:0039518","GO:0046774"}
STABLE=re.compile(r"R-[A-Z]+-[0-9]+")
def emit(*args):
    print(*args,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def local(s): return s.rsplit("}",1)[-1]
def writer(stack,path,header):
    f=stack.enter_context(path.open("w",encoding="utf-8",newline=""))
    w=csv.writer(f,delimiter="\t",lineterminator="\n")
    w.writerow(header.split("\t"))
    return w
def bp_role(cls,p):
    # Return typed source incidence; no flattening of controls, sets or steps.
    if cls=="BiochemicalReaction" and p in ("left","right"):
        return ("material_input","reverse") if p=="left" else ("material_output","forward")
    if cls=="Degradation" and p=="left": return "degradation_input","reverse"
    if cls=="TemplateReaction" and p in ("template","product"):
        return ("template","reverse") if p=="template" else ("template_product","forward")
    if cls in CONTROL and p in ("controller","controlled"):
        return ("controller_of_record","reverse") if p=="controller" else ("controlled_by_record","forward")
    if cls in ENTITY and p in ("component","memberPhysicalEntity","entityReference","cellularLocation","feature","notFeature","componentStoichiometry"):
        return {"component":"complex_component","memberPhysicalEntity":"set_member",
                "entityReference":"entity_reference","cellularLocation":"cellular_location",
                "feature":"feature","notFeature":"not_feature","componentStoichiometry":"component_stoichiometry"}[p],"forward"
    if cls in EVENT and p=="participantStoichiometry": return "participant_stoichiometry","forward"
    if cls=="Stoichiometry" and p=="physicalEntity": return "stoichiometry_entity","forward"
    if cls=="PathwayStep" and p in ("nextStep","stepProcess"):
        return ("curated_step_order" if p=="nextStep" else "step_process"),"forward"
    return None
def self_test():
    need(bp_role("BiochemicalReaction","left")==("material_input","reverse"),"input direction")
    need(bp_role("Catalysis","controller")==("controller_of_record","reverse"),"control direction")
    need(bp_role("Complex","component")!=bp_role("Complex","memberPhysicalEntity"),"set distinction")
    need(RO["RO:0002413"][2]=="" and RO["RO:0012010"][2]=="-","flow and removal")
    need(sum(EXPECTED.values())==15491 and len(RO)==15,"RO inventory")
    emit("SELF_TEST_PASS")
def run(out):
    from contextlib import ExitStack
    need(DB.stat().st_size==1451589632,"staged database size mismatch")
    for path,h in [
      (RAW/"ontology_interpretation/RO/2026-09-04/ro.obo","f2145eaecd649b395f5a908906f429fd6ac6dca2bc6816045b5c8a80da9a6b10"),
      (RAW/"ontology_interpretation/ECO/2026-07-10/eco.obo","c8b52a2da2b0f92224b50a1d8773fd7e3bca0977d920d4e31ef9fde960c84fd4")]:
        need(digest(path)==h,"ontology hash mismatch "+str(path))
    with ExitStack() as stack:
        db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
        stack.callback(db.close)
        db.execute("PRAGMA query_only=ON")
        need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","schema mismatch")
        edges=writer(stack,out/"biopax_incidence.tsv","source_id\tsource_class\ttarget_id\ttarget_class\tlayer\toriginal_predicate\trecord_ordinal\tproperty_ordinal")
        ops=writer(stack,out/"biopax_operations.tsv","rdf_id\tstable_id\tclass\tconversion_direction\tname")
        controls=writer(stack,out/"biopax_control_records.tsv","rdf_id\tclass\tcontrol_type\tcatalysis_direction")
        go_ops=writer(stack,out/"gocam_activities.tsv","model_id\tactivity_id\tactivity_ordinal\tmolecular_function\tenabler_type\tenabler_term\tprimary_model_taxon\tindexed_taxa_json")
        go_edges=writer(stack,out/"gocam_incidence.tsv","model_id\tsource_id\tsource_kind\ttarget_id\ttarget_kind\tlayer\tqualitative_effect\toriginal_predicate\tactivity_ordinal\tassociation_ordinal\tevidence_n")
        old=writer(stack,out/"gocam_obsolete_term_occurrences.tsv","model_id\tfield\trecord_ordinal\tsource_id\tjson_path\tterm")
        by_ord={}
        by_id={}
        for o,rid,uri in db.execute("SELECT ordinal,rdf_id,class_uri FROM bp_record"):
            c=local(uri); by_ord[o]=(rid,c)
            if rid is not None: by_id[rid]=(o,c)
        emit("BIOPAX_RECORDS",len(by_ord))
        need(len(by_ord)==548688,"record count")
        refs=defaultdict(list)
        literals=defaultdict(lambda:defaultdict(list))
        bc=Counter()
        for i,(o,po,uri,target,txt) in enumerate(db.execute(
            "SELECT record_ordinal,property_ordinal,predicate_uri,target_local_id,literal_text FROM bp_property"),1):
            rid,c=by_ord[o]; p=local(uri)
            if p=="xref" and c in EVENT: refs[o].append(target)
            if (c in EVENT and p in ("displayName","conversionDirection")) or (c in CONTROL and p in ("controlType","catalysisDirection")) or (c=="UnificationXref" and p in ("db","id")):
                literals[o][p].append(txt)
            role=bp_role(c,p)
            if role is not None:
                need(rid is not None and target in by_id,"unresolved typed BP incidence")
                tc=by_id[target][1]
                left,right,lc,rc=(target,rid,tc,c) if role[1]=="reverse" else (rid,target,c,tc)
                edges.writerow([left,lc,right,rc,role[0],p,o,po])
                bc[(c,p)]+=1
            if i%500000==0: emit("BIOPAX_PROPERTIES_SCANNED",i)
        need(i==2705676,"property count")
        stable_ids=set()
        for o,(rid,c) in by_ord.items():
            if c in EVENT:
                ids=set()
                for ref in refs[o]:
                    if ref not in by_id or by_id[ref][1]!="UnificationXref": continue
                    lp=literals[by_id[ref][0]]
                    if "Reactome" in lp.get("db",[]):
                        ids.update(x for x in lp.get("id",[]) if x and STABLE.fullmatch(x))
                need(len(ids)==1,"stable event identity "+str(rid))
                sid=next(iter(ids)); need(sid not in stable_ids,"duplicate stable event "+sid)
                stable_ids.add(sid)
                direction=literals[o].get("conversionDirection",[])
                if c=="BiochemicalReaction": need(direction==["LEFT-TO-RIGHT"],"unreviewed conversion direction")
                ops.writerow([rid,sid,c,json.dumps(direction),json.dumps(literals[o].get("displayName",[]),ensure_ascii=False)])
            elif c in CONTROL:
                controls.writerow([rid,c,json.dumps(literals[o].get("controlType",[])),json.dumps(literals[o].get("catalysisDirection",[]))])
        need(len(stable_ids)==16187,"event count")
        for k,n in {("BiochemicalReaction","left"):31282,("BiochemicalReaction","right"):24986,
            ("Complex","component"):40083,("Complex","componentStoichiometry"):40083,
            ("Catalysis","controller"):6871,("Catalysis","controlled"):6871,
            ("Control","controller"):3271,("Control","controlled"):3271,
            ("PathwayStep","nextStep"):15528,("PathwayStep","stepProcess"):29192}.items():
            need(bc[k]==n,"incidence count "+str(k))
        emit("BIOPAX_TYPED_INCIDENCE",sorted(bc.items()),"TOTAL",sum(bc.values()))
        emit("BIOPAX_OPERATIONS",len(stable_ids))
        models={m:(t,a) for m,t,a in db.execute("SELECT model_id,primary_taxon,indexed_taxa_json FROM go_model")}
        activities=set(); molecules=set(); obsolete=Counter()
        def scan(x,path,context):
            if isinstance(x,dict):
                for k,v in x.items():
                    pp=path+"."+k
                    if k=="term" and isinstance(v,str) and v in OBSOLETE:
                        old.writerow([*context,pp,v]); obsolete[v]+=1
                    scan(v,pp,context)
            elif isinstance(x,list):
                for i,v in enumerate(x): scan(v,path+"["+str(i)+"]",context)
        for mid,field,ordinal,sid,payload in db.execute("SELECT model_id,field,ordinal,source_id,payload_json FROM go_record"):
            obj=json.loads(payload)
            scan(obj,"$",[mid,field,ordinal,sid])
            if field=="activities":
                activities.add((mid,sid))
                eb=obj.get("enabled_by") or {}; mf=obj.get("molecular_function") or {}
                go_ops.writerow([mid,sid,ordinal,mf.get("term"),eb.get("type"),eb.get("term"),*models[mid]])
            elif field=="molecules": molecules.add((mid,sid))
        need(len(activities)==7385 and len(molecules)==7322,"GO scoped ID count")
        counts=Counter()
        for mid,ao,aid,field,ordinal,pred,tid,payload in db.execute(
            "SELECT model_id,activity_ordinal,activity_id,field,ordinal,predicate,target_id,payload_json FROM go_association"):
            need(pred in RO,"unreviewed RO predicate "+str(pred))
            mol=field=="molecular_associations"
            need(field in ("molecular_associations","causal_associations") and mol==(pred in MOL),"RO endpoint type mismatch")
            need((mid,aid) in activities and (mid,tid) in (molecules if mol else activities),"GO endpoint mismatch")
            layer,direction,effect=RO[pred]
            target_kind="molecule" if mol else "activity"
            left,right,lk,rk=(tid,aid,target_kind,"activity") if direction=="reverse" else (aid,tid,"activity",target_kind)
            ev=json.loads(payload).get("evidence")
            go_edges.writerow([mid,left,lk,right,rk,layer,effect,pred,ao,ordinal,len(ev) if isinstance(ev,list) else ""])
            counts[pred]+=1
        need(dict(counts)==EXPECTED,"GO predicate counts changed")
        emit("GOCAM_TYPED_INCIDENCE",sorted(counts.items()),"TOTAL",sum(counts.values()))
        emit("GOCAM_OBSOLETE_TERM_OCCURRENCES_IN_RECORDS",dict(obsolete))
    return {"biopax_typed_incidence":sum(bc.values()),"biopax_operations":len(stable_ids),
            "gocam_activities":len(activities),"gocam_typed_incidence":sum(counts.values()),
            "obsolete_term_occurrences":dict(obsolete)}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--execute",action="store_true")
    a=ap.parse_args()
    if a.self_test: self_test(); return
    if not a.execute: emit("PLAN_ONLY",OUT.relative_to(ROOT)); return
    need(not OUT.exists(),"output exists; do not overwrite")
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".typed-source-",dir=OUT.parent) as tmp:
        p=Path(tmp)/"payload"; p.mkdir()
        summary=run(p)
        manifest={"schema_version":"P1A04F_v1","scope":"source-local typed incidence only; no canonical merges",
                  "expected_staging_sha256":"9fc75c327a76de4ceb82c9095e9f5c82e34c6a202b1468f0836a53a498cdd63b",
                  "staging_hash_recomputed":False,"summary":summary,
                  "files":{f.name:{"bytes":f.stat().st_size,"sha256":digest(f)} for f in sorted(p.iterdir())}}
        (p/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared")
        p.rename(OUT)
    emit("OUTPUT",OUT.relative_to(ROOT),summary)
    emit("P1A04F_TYPED_SOURCE_LAYER_PASS; NO_CANONICAL_MERGES")
if __name__=="__main__": main()
