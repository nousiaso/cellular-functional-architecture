#!/usr/bin/env python3
"""Source-local state catalogue and structural fingerprint candidates, not canonical merges."""
import argparse
import csv
from collections import Counter, defaultdict
from contextlib import ExitStack
from decimal import Decimal
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
OUT=ROOT/"results/P1A04H_STATE_CATALOGUE_V1"
ENTITY={"Protein","Complex","PhysicalEntity","SmallMolecule","Dna","Rna"}
EVENT={"BiochemicalReaction","Degradation","TemplateReaction"}
START=time.monotonic()
def emit(*args): print(*args,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def js(x): return json.dumps(x,sort_keys=True,separators=(",",":"),ensure_ascii=False)
def fingerprint(x): return hashlib.sha256(js(x).encode("utf-8")).hexdigest()
def filehash(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def coefficient(v):
    d=Decimal(v)
    need(d.is_finite() and d>0,"invalid coefficient")
    return format(d.normalize(),"f")
def self_test():
    need(fingerprint({"component":[["a","2"]],"members":[]})!=fingerprint({"component":[],"members":["a"]}),"AND vs OR")
    need(fingerprint({"ref":["a"],"location":["x"]})!=fingerprint({"ref":["a"],"location":["y"]}),"location")
    need(fingerprint({"component":[["a","1"]]})!=fingerprint({"component":[["a","2"]]}),"stoichiometry")
    need(coefficient("2.00")=="2","coefficient normalization")
    emit("SELF_TEST_PASS")
def run(out):
    need(DB.stat().st_size==1451589632,"staging size mismatch")
    with ExitStack() as stack:
        db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
        stack.callback(db.close); db.execute("PRAGMA query_only=ON")
        need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","schema mismatch")
        def writer(name,cols):
            f=stack.enter_context((out/name).open("w",encoding="utf-8",newline=""))
            w=csv.writer(f,delimiter="\t",lineterminator="\n")
            w.writerow(cols.split("\t")); return w
        by_ord={}; cls={}
        for o,rid,uri in db.execute("SELECT ordinal,rdf_id,class_uri FROM bp_record"):
            by_ord[o]=rid
            if rid is not None: cls[rid]=uri.rsplit("}",1)[-1]
        need(len(by_ord)==548688,"record count"); emit("RECORDS",len(by_ord))
        keep={"xref","db","id","idVersion","organism","entityReference","cellularLocation",
              "feature","notFeature","memberPhysicalEntity","component","componentStoichiometry",
              "participantStoichiometry","physicalEntity","stoichiometricCoefficient",
              "left","right","template","product","displayName","modificationType",
              "featureLocation","sequencePosition","positionStatus","sequenceIntervalBegin","sequenceIntervalEnd","term"}
        props=defaultdict(lambda:defaultdict(list))
        for i,(o,p,t,v) in enumerate(db.execute("SELECT record_ordinal,predicate_uri,target_local_id,literal_text FROM bp_property"),1):
            k=p.rsplit("}",1)[-1]
            if k in keep and by_ord[o] is not None: props[by_ord[o]][k].append(t if t is not None else v)
            if i%500000==0: emit("PROPERTIES_SCANNED",i)
        need(i==2705676,"property count")
        def vals(r,k): return props.get(r,{}).get(k,[])
        @lru_cache(None)
        def xrefs(r):
            return tuple(sorted({(ns,acc) for x in vals(r,"xref") for ns in vals(x,"db") for acc in vals(x,"id") if ns and acc}))
        def stoich(r,role):
            out=defaultdict(list)
            for s in vals(r,role):
                e=vals(s,"physicalEntity"); v=vals(s,"stoichiometricCoefficient")
                need(len(e)==len(v)==1,"stoichiometry cardinality "+s)
                out[e[0]].append((s,coefficient(v[0])))
            return out
        @lru_cache(None)
        def state(r):
            c=cls[r]; need(c in ENTITY,"nonentity "+r)
            references=sorted(vals(r,"entityReference"))
            locations=sorted(vals(r,"cellularLocation"))
            features=sorted(vals(r,"feature"))
            absent_features=sorted(vals(r,"notFeature"))
            components=sorted(vals(r,"component"))
            members=sorted(vals(r,"memberPhysicalEntity"))
            st=stoich(r,"componentStoichiometry")
            flags=set(); taxa=set(); unresolved=False
            if not locations: flags.add("location_unrecorded")
            if not references and not components and not members:
                flags.add("opaque_leaf")
                unresolved=c!="SmallMolecule"
            for ref in references:
                t={acc for o in vals(ref,"organism") for ns,acc in xrefs(o) if ns=="NCBI Taxonomy"}
                taxa.update(t)
                if not t and c!="SmallMolecule": unresolved=True
            for f in features+absent_features:
                if cls[f]=="ModificationFeature":
                    if not vals(f,"modificationType"): flags.add("modification_type_unrecorded")
                    if not vals(f,"featureLocation"): flags.add("modification_position_unrecorded")
            children=[]
            for child in components:
                need(len(st[child])==1,"ambiguous component coefficient "+r)
                h,t,u,f=state(child)
                children.append((h,st[child][0][1]))
                taxa.update(t); unresolved|=u; flags.update(f)
            need(set(st)==set(components),"component stoich owner mismatch "+r)
            alternatives=[]
            for child in members:
                h,t,u,f=state(child)
                alternatives.append(h); taxa.update(t); unresolved|=u; flags.update(f)
            # Exact source reference, feature and location IDs, rather than unvalidated
            # biochemical equivalence of crossrefs or modification descriptions.
            # Opaque leaves are identity-anchored to avoid merging unknown entities.
            payload={"class":c,"reference_ids":references,"location_ids":locations,
                "feature_ids":features,"not_feature_ids":absent_features,
                "components":sorted(children),"members":sorted(alternatives),
                "opaque_anchor":r if not references and not children and not alternatives else None}
            return fingerprint(payload),frozenset(taxa),unresolved,frozenset(flags)
        entities=sorted(r for r,c in cls.items() if c in ENTITY)
        need(len(entities)==61418,"entity count")
        states=writer("biopax_states.tsv","rdf_id\tclass\tname_json\tfingerprint_candidate\treference_ids_json\treference_xrefs_json\tlocation_ids_json\tlocation_xrefs_json\tfeature_ids_json\tnot_feature_ids_json\ttaxon_union_json\tunresolved_taxon\tflags_json")
        comp=writer("biopax_composition.tsv","owner_id\tchild_id\trelation\tstoichiometry_id\tcoefficient")
        groups=defaultdict(list); flag_counts=Counter(); relations=Counter()
        for r in entities:
            h,t,u,f=state(r); groups[h].append(r); flag_counts.update(f)
            references=vals(r,"entityReference"); locations=vals(r,"cellularLocation")
            states.writerow([r,cls[r],js(vals(r,"displayName")),h,js(references),js([(x,xrefs(x)) for x in references]),
                js(locations),js([(x,xrefs(x)) for x in locations]),js(vals(r,"feature")),js(vals(r,"notFeature")),js(sorted(t)),int(u),js(sorted(f))])
            st=stoich(r,"componentStoichiometry")
            for child in vals(r,"component"):
                sid,v=st[child][0]; comp.writerow([r,child,"component",sid,v]); relations["component"]+=1
            for child in vals(r,"memberPhysicalEntity"):
                comp.writerow([r,child,"memberPhysicalEntity","",""]); relations["memberPhysicalEntity"]+=1
        need(relations=={"component":40083,"memberPhysicalEntity":36113},"composition count")
        dup=writer("state_fingerprint_candidates.tsv","fingerprint\tgroup_n\trdf_id\tclass\tname_json\tflags_json")
        for h,rs in sorted(groups.items()):
            if len(rs)>1:
                for r in rs: dup.writerow([h,len(rs),r,cls[r],js(vals(r,"displayName")),js(sorted(state(r)[3]))])
        emit("STATE_CATALOGUE",len(entities),"FINGERPRINTS",len(groups),
             "MULTIRECORD_GROUPS",sum(len(rs)>1 for rs in groups.values()),
             "GROUPED_RECORDS",sum(len(rs) for rs in groups.values() if len(rs)>1),"FLAGS",dict(flag_counts))
        emit("COMPOSITION_RELATIONS",dict(relations))
        # Preserve per-role incidence and coefficient records; if same entity occurs on
        # both sides, do not assign an owner-level stoichiometry to a particular side.
        part=writer("biopax_event_participants.tsv","event_rdf_id\tclass\trole\tentity_rdf_id\towner_stoichiometry_json\tcoefficient_status")
        review=writer("reaction_stoichiometry_review.tsv","event_rdf_id\tentity_rdf_id\troles_json\tstoichiometry_json\treason")
        pc=Counter(); inputs=defaultdict(set); outputs=defaultdict(set)
        for r,c in sorted(cls.items()):
            if c not in EVENT: continue
            st=stoich(r,"participantStoichiometry")
            roles=defaultdict(list)
            for role in ("left","right","template","product"):
                for e in vals(r,role): roles[e].append(role)
            for e,rs in roles.items():
                need(e in cls and cls[e] in ENTITY,"invalid participant "+str(e))
                s=st.get(e,[])
                status="multiple_records" if len(s)>1 else "side_unresolved" if s and len(set(rs))>1 else "explicit_owner_coefficient" if s else "unrecorded"
                if status in ("multiple_records","side_unresolved"):
                    review.writerow([r,e,js(rs),js(s),status])
                for role in rs:
                    part.writerow([r,c,role,e,js(s),status]); pc[(c,role,status)]+=1
                    if c=="BiochemicalReaction":
                        if role=="left": inputs[e].add(r)
                        if role=="right": outputs[e].add(r)
        emit("EVENT_PARTICIPANTS",sorted(pc.items()))
        # Exact source-ID bridges, no projection or removal. Record candidate weighting
        # inputs, not a selected network model.
        bridge=writer("raw_flow_bridges.tsv","entity_rdf_id\tname_json\tn_input_reactions\tn_output_reactions\tn_self_pairs\tn_ordered_nonself_pair_occurrences")
        total=0; nbridge=0
        for e in sorted(inputs.keys()&outputs.keys()):
            ni,no=len(inputs[e]),len(outputs[e]); ns=len(inputs[e]&outputs[e]); n=ni*no-ns
            bridge.writerow([e,js(vals(e,"displayName")),ni,no,ns,n]); total+=n; nbridge+=1
        need((nbridge,total)==(8098,345409),"raw bridge count mismatch")
        emit("RAW_FLOW_BRIDGES",nbridge,total)
        # GO model-scoped enablers: full structured object retained; no cross-source joins.
        go=writer("gocam_enablers.tsv","model_id\tactivity_ordinal\tactivity_id\tenabler_type\tenabler_term\tenabler_json")
        gn=0
        for mid,ordinal,sid,payload in db.execute("SELECT model_id,ordinal,source_id,payload_json FROM go_record WHERE field='activities'"):
            eb=json.loads(payload).get("enabled_by") or {}
            go.writerow([mid,ordinal,sid,eb.get("type"),eb.get("term"),js(eb)]); gn+=1
        need(gn==7385,"GO activity count")
        emit("GOCAM_ENABLERS",gn)
        return {"states":len(entities),"fingerprint_candidates":len(groups),"composition":dict(relations),"raw_bridge_entities":nbridge,"raw_bridge_pair_occurrences":total,"go_enablers":gn}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true"); ap.add_argument("--execute",action="store_true")
    a=ap.parse_args()
    if a.self_test: self_test(); return
    if not a.execute: emit("PLAN_ONLY",OUT.relative_to(ROOT)); return
    need(not OUT.exists(),"output exists; no overwrite")
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".state-catalogue-",dir=OUT.parent) as tmp:
        p=Path(tmp)/"payload"; p.mkdir()
        summary=run(p)
        manifest={"schema_version":"P1A04H_v1","scope":"source-local catalogue and candidate structural fingerprints; no merges",
                  "staging_hash_recomputed":False,"expected_staging_sha256":"9fc75c327a76de4ceb82c9095e9f5c82e34c6a202b1468f0836a53a498cdd63b",
                  "summary":summary,"files":{f.name:{"bytes":f.stat().st_size,"sha256":filehash(f)} for f in sorted(p.iterdir())}}
        (p/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared"); p.rename(OUT)
    emit("OUTPUT",OUT.relative_to(ROOT),summary)
    emit("P1A04H_STATE_CATALOGUE_PASS; NO_CANONICAL_MERGES")
if __name__=="__main__": main()
