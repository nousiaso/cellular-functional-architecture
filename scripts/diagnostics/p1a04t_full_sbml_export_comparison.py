#!/usr/bin/env python3
"""Full SBML/BioPAX export comparison; no merging or source changes."""
from collections import Counter,defaultdict
from decimal import Decimal
import json
from pathlib import Path
import re
import sqlite3
import tarfile
import xml.etree.ElementTree as ET
from p1a04s_remaining_identity_and_format_preflight import core,ROOT,STATE,DB,rows,sha,need
from p1a04k_scope_feasibility_audit import TYPED,check_manifest

S="{http://www.sbml.org/sbml/level3/version1/core}"
R="{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
B="{http://biomodels.net/biology-qualifiers/}"
STABLE=re.compile(r"R-[A-Z]+-[0-9]+")
RELEASE_NS="Reactome Database ID Release 98"
STAGING_HASH="9fc75c327a76de4ceb82c9095e9f5c82e34c6a202b1468f0836a53a498cdd63b"
def emit(*x):print(*x,flush=True)
def numeric(s,prefix):
    m=re.fullmatch(prefix+r"_([0-9]+)",s or "")
    return m.group(1) if m else None
def annotations(x):
    # Direct is annotations only; hasPart/isHomologTo/isDescribedBy are not identity.
    a=x.find(S+"annotation")
    if a is None:return set()
    meta=x.get("metaid")
    if not meta:return set()
    out=set()
    for desc in a.iter(R+"Description"):
        if desc.get(R+"about")!="#"+meta:continue
        for q in desc.findall(B+"is"):
            for li in q.iter(R+"li"):
                url=li.get(R+"resource","")
                m=re.fullmatch(r"https?://identifiers\.org/(?:reactome[:/])(R-[A-Z]+-[0-9]+)(?:\.[0-9]+)?",url)
                if m:out.add(m.group(1))
    return out
def bp_index():
    need(sha(DB)==STAGING_HASH,"staged DB hash mismatch")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True);db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","staged schema")
    records={o:(i,c.rsplit("}",1)[-1]) for o,i,c in db.execute("SELECT ordinal,rdf_id,class_uri FROM bp_record")}
    xp=defaultdict(lambda:defaultdict(set));links=defaultdict(set)
    for o,k,t,v in db.execute("SELECT record_ordinal,predicate_uri,target_local_id,literal_text FROM bp_property"):
        k=k.rsplit("}",1)[-1]
        if records[o][1]=="UnificationXref" and k in ("db","id") and v is not None:xp[records[o][0]][k].add(v)
        elif k=="xref" and t is not None:links[records[o][0]].add(t)
    db.close()
    stable=defaultdict(set);internal=defaultdict(set);ns=Counter();examples=defaultdict(list)
    for o,targets in links.items():
        for t in targets:
            for n in xp[t]["db"]:
                if n not in ("Reactome",RELEASE_NS):continue
                for ident in xp[t]["id"]:
                    ns[n]+=1
                    if len(examples[n])<5:examples[n].append((o,ident))
                    if n=="Reactome" and STABLE.fullmatch(ident):stable[ident].add(o)
                    if n==RELEASE_NS and re.fullmatch(r"[0-9]+",ident):internal[ident].add(o)
    emit("BP_DIRECT_XREF_NAMESPACES",dict(ns),"EXAMPLES",dict(examples))
    emit("BP_DIRECT_XREF_KEYS",len(stable),len(internal))
    return stable,internal
def resolve(ids,num,stable,internal,allowed):
    st=set().union(*(stable.get(x,set()) for x in ids))&allowed
    it=internal.get(num,set())&allowed if num else set()
    if ids and not st:return None,"direct_stable_unmatched"
    if len(ids)>1:return None,"multiple_stable_annotations"
    if st and it and not st&it:return None,"stable_numeric_conflict"
    hit=st&it if st and it else st or it
    if len(hit)==1:return next(iter(hit)),"direct_stable" if st else "release_numeric_candidate"
    return None,"ambiguous" if hit else "unmapped"
def main():
    check_manifest(TYPED,"P1A04F_v1");check_manifest(STATE,"P1A04H_v1")
    ops={x["rdf_id"]:x for x in rows(TYPED/"biopax_operations.tsv")}
    states={x["rdf_id"]:x for x in rows(STATE/"biopax_states.tsv")}
    need(len(ops)==16187 and len(states)==61418,"BP inventory")
    stable,internal=bp_index()
    # Confirm that direct stable xrefs identify the same operation catalogue.
    need(all(stable.get(x["stable_id"],set())&ops.keys()=={e} for e,x in ops.items()),"BP event stable index")
    allowed_events=set(ops);allowed_states=set(states)
    # Restrict indexes up front to avoid building large sets per participant.
    st_ev={k:v&allowed_events for k,v in stable.items()};st_ent={k:v&allowed_states for k,v in stable.items()}
    in_ev={k:v&allowed_events for k,v in internal.items()};in_ent={k:v&allowed_states for k,v in internal.items()}
    part=defaultdict(lambda:defaultdict(set));coeff={}
    for x in rows(STATE/"biopax_event_participants.tsv"):
        e=x["event_rdf_id"];role={"left":"input","right":"output","template":"template","product":"product"}[x["role"]];s=x["entity_rdf_id"]
        part[e][role].add(s)
        if x["coefficient_status"]=="explicit_owner_coefficient":
            vals=json.loads(x["owner_stoichiometry_json"]);need(len(vals)==1,"BP coefficient schema")
            coeff[e,role,s]=Decimal(vals[0][1])
    controls={x["rdf_id"]:x for x in rows(TYPED/"biopax_control_records.tsv")}
    ctr=defaultdict(set);target=defaultdict(set)
    for x in rows(TYPED/"biopax_incidence.tsv"):
        if x["layer"]=="controller_of_record":ctr[x["target_id"]].add(x["source_id"])
        elif x["layer"]=="controlled_by_record":target[x["source_id"]].add(x["target_id"])
    for c,x in controls.items():
        role="catalyst" if x["class"]=="Catalysis" else "positiveregulator" if x["control_type"]=="ACTIVATION" else "negativeregulator" if x["control_type"]=="INHIBITION" else "other_modifier"
        for e in target[c]&allowed_events:part[e][role].update(ctr[c])
    # Aggregate repeated pathway appearances separately, and compare distinct
    # event signatures only once. Do not treat export repetitions as evidence.
    counts=Counter();maps=Counter();species_maps=Counter();sb_roles=Counter();examples=defaultdict(list)
    seen=defaultdict(set);mapped=defaultdict(set);comparison=Counter();comparison_events=defaultdict(set);missing=[]
    def example(key,*v):
        if len(examples[key])<8:examples[key].append(v)
    def resolve_fast(ids,num,st,it):
        # Index already restricted to the appropriate entity type.
        st_hits=set().union(*(st.get(x,set()) for x in ids))
        it_hits=it.get(num,set())
        if ids and not st_hits:return None,"direct_stable_unmatched"
        if len(ids)>1:return None,"multiple_stable_annotations"
        if st_hits and it_hits and not st_hits&it_hits:return None,"stable_numeric_conflict"
        hit=st_hits&it_hits if st_hits and it_hits else st_hits or it_hits
        if len(hit)==1:return next(iter(hit)),"direct_stable" if st_hits else "release_numeric_candidate"
        return None,"ambiguous" if hit else "unmapped"
    path=core("R03")
    with tarfile.open(path,"r|gz") as tf:
        for m in tf:
            if not m.isfile():continue
            counts["files"]+=1
            need(m.name.endswith(".sbml"),"unexpected archive member "+m.name)
            with tf.extractfile(m) as f:root=ET.parse(f).getroot()
            need(root.tag==S+"sbml" and root.get("level")=="3" and root.get("version")=="1","unexpected SBML schema "+m.name)
            model=root.find(S+"model");need(model is not None,"missing model")
            sp={}
            for x in model.findall(S+"listOfSpecies/"+S+"species"):
                sid=x.get("id");need(sid and sid not in sp,"duplicate species")
                ids=annotations(x);num=numeric(sid,"species")
                e,k=resolve_fast(ids,num,st_ent,in_ent)
                sp[sid]=(e,k)
                species_maps[k]+=1
                if e is None:example("species_"+k,m.name,sid,sorted(ids),x.get("name"))
            for x in model.findall(S+"listOfReactions/"+S+"reaction"):
                rid=x.get("id");need(rid,"reaction ID")
                ids=annotations(x);num=numeric(rid,"reaction")
                e,k=resolve_fast(ids,num,st_ev,in_ev)
                counts["reaction_occurrences"]+=1;maps[k]+=1
                roles=defaultdict(list)
                for listtag,tag,role in (("listOfReactants","speciesReference","input"),("listOfProducts","speciesReference","output"),("listOfModifiers","modifierSpeciesReference","modifier")):
                    for y in x.findall(S+listtag+"/"+S+tag):
                        sid=y.get("species");need(sid in sp,"unknown species reference "+str(sid))
                        z=role
                        if z=="modifier":
                            sm={"SBO:0000013":"catalyst","SBO:0000459":"positiveregulator","SBO:0000020":"negativeregulator"}
                            z=sm.get(y.get("sboTerm"),"other_modifier:"+str(y.get("sboTerm")))
                        roles[z].append((sid,y.get("stoichiometry"),sp[sid][0],sp[sid][1]))
                        sb_roles[z]+=1
                signature=(x.get("reversible"),tuple(sorted((z,tuple(sorted((sid,str(v)) for sid,v,_,_ in vs))) for z,vs in roles.items())))
                key=("mapped",e) if e else ("unmapped",rid,tuple(sorted(ids)))
                fresh=signature not in seen[key];seen[key].add(signature)
                if e is None:example("reaction_"+k,m.name,rid,sorted(ids),x.get("name"));continue
                mapped[e].add(k)
                if not fresh:continue
                cl=ops[e]["class"]
                comparison[(cl,"unique_SBML_signature")]+=1
                for z in ("input","output","catalyst","positiveregulator","negativeregulator"):
                    vs=roles[z]
                    if cl=="TemplateReaction" and z in ("input","output"):
                        tag="template_semantics_not_comparable"
                    elif any(v[2] is None for v in vs):tag="SBML_species_unresolved"
                    else:
                        ss={v[2] for v in vs};bs=part[e][z]
                        tag="equal" if ss==bs else "SBML_subset" if ss<bs else "SBML_superset" if ss>bs else "different"
                        if tag!="equal":example("role_"+z+"_"+tag,ops[e]["stable_id"],sorted(bs-ss)[:10],sorted(ss-bs)[:10])
                        if z in ("input","output"):
                            for sid,v,bid,_ in vs:
                                if bid not in bs:continue
                                if v is None:ck="SBML_coefficient_unrecorded"
                                elif (e,z,bid) not in coeff:ck="BP_coefficient_not_explicit_or_side_unresolved"
                                else:ck="equal_explicit" if Decimal(v)==coeff[e,z,bid] else "different_explicit"
                                comparison[(cl,z,ck)]+=1
                                if ck=="different_explicit":example("coefficient_difference",ops[e]["stable_id"],z,bid,v,str(coeff[e,z,bid]))
                    comparison[(cl,z,tag)]+=1;comparison_events[(cl,z,tag)].add(e)
                comparison[(cl,"reversible",str(x.get("reversible")))]+=1
            if counts["files"]%250==0:emit("SBML_PROGRESS",counts["files"],counts["reaction_occurrences"])
    need(counts["files"]==2889,"SBML archive file count")
    emit("SBML_TOTALS",dict(counts),"UNIQUE_REACTION_KEYS",len(seen))
    emit("SBML_REACTION_MAPPING_OCCURRENCES",dict(maps))
    emit("SBML_SPECIES_MAPPING_OCCURRENCES",dict(species_maps))
    emit("SBML_ROLE_OCCURRENCES",dict(sb_roles))
    emit("BP_EVENT_COVERAGE",len(mapped),"BP_ONLY",len(ops)-len(mapped),"UNMAPPED_SBML_KEYS",sum(k[0]=="unmapped" for k in seen))
    emit("SBML_EVENT_MATCH_ROUTE",dict(Counter(tuple(sorted(v)) for v in mapped.values())))
    emit("SBML_SIGNATURE_VARIANTS",dict(Counter(len(v) for v in seen.values())))
    emit("SBML_ROLE_AND_COEFFICIENT_COMPARISON_SIGNATURES",sorted(comparison.items()))
    emit("SBML_ROLE_COMPARISON_DISTINCT_EVENTS",sorted((k,len(v)) for k,v in comparison_events.items()))
    emit("SBML_BP_ONLY_BY_CLASS",dict(Counter(ops[e]["class"] for e in allowed_events-mapped.keys())))
    emit("SBML_BP_ONLY_EXAMPLES",[(ops[e]["stable_id"],ops[e]["name"]) for e in sorted(allowed_events-mapped.keys())[:15]])
    for k,v in sorted(examples.items()):emit("REVIEW_EXAMPLES",k,v)
    emit("P1A04T_FULL_SBML_EXPORT_COMPARISON_COMPLETE; NO_MERGES; REVIEW_DISCREPANCIES")
if __name__=="__main__":main()
