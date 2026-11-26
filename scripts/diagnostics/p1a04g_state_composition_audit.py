#!/usr/bin/env python3
"""Read-only molecular state, nested composition and stoichiometry audit."""
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from functools import lru_cache
import json
from pathlib import Path
import sqlite3
import time

ROOT=Path(__file__).resolve().parents[2]
DB=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
ENTITY={"Protein","Complex","PhysicalEntity","SmallMolecule","Dna","Rna"}
EVENT={"BiochemicalReaction","Degradation","TemplateReaction"}
CONTROL={"Catalysis","Control","TemplateReactionRegulation"}
START=time.monotonic()
def emit(*args):
    print(*args,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def main():
    need(DB.stat().st_size==1451589632,"unexpected DB size")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","schema mismatch")
    by_ord={}
    cls={}
    for o,rid,uri in db.execute("SELECT ordinal,rdf_id,class_uri FROM bp_record"):
        by_ord[o]=rid
        if rid is not None: cls[rid]=uri.rsplit("}",1)[-1]
    emit("RECORDS",len(by_ord))
    keep={"xref","db","id","idVersion","organism","entityReference","cellularLocation",
          "feature","notFeature","memberPhysicalEntity","component","componentStoichiometry",
          "participantStoichiometry","physicalEntity","stoichiometricCoefficient",
          "left","right","controller","controlled","controlType","template","product",
          "displayName","modificationType","featureLocation","sequencePosition",
          "positionStatus","sequenceIntervalBegin","sequenceIntervalEnd","term"}
    prop=defaultdict(lambda:defaultdict(list))
    for i,(o,p,t,v) in enumerate(db.execute(
        "SELECT record_ordinal,predicate_uri,target_local_id,literal_text FROM bp_property"),1):
        k=p.rsplit("}",1)[-1]
        if k in keep and by_ord[o] is not None:
            prop[by_ord[o]][k].append(t if t is not None else v)
        if i%500000==0: emit("PROPERTIES_SCANNED",i)
    need(i==2705676,"property count")
    def vals(r,k): return prop.get(r,{}).get(k,[])
    @lru_cache(None)
    def xrefs(r):
        return tuple(sorted({(ns,acc) for x in vals(r,"xref")
                             for ns in vals(x,"db") for acc in vals(x,"id")
                             if ns and acc}))
    @lru_cache(None)
    def reference_taxa(r):
        return frozenset(acc for o in vals(r,"organism")
                         for ns,acc in xrefs(o) if ns=="NCBI Taxonomy")
    entities=sorted(r for r,c in cls.items() if c in ENTITY)
    composition={r:tuple(vals(r,"component")+vals(r,"memberPhysicalEntity")) for r in entities}
    comp_counts=Counter()
    for r,children in composition.items():
        for child in children:
            comp_counts[(cls[r],cls.get(child,"UNRESOLVED"))]+=1
            need(child in composition,"unexpected composition target "+str(child))
    emit("COMPOSITION_ENDPOINT_CLASSES",sorted(comp_counts.items()))
    # Detect cycles without assuming tree structure. Shared children are not cycles.
    visiting=set()
    done={}
    cycles=set()
    def depth(r):
        if r in done: return done[r]
        if r in visiting:
            cycles.add(r); return None
        visiting.add(r)
        ds=[depth(c) for c in composition[r]]
        visiting.remove(r)
        ans=None if any(d is None for d in ds) else 1+max(ds,default=-1)
        done[r]=ans
        return ans
    for r in entities: depth(r)
    emit("COMPOSITION_CYCLE_BACKEDGE_IDS",len(cycles),sorted(cycles)[:15])
    emit("COMPOSITION_DEPTH",dict(Counter(str(d) for d in done.values())))
    # Union of known leaf taxon labels is descriptive; not logical feasibility of a set.
    @lru_cache(None)
    def taxa(r):
        if done[r] is None: return frozenset(),True
        known=set()
        unknown=False
        for ref in vals(r,"entityReference"):
            t=reference_taxa(ref)
            known.update(t)
            if not t and cls[r]!="SmallMolecule": unknown=True
        if not vals(r,"entityReference") and not composition[r] and cls[r]!="SmallMolecule":
            unknown=True
        for child in composition[r]:
            t,u=taxa(child); known.update(t); unknown|=u
        return frozenset(known),unknown
    def tax_summary(r):
        t,u=taxa(r)
        label=("mixed_known" if "9606" in t and len(t)>1 else "human_known" if t=={"9606"}
               else "nonhuman_known" if t else "no_known_taxon")
        return label+("_with_unresolved" if u else "")
    emit("ENTITY_TAXON_COMPOSITION",sorted(Counter((cls[r],tax_summary(r)) for r in entities).items()))
    ec=Counter()
    for r,c in cls.items():
        if c not in EVENT|CONTROL: continue
        for p in ("left","right","template","product","controller"):
            for e in vals(r,p):
                ec[(c,p,tax_summary(e) if e in composition else "other_target_class:"+cls.get(e,"UNKNOWN"))]+=1
    emit("EVENT_ROLE_TAXON_COMPOSITION",sorted(ec.items()))
    # Stoichiometric records: validate explicit coefficients and owner membership.
    coeff=Counter(); st_status=Counter(); st_examples=defaultdict(list)
    stoich={r for r,c in cls.items() if c=="Stoichiometry"}
    for r in stoich:
        es=vals(r,"physicalEntity"); vs=vals(r,"stoichiometricCoefficient")
        status="valid_positive"
        if len(es)!=1 or len(vs)!=1: status="cardinality"
        else:
            try:
                v=Decimal(vs[0])
                if not v.is_finite() or v<=0: status="nonpositive_or_nonfinite"
                else: coeff[str(v)]+=1
            except (InvalidOperation,TypeError): status="non_numeric"
        st_status[status]+=1
        if status!="valid_positive" and len(st_examples[status])<5:
            st_examples[status].append((r,es,vs))
    emit("STOICHIOMETRY_RECORDS",len(stoich),dict(st_status),"COEFFICIENTS",coeff.most_common(30),"EXAMPLES",dict(st_examples))
    owner_stats=Counter(); ownership=Counter(); owner_examples=defaultdict(list)
    for r,c in cls.items():
        if c not in ENTITY|EVENT: continue
        for p,role in (("componentStoichiometry","component"),("participantStoichiometry","participant")):
            members=vals(r,"component") if role=="component" else vals(r,"left")+vals(r,"right")
            refs=vals(r,p)
            if not members and not refs: continue
            covered=Counter()
            for s in refs:
                ownership[s]+=1
                targets=vals(s,"physicalEntity")
                if s not in stoich or len(targets)!=1 or targets[0] not in members:
                    owner_stats[(c,p,"invalid_owner_target")]+=1
                    if len(owner_examples["invalid_owner_target"])<8:
                        owner_examples["invalid_owner_target"].append((r,s,targets))
                else:
                    covered[targets[0]]+=1
            owner_stats[(c,p,"owners")]+=1
            owner_stats[(c,p,"members_distinct")]+=len(set(members))
            owner_stats[(c,p,"explicitly_covered_members")]+=len(covered)
            owner_stats[(c,p,"members_without_explicit_coefficient")]+=len(set(members)-covered.keys())
            owner_stats[(c,p,"multiple_stoich_refs_same_member")]+=sum(n>1 for n in covered.values())
    emit("STOICHIOMETRY_OWNER_COVERAGE",sorted(owner_stats.items()))
    emit("STOICHIOMETRY_OWNERSHIP","unreferenced",len(stoich-set(ownership)),"multi_owner",sum(v>1 for v in ownership.values()),"examples",dict(owner_examples))
    fc=Counter()
    for r,c in cls.items():
        if c not in ("ModificationFeature","FragmentFeature"): continue
        fc[(c,"records")]+=1
        locs=vals(r,"featureLocation")
        fc[(c,"has_location")]+=bool(locs)
        fc[(c,"has_modification_type")]+=bool(vals(r,"modificationType"))
        for l in locs:
            fc[(c,"location_class:"+cls.get(l,"UNKNOWN"))]+=1
    emit("FEATURE_STRUCTURE",sorted(fc.items()))
    lc=Counter()
    for r in entities:
        ls=vals(r,"cellularLocation")
        loc_ids={acc for l in ls for ns,acc in xrefs(l) if ns=="GENE ONTOLOGY"}
        lc[(cls[r],"GO_location_count",len(loc_ids))]+=1
    emit("ENTITY_GO_LOCATION_COVERAGE",sorted(lc.items()))
    # Raw exact state ID reuse only: not a canonical or projected reaction graph.
    inputs=defaultdict(set); outputs=defaultdict(set)
    for r,c in cls.items():
        if c!="BiochemicalReaction": continue
        for e in vals(r,"left"): inputs[e].add(r)
        for e in vals(r,"right"): outputs[e].add(r)
    shared=inputs.keys() & outputs.keys()
    ranks=[]
    for e in shared:
        n=len(inputs[e])*len(outputs[e])-len(inputs[e]&outputs[e])
        ranks.append((n,e,len(inputs[e]),len(outputs[e]),vals(e,"displayName")))
    ranks.sort(reverse=True)
    emit("RAW_FLOW_BRIDGE_ENTITY_COUNT",len(shared))
    emit("RAW_FLOW_PAIR_OCCURRENCES_NO_SELF",sum(x[0] for x in ranks))
    emit("RAW_FLOW_TOP20_BRIDGES",ranks[:20])
    # GO-CAM enabler structure; preserve all nested evidence, no protein set expansion.
    gc=Counter(); ex=defaultdict(list)
    for mid,ordinal,sid,payload in db.execute(
        "SELECT model_id,ordinal,source_id,payload_json FROM go_record WHERE field='activities'"):
        obj=json.loads(payload); eb=obj.get("enabled_by") or {}
        typ=eb.get("type","MISSING"); term=eb.get("term") or ""
        gc[(typ,"records")]+=1
        gc[(typ,"term_namespace:"+term.split(":",1)[0])]+=1
        for k in ("has_part","part_of"):
            if k in eb:
                gc[(typ,k,"present")]+=1
                gc[(typ,k,"nonempty")]+=bool(eb[k])
                if len(ex[(typ,k)])<2:
                    ex[(typ,k)].append((mid,sid,json.dumps(eb[k],ensure_ascii=False)[:700]))
    emit("GO_ENABLER_STRUCTURE",sorted(gc.items()))
    for k,v in sorted(ex.items()): emit("GO_ENABLER_EXAMPLES",k,v)
    db.close()
    emit("P1A04G_STATE_COMPOSITION_AUDIT_COMPLETE; NO_CANONICAL_MERGES")
if __name__=="__main__": main()
