#!/usr/bin/env python3
"""P1B04: cross-table integrity certificate; no graph fitting or scope lock."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import csv
import json
from pathlib import Path
import subprocess
import tempfile
from p1b03_scope_provenance_benchmark import (
    ROOT, B01, B02, OUT as B03, B01_HASH, B02_HASH, verify, rows, sha,
    need, js, emit, file_for, references, norm_pmid, evidence_terms
)
from p1b02_molecular_compatibility import loc_profile, STRICT

OUT=ROOT/"results/P1B04_HARMONIZATION_CERTIFICATE_V1"
B03_HASH="953f8eae0970ea86b7f4fae300d7b9015a0956e98cd3b61647a8882a6a503d28"

def table(path,key):
    d={}
    for r in rows(path):
        k=tuple(r[x] for x in key) if isinstance(key,tuple) else r[key]
        need(k not in d,"duplicate "+path.name+" "+repr(k))
        d[k]=r
    return d

def arr(r,k):return set(json.loads(r[k]))

def verify_chain():
    idx,lh=verify()
    p=B03/"manifest.json"
    need(sha(p)==B03_HASH,"P1B03 manifest hash")
    m=json.loads(p.read_text())
    need(m["schema_version"]=="P1B03_v1" and m["source_lock_sha256"]==lh
         and m["P1B01_manifest_sha256"]==B01_HASH and m["P1B02_manifest_sha256"]==B02_HASH,"P1B03 lineage")
    for n,x in m["files"].items():
        need(Path(n).name==n,"unsafe name")
        f=B03/n
        need(f.stat().st_size==x["bytes"] and sha(f)==x["sha256"],"P1B03 file "+n)
    need(len(m["files"])==9,"P1B03 output inventory")
    emit("P1B03_MANIFEST_VERIFIED",len(m["files"]))
    return idx,lh

def certify(out,idx):
    import re
    nodes=table(B01/"nodes.tsv","node_id")
    kinds=Counter((r["source"],r["kind"]) for r in nodes.values())
    need(len(nodes)==235140,"node count")
    state={e for e,r in nodes.items() if r["source"]=="BP" and r["kind"]=="state"}
    ops={e for e,r in nodes.items() if r["source"]=="BP" and r["kind"]=="operation"}
    acts={e for e,r in nodes.items() if r["source"]=="GO" and r["kind"]=="activity"}
    need((len(state),len(ops),len(acts))==(61418,16187,7385),"primary nodes")
    seen=set(); counts=Counter();go={}; ins=defaultdict(set);outs=defaultdict(set)
    for r in rows(B01/"assertions.tsv"):
        a,b=r["from_id"],r["to_id"];aid=r["assertion_id"]
        need(a in nodes and b in nodes and aid not in seen,"assertion integrity "+aid)
        seen.add(aid);counts[r["source"]]+=1
        if r["source"]=="GO":
            go[aid]=r
            need(nodes[a]["model_id"]==nodes[b]["model_id"],"GO cross-model assertion")
        if r["source"]=="BP" and r["layer"]=="material_input":
            need(a in state and b in ops,"material input endpoints");ins[a].add(b)
        if r["source"]=="BP" and r["layer"]=="material_output":
            need(a in ops and b in state,"material output endpoints");outs[b].add(a)
    need(dict(counts)=={"BP":441156,"GO":15491,"GO_ENABLER":7385},"assertion counts")
    del seen
    emit("ASSERTION_ENDPOINTS_PASS",sum(counts.values()))
    def bridge(a,b):return len(outs[a])*len(ins[b])-len(outs[a]&ins[b])
    need(sum(bridge(e,e) for e in state)==345409,"exact bridge count")
    payload={e:json.loads(nodes[e]["payload_json"]) for e in state}
    loc={e:arr(r,"location_ids_json") for e,r in payload.items()}
    gloc={e:{t for _,xs in json.loads(r["location_xrefs_json"]) for _,t in xs if re.fullmatch(r"GO:[0-9]+",t)} for e,r in payload.items()}
    feat={e:(arr(r,"feature_ids_json"),arr(r,"not_feature_ids_json")) for e,r in payload.items()}
    comp=table(B01/"composition.tsv",("owner_id","child_id","relation"))
    interfaces=table(B02/"composition_interfaces.tsv",("owner_id","child_id","relation"))
    need(comp.keys()==interfaces.keys() and len(comp)==76196,"composition conservation")
    members=defaultdict(set);components=defaultdict(set);ctx={}
    for k,r in interfaces.items():
        a,b,rel=k;orig=comp[k]
        need(a in state and b in state and a!=b,"composition endpoints")
        need(rel in ("component","memberPhysicalEntity"),"composition kind")
        need(all(r[x]==orig[x] for x in ("stoichiometry_id","coefficient")),"coefficient preservation")
        need(r["location_status"]==loc_profile(loc[a],loc[b],gloc[a],gloc[b]),"location profile")
        need(r["feature_id_status"]==("equal_recorded_feature_ids" if feat[a]==feat[b] else "different_recorded_feature_ids"),"feature profile")
        need(int(r["child_to_owner_material_occurrences"])==bridge(b,a) and int(r["owner_to_child_material_occurrences"])==bridge(a,b),"bridge count")
        if rel=="component":
            need(r["interpretation"]=="required_component_not_material_transfer","AND interpretation")
            components[a].add(b)
        else:
            need(r["interpretation"]=="alternative_realization","OR interpretation")
            members[a].add(b);ctx[a,b]=r["location_status"]
    need(not set(members)&set(components),"mixed composition")
    @lru_cache(None)
    def realizations(e):
        if e not in members:return frozenset([(e,True)])
        return frozenset((leaf,good and ctx[e,c] in STRICT) for c in members[e] for leaf,good in realizations(c))
    expected={(e,leaf,"location_supported" if good else "location_uncertain") for e in members for leaf,good in realizations(e)}
    got=set()
    for r in rows(B02/"or_realizations.tsv"):
        a,b,status=r["set_id"],r["terminal_id"],r["location_path_class"]
        k=(a,b,status);path=json.loads(r["witness_path_json"])
        need(k not in got and k in expected,"OR witness key")
        need(path[0]==a and path[-1]==b and len(path)-1==int(r["depth"]) and b not in members,"OR path endpoints")
        need(all(y in members.get(x,set()) for x,y in zip(path,path[1:])),"OR path edge")
        good=all(ctx[x,y] in STRICT for x,y in zip(path,path[1:]))
        need(status==("location_supported" if good else "location_uncertain"),"OR witness context")
        need(int(r["feature_ids_differ_along_witness"])==int(any(feat[x]!=feat[y] for x,y in zip(path,path[1:]))),"OR witness feature")
        got.add(k)
    need(got==expected and len({(a,b) for a,b,_ in got})==40324,"OR completeness")
    assemblies=table(B02/"complex_production_candidates.tsv",("complex_id","operation_id"))
    need(assemblies.keys()=={(e,r) for e in components for r in outs[e]},"complex production coverage")
    for (e,r),x in assemblies.items():
        matched={c for c in components[e] if r in ins[c]}
        missing=components[e]-matched
        need(arr(x,"matched_direct_components_json")==matched and arr(x,"unmatched_direct_components_json")==missing,"complex input witnesses")
    emit("MOLECULAR_WITNESSES_PASS",len(got),len(assemblies))
    # Validate factored identity endpoints. Reference identity is not state equivalence.
    for r in rows(B02/"identifier_memberships.tsv"):
        need(r["node_id"] in nodes,"identifier node")
    for r in rows(B02/"cross_source_molecule_candidates.tsv"):
        need(r["go_node_id"] in nodes and nodes[r["go_node_id"]]["kind"]=="molecule","GO molecule identity")
        groups=[arr(r,k) for k in ("BP_exact_location_ids_json","BP_missing_location_ids_json","BP_other_location_ids_json")]
        need(set.union(*groups)<=state and sum(map(len,groups))==len(set.union(*groups)),"cross-source candidate partition")
    for r in rows(B02/"gene_candidates.tsv"):
        if r["source_kind"]=="BP_reference":
            need(r["reference_or_token"] in nodes,"BP gene reference endpoint")
    pathways=table(B03/"pathways.tsv","pathway_id")
    need(len(pathways)==2889,"pathway count")
    stable={e:next(iter(arr(r,"stable_id_candidates_json"))) for e,r in pathways.items() if len(arr(r,"stable_id_candidates_json"))==1}
    need(len(stable)==len(set(stable.values()))==2889,"pathway stable IDs")
    parents=defaultdict(set);direct=defaultdict(set);step_owners=defaultdict(set);step_ops=defaultdict(set);native_counts=Counter()
    prov=set()
    for r in rows(B03/"native_pathway_assertions.tsv"):
        a,b,k=r["from_id"],r["to_id"],r["predicate"]
        key=(r["record_ordinal"],r["property_ordinal"])
        need(key not in prov,"native pathway assertion provenance");prov.add(key);native_counts[k]+=1
        if k=="pathwayComponent":
            need(a in pathways and (b in pathways or b in ops),"pathway component")
            if b in pathways:parents[b].add(a)
            else:direct[b].add(a)
        elif k=="pathwayOrder":
            need(a in pathways,"pathway step owner");step_owners[b].add(a)
        elif k=="stepProcess":
            if b in ops:step_ops[a].add(b)
        elif k=="nextStep":pass
        else:raise RuntimeError("unknown native pathway predicate "+k)
    need(dict(native_counts)=={"pathwayComponent":19996,"pathwayOrder":19996,"stepProcess":29192,"nextStep":15528},"native path counts")
    need(len(direct)==16187,"native event coverage")
    @lru_cache(None)
    def ancestors(e):
        return frozenset(parents[e]).union(*(ancestors(x) for x in parents[e]))
    for e,r in pathways.items():
        need(arr(r,"parent_ids_json")==parents[e],"pathway parents")
        roots={x for x in ancestors(e)|{e} if not parents[x]}
        need(arr(r,"root_ids_json")==roots,"pathway roots")
    steps=defaultdict(set)
    for s,es in step_ops.items():
        for e in es:steps[e].update(step_owners[s])
    membership=table(B03/"event_pathway_benchmark.tsv","operation_id")
    need(membership.keys()==ops,"event pathway coverage")
    # Native disease context is an additional benchmark, not an eligibility filter.
    with file_for(idx,"R11").open(newline="",encoding="utf-8") as f:
        disease={a for a,b in csv.reader(f,delimiter="\t") if re.fullmatch(r"R-[A-Z]+-[0-9]+",a)}
    def dclass(ps):
        return "no_membership" if not ps else "disease_only" if ps<=disease else "mixed" if ps&disease else "no_listed_disease_context"
    changes=Counter();native_disease=Counter()
    with (out/"native_disease_context.tsv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["operation_id","native_direct_disease_context","native_direct_pathways_with_disease_ancestor_json","R04_disease_context"])
        for e,r in sorted(membership.items()):
            ds=direct[e];anc=set().union(*(ancestors(x) for x in ds))
            ss={stable[x] for x in ds}
            need(arr(r,"native_direct_json")==ds and arr(r,"native_ancestor_json")==anc,"native membership")
            need(arr(r,"native_root_json")=={x for x in ds|anc if not parents[x]},"event roots")
            need(arr(r,"step_process_candidate_json")==steps[e],"step membership")
            need(arr(r,"native_direct_stable_json")==ss and arr(r,"R04_stable_pathways_json")<=ss,"R04 containment")
            d=dclass(ss);native_disease[d]+=1;changes[r["R04_disease_context"],d]+=1
            inherited={stable[x] for x in ds if {stable[y] for y in ancestors(x)|{x}}&disease}
            w.writerow([e,d,js(sorted(inherited)),r["R04_disease_context"]])
    emit("PATHWAY_BENCHMARK_PASS",dict(native_disease))
    sc=table(B03/"state_scope.tsv","state_id");ec=table(B03/"event_scope.tsv","operation_id")
    need(sc.keys()==state and ec.keys()==ops,"scope coverage")
    need(Counter(r["composition_feasibility"] for r in sc.values())=={"H":55616,"U":3275,"F":2527},"state scope counts")
    need(Counter(r["catalyst_profile"] for r in ec.values())=={"H":6092,"none_recorded":9323,"U":516,"F":256},"catalyst scope counts")
    pm=defaultdict(set);pub_pairs=set();pub_prov=set();pub_n=Counter()
    for r in rows(B03/"bp_publication_links.tsv"):
        e,x=r["owner_id"],r["xref_id"];key=(r["record_ordinal"],r["property_ordinal"])
        need(e in ops and key not in pub_prov and (e,x) not in pub_pairs,"BP publication provenance")
        pub_prov.add(key);pub_pairs.add((e,x))
        ns=arr(r,"namespaces_json");ii=arr(r,"identifiers_json")
        pub_n[r["publication_identity_status"]]+=1
        if any(isinstance(a,str) and a.lower()=="pubmed" for a in ns):
            pm[e].update(a for a in ii if isinstance(a,str) and re.fullmatch(r"[0-9]+",a))
    need(len(pub_pairs)==60015 and pub_n=={"namespace_and_identifier_recorded":59951,"missing_namespace_and_identifier":64},"publication identity")
    need(sum(bool(x) for x in pm.values())==15701,"PMID coverage")
    all_bp=set().union(*pm.values())
    ga=table(B03/"go_activity_provenance.tsv","activity_id");ge=table(B03/"go_assertion_provenance.tsv","assertion_id")
    need(ga.keys()==acts and ge.keys()==go.keys(),"GO provenance coverage")
    for e,r in ga.items():
        obj=json.loads(nodes[e]["payload_json"]);refs=set(references(obj))
        ps={x for s in refs if (x:=norm_pmid(s))}
        need(arr(r,"all_references_json")==refs and arr(r,"PMIDs_json")==ps and arr(r,"BP_publication_overlap_json")==ps&all_bp,"GO activity evidence")
    strata=Counter();candidate_layers=Counter();candidate_models=set();candidate_pmids=set()
    for e,r in ge.items():
        original=go[e];ev=json.loads(original["evidence_json"]);pr=json.loads(original["provenance_json"])
        refs=set(references(ev));ps={x for s in refs if (x:=norm_pmid(s))}
        need(r["model_id"]==pr["model_id"] and r["layer"]==original["layer"],"GO assertion provenance identity")
        need(arr(r,"ECO_terms_json")==evidence_terms(ev) and arr(r,"references_json")==refs and arr(r,"PMIDs_json")==ps,"GO assertion evidence")
        need(arr(r,"BP_publication_overlap_json")==ps&all_bp,"GO publication overlap")
        ap=set().union(*(arr(ga[a],"PMIDs_json") for a in (original["from_id"],original["to_id"]) if a in ga))&all_bp
        need(arr(r,"endpoint_activity_PMIDs_overlap_json")==ap,"GO activity publication overlap")
        tag=r["corroboration_stratum"];strata[tag]+=1
        if tag=="candidate_corroboration_no_observed_dependence_marker":
            need(all(r[k]=="0" for k in ("import_marker","automatic_inference_marker","endpoint_or_model_Reactome_dependency")) and ps and not ps&all_bp and not ap,"corroboration exclusion marker")
            candidate_layers[pr["field"],r["layer"]]+=1
            candidate_models.add(r["model_id"]);candidate_pmids.update(ps)
    need(strata["candidate_corroboration_no_observed_dependence_marker"]==5573 and sum(strata.values())==15491,"GO stratum counts")
    emit("PROVENANCE_CROSS_TABLE_PASS",len(pub_pairs),len(ge))
    return {"node_count":len(nodes),"assertions":dict(counts),"composition_links":len(comp),
        "OR_witness_rows":len(got),"complex_production_candidates":len(assemblies),
        "native_pathways":len(pathways),"native_event_memberships":len(membership),
        "native_disease_context":dict(native_disease),"R04_native_disease_comparison":[[list(k),v] for k,v in sorted(changes.items())],
        "BP_publication_pairs":len(pub_pairs),"GO_strata":dict(strata),
        "GO_candidate_layers":[[list(k),v] for k,v in sorted(candidate_layers.items())],
        "GO_candidate_models":len(candidate_models),"GO_candidate_PMIDs":len(candidate_pmids),
        "scope_filters":0,"graph_edges_created":0}

def main():
    p=argparse.ArgumentParser();p.add_argument("--execute",action="store_true");a=p.parse_args()
    idx,lh=verify_chain()
    if not a.execute:emit("PLAN_ONLY",OUT.relative_to(ROOT));return
    need(not OUT.exists(),"certificate exists; no overwrite")
    with tempfile.TemporaryDirectory(prefix=".p1b04_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir()
        summary=certify(tmp,idx);emit("P1B04_COUNTS",summary)
        m={"schema_version":"P1B04_v1","status":"TECHNICAL_PASS_REVIEW_REQUIRED",
           "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
           "source_lock_sha256":lh,"P1B01_manifest_sha256":B01_HASH,"P1B02_manifest_sha256":B02_HASH,
           "P1B03_manifest_sha256":B03_HASH,"summary":summary,
           "files":{p.name:{"bytes":p.stat().st_size,"sha256":sha(p)} for p in sorted(tmp.iterdir())},
           "limits":["Cross-table consistency, not independent biological validation","Native disease context is annotation only","GO candidate flags do not establish independent evidence","Evidence split feasibility, domain coverage, projection weights and primary scope choice remain open"]}
        (tmp/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"certificate appeared");tmp.rename(OUT)
    emit("P1B04_MANIFEST_SHA256",sha(OUT/"manifest.json"))
    emit("P1B04_TECHNICAL_PASS_REVIEW_REQUIRED")
if __name__=="__main__":main()
