#!/usr/bin/env python3
"""Read-only Part 1D feasibility. No fitted model, splits or null samples."""
from collections import Counter, defaultdict
import importlib.util
import json
import random
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"staging"))
from p1c02_expanded_compatibility import ROOT, B01, B02, C01, OUT as C02, C01_HASH, verify, rows, sha, need, emit
from p1c01_native_exact_network import B03

C02_HASH="c7bc4a83b0ec28dd450e16f6cd69cc3b6de0b6356285a8d4e04284605e7f0534"
def arr(r,k):return set(json.loads(r[k]))
class DSU:
    def __init__(self,n):self.p=list(range(n));self.sz=[1]*n
    def find(self,i):
        while self.p[i]!=i:
            self.p[i]=self.p[self.p[i]];i=self.p[i]
        return i
    def union(self,i,j):
        i,j=self.find(i),self.find(j)
        if i==j:return
        if self.sz[i]<self.sz[j]:i,j=j,i
        self.p[j]=i;self.sz[i]+=self.sz[j]
    def report(self,ids):
        c=Counter(self.find(i) for i in ids)
        return {"components":len(c),"largest_10":sorted(c.values(),reverse=True)[:10],"singletons":sum(v==1 for v in c.values())}

def main():
    lh,bh=verify()
    p=C02/"manifest.json";need(sha(p)==C02_HASH,"C02 manifest")
    m=json.loads(p.read_text())
    need(m["schema_version"]=="P1C02_v1" and m["source_lock_sha256"]==lh and m["harmonization_lock_sha256"]==bh and m["P1C01_manifest_sha256"]==C01_HASH,"C02 lineage")
    for name,x in m["files"].items():
        need(Path(name).name==name,"unsafe filename");f=C02/name
        need(f.stat().st_size==x["bytes"] and sha(f)==x["sha256"],"C02 file "+name)
    need(len(m["files"])==5,"C02 inventory")
    emit("P1C02_MANIFEST_VERIFIED",len(m["files"]))
    E={r["state_id"]:int(r["state_index"]) for r in rows(C01/"state_index.tsv")}
    R={r["operation_id"]:int(r["operation_index"]) for r in rows(C01/"operation_index.tsv")}
    n=len(R);ne=len(E);need((n,ne)==(16187,61418),"universe")
    rclass={R[r["operation_id"]]:r["class"] for r in rows(C01/"operation_index.tsv")}
    scope={E[r["state_id"]]:r["composition_feasibility"] for r in rows(B03/"state_scope.tsv")}
    escope={R[r["operation_id"]]:r["candidate_scope_stratum"] for r in rows(B03/"event_scope.tsv")}
    stateprof={};state_name={}
    for r in rows(B01/"nodes.tsv"):
        if r["node_id"] not in E:continue
        e=E[r["node_id"]];x=json.loads(r["payload_json"])
        stateprof[e]=(r["class"],scope[e],tuple(sorted(arr(x,"location_ids_json"))))
        state_name[e]=json.loads(x["name_json"])
    need(len(stateprof)==ne,"state metadata")
    paths={r["pathway_id"]:json.loads(r["name_json"]) for r in rows(B03/"pathways.tsv")}
    root={};direct={};ancestor={}
    for r in rows(B03/"event_pathway_benchmark.tsv"):
        i=R[r["operation_id"]];root[i]=arr(r,"native_root_json");direct[i]=arr(r,"native_direct_json")
        ancestor[i]=direct[i]|arr(r,"native_ancestor_json")
    need(len(root)==n and all(root.values()) and all(direct.values()),"pathway profiles")
    rootops=defaultdict(set)
    for i,rs in root.items():
        for a in rs:rootops[a].add(i)
    emit("ROOT_PROFILE_COUNTS",{"roots":len(rootops),"distinct_root_profiles":len({tuple(sorted(x)) for x in root.values()}),"distinct_direct_profiles":len({tuple(sorted(x)) for x in direct.values()}),"multi_root_operations":sum(len(x)>1 for x in root.values()),"multi_direct_operations":sum(len(x)>1 for x in direct.values())})
    edges=defaultdict(set);by_state=defaultdict(set);by_op=defaultdict(set);native=DSU(n+ne);incident=set()
    for r in rows(C01/"typed_incidence.tsv"):
        e,i,k=int(r["state_index"]),int(r["operation_index"]),r["role"]
        need((e,i) not in edges[k] and r["binary_value"]=="1","native duplicate")
        edges[k].add((e,i));native.union(i,n+e);incident.add(e)
        by_state[e].add(i);by_op[i].add(e)
    emit("NATIVE_ROLE_CONNECTIVITY",native.report(range(n)),"used_states",len(incident))
    for r in rows(B01/"composition.tsv"):
        native.union(n+E[r["owner_id"]],n+E[r["child_id"]])
    emit("NATIVE_WITH_COMPOSITION_CONNECTIVITY",native.report(range(n)))
    locations=defaultdict(set)
    for r in rows(C02/"context_incidence.tsv"):
        locations[int(r["operation_index"])].add(r["location_id"])
    exact=defaultdict(set);weight=defaultdict(float);summary=defaultdict(Counter)
    for r in rows(C01/"exact_projection_candidates.tsv"):
        k=r["layer"];i,j=int(r["from_operation_index"]),int(r["to_operation_index"])
        need(i!=j and (i,j) not in exact[k],"exact edge duplicate");exact[k].add((i,j))
        summary[k]["pairs"]+=1
        summary[k]["same_direct"]+=bool(direct[i]&direct[j])
        summary[k]["same_ancestor"]+=bool(ancestor[i]&ancestor[j])
        summary[k]["same_root"]+=bool(root[i]&root[j])
        summary[k]["same_context"]+=bool(locations[i]&locations[j])
        summary[k]["witness_occurrences"]+=int(r["unweighted_count"])
    need(len(exact["flow"])==339542,"flow")
    expanded=defaultdict(set)
    for r in rows(C02/"expanded_operation_pairs.tsv"):
        k=r["layer"];i,j=int(r["from_operation_index"]),int(r["to_operation_index"])
        need(i!=j and (i,j) not in expanded[k],"expanded duplicate");expanded[k].add((i,j))
        ws=[tuple(x) for x in json.loads(r["state_pair_witnesses_json"])]
        need(len(ws)==len(set(ws))==int(r["n_unique_state_pairs"]),"expanded witness multiplicity")
        need(all((a,i) in edges["material_output"] and (b,j) in edges["material_input"] for a,b in ws),"expanded witness incidence")
        need(int(r["exact_flow_pair"])==int((i,j) in exact["flow"]),"exact overlap")
        summary[k]["pairs"]+=1;summary[k]["same_direct"]+=bool(direct[i]&direct[j])
        summary[k]["same_ancestor"]+=bool(ancestor[i]&ancestor[j])
        summary[k]["same_root"]+=bool(root[i]&root[j])
        summary[k]["same_context"]+=bool(locations[i]&locations[j])
        summary[k]["witness_occurrences"]+=len(ws)
    need(len(set().union(*expanded.values())-exact["flow"])==77716,"expanded added-only union")
    for k,z in sorted({**exact,**expanded}.items()):
        d=DSU(n)
        for i,j in z:d.union(i,j)
        emit("LAYER_FEASIBILITY",k,dict(summary[k]),"weak_components",d.report(range(n)))
    # Root coverage is annotation only. Overlapping counts are intentional.
    for a,os in sorted(rootops.items(),key=lambda x:(-len(x[1]),x[0])):
        cv={k:len(os&{i for ij in z for i in ij}) for k,z in {**exact,**expanded}.items()}
        emit("ROOT_COVERAGE",a,paths[a],len(os),cv)
    # Evidence group upper bound: all event records sharing a numeric PMID
    # or recorded namespace/id pair. Unknown references do not connect events.
    pm=defaultdict(set);pub=defaultdict(set)
    for r in rows(B03/"bp_publication_links.tsv"):
        i=R[r["owner_id"]];ns=arr(r,"namespaces_json");ids=arr(r,"identifiers_json")
        for a in ns:
            for b in ids:
                if isinstance(a,str) and a and isinstance(b,str) and b:
                    pub[a.lower(),b].add(i)
                    if a.lower()=="pubmed" and b.isascii() and b.isdigit():pm[b].add(i)
    for name,groups in (("numeric_PMID",pm),("recorded_publication_ID",pub)):
        d=DSU(n);observed=set()
        for es in groups.values():
            a=min(es);observed.update(es)
            for b in es:d.union(a,b)
        emit("EVIDENCE_COMPONENTS",name,{"identities":len(groups),"events_with_identity":len(observed),"events_without_identity":n-len(observed),"largest_identity_groups":sorted(map(len,groups.values()),reverse=True)[:10],"event_components_all":d.report(range(n)),"event_components_known":d.report(observed)})
    # Local proposal feasibility on observed incidence, not a Markov chain.
    # Same-role swaps (e,i),(f,j)->(e,j),(f,i) preserve binary margins.
    # N1 adds exact state class/location/HFU and event class/catalyst scope;
    # N2 adds exact full root profile; N3 adds exact direct membership profile.
    # These are deliberately restrictive candidate strata; not frozen nulls.
    rng=random.Random(20261127);samples=5000
    for level in range(4):
        buckets=defaultdict(list);all_edges=[]
        for k,z in sorted(edges.items()):
            for e,i in sorted(z):
                key=(k,)
                if level>=1:key+=(stateprof[e],rclass[i],escope[i])
                if level>=2:key+=(tuple(sorted(root[i])),)
                if level>=3:key+=(tuple(sorted(direct[i])),)
                buckets[key].append((e,i));all_edges.append((key,e,i))
        c=Counter()
        for _ in range(samples):
            key,e,i=rng.choice(all_edges);pool=buckets[key]
            if len(pool)==1:c["singleton_stratum"]+=1;continue
            f,j=rng.choice(pool)
            if e==f or i==j:c["same_state_or_operation"]+=1;continue
            if (e,j) in edges[key[0]] or (f,i) in edges[key[0]]:c["existing_cross_edge"]+=1;continue
            c["admissible_local_swap"]+=1
        emit("NULL_LOCAL_PROPOSAL_PILOT","N"+str(level),{"strata":len(buckets),"singleton_edges":sum(len(z)==1 for z in buckets.values()),"largest_strata":sorted(map(len,buckets.values()),reverse=True)[:5],"attempts":samples,"outcomes":dict(c)})
    emit("ENVIRONMENT_MODULE_DISCOVERY",{x:importlib.util.find_spec(x) is not None for x in ("numpy","scipy","pandas","networkx","igraph","leidenalg","graph_tool","infomap","sklearn")})
    emit("P1D01_READONLY_PREFLIGHT_COMPLETE_NO_FITTING_OR_SPLITS")

if __name__=="__main__":main()
