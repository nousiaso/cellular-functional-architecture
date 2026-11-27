#!/usr/bin/env python3
"""Pilot only: finite-sweep native-incidence switch chains and evidence split feasibility.
No inferential p-values, no persisted randomized graphs or chosen split.
"""
from collections import Counter, defaultdict
import itertools
import json
import random
import time
from p1d01_network_feasibility import (
    ROOT, B01, B03, C01, C02, C01_HASH, C02_HASH,
    verify, rows, sha, need, emit, arr, DSU
)

def load():
    lh,bh=verify()
    p=C02/"manifest.json";need(sha(p)==C02_HASH,"C02 manifest")
    m=json.loads(p.read_text())
    need(m["P1C01_manifest_sha256"]==C01_HASH and m["source_lock_sha256"]==lh and m["harmonization_lock_sha256"]==bh,"C02 lineage")
    for name,x in m["files"].items():
        need(sha(C02/name)==x["sha256"] and (C02/name).stat().st_size==x["bytes"],"C02 file "+name)
    E={r["state_id"]:int(r["state_index"]) for r in rows(C01/"state_index.tsv")}
    R={r["operation_id"]:int(r["operation_index"]) for r in rows(C01/"operation_index.tsv")}
    n=len(R);need((len(E),n)==(61418,16187),"universe")
    rc={R[r["operation_id"]]:r["class"] for r in rows(C01/"operation_index.tsv")}
    scope={E[r["state_id"]]:r["composition_feasibility"] for r in rows(B03/"state_scope.tsv")}
    escope={R[r["operation_id"]]:r["candidate_scope_stratum"] for r in rows(B03/"event_scope.tsv")}
    sp={}
    for r in rows(B01/"nodes.tsv"):
        if r["node_id"] in E:
            e=E[r["node_id"]];x=json.loads(r["payload_json"])
            sp[e]=(r["class"],scope[e],tuple(sorted(arr(x,"location_ids_json"))))
    root={};direct={}
    for r in rows(B03/"event_pathway_benchmark.tsv"):
        i=R[r["operation_id"]];root[i]=tuple(sorted(arr(r,"native_root_json")));direct[i]=tuple(sorted(arr(r,"native_direct_json")))
    edge=set()
    for r in rows(C01/"typed_incidence.tsv"):
        z=(r["role"],int(r["state_index"]),int(r["operation_index"]))
        need(z not in edge and r["binary_value"]=="1","native incidence")
        edge.add(z)
    return R,edge,sp,rc,escope,root,direct

def evidence(R,edge,root,direct):
    n=len(R);pub=defaultdict(set)
    for r in rows(B03/"bp_publication_links.tsv"):
        i=R[r["owner_id"]]
        for a in arr(r,"namespaces_json"):
            for b in arr(r,"identifiers_json"):
                if isinstance(a,str) and a and isinstance(b,str) and b:pub[a.lower(),b].add(i)
    known=set().union(*pub.values());d=DSU(n)
    for es in pub.values():
        a=min(es)
        for b in es:d.union(a,b)
    groups=defaultdict(set)
    for i in known:groups[d.find(i)].add(i)
    g=sorted(groups.values(),key=lambda z:(-len(z),min(z)))
    need(len(g)==4806 and len(g[0])==4962 and len(known)==15739,"evidence components")
    rootops=defaultdict(set)
    for i,rs in root.items():
        for a in rs:rootops[a].add(i)
    paths={r["pathway_id"]:json.loads(r["name_json"]) for r in rows(B03/"pathways.tsv")}
    emit("EVIDENCE_GIANT_PROFILE",{"operations":len(g[0]),"unknown_publication_operations":n-len(known),
        "root_counts":[[a,paths[a],len(os&g[0]),len(os)] for a,os in sorted(rootops.items())]})
    byop=defaultdict(set)
    for k,e,i in edge:byop[i].add((k,e))
    # Illustrative capacity check: giant component stays in training,
    # unknown publication operations are quarantined, no final assignment.
    # Holding out an entire operation leaves no native incidences for it.
    # A separate frozen mapping rule would be required for new-node inference.
    target=round(.2*len(known));all_direct=set().union(*(set(direct[i]) for i in known))
    for seed in (20261127,20261128,20261129):
        rng=random.Random(seed);pool=list(g[1:]);rng.shuffle(pool)
        test=set();chosen=0
        for z in pool:
            if len(test)>=target:break
            if abs(target-len(test)-len(z))<=abs(target-len(test)):
                test.update(z);chosen+=1
        train=known-test
        need(not any(es&train and es&test for es in pub.values()),"publication crossing")
        tr_features=set().union(*(byop[i] for i in train))
        te_features=set().union(*(byop[i] for i in test))
        tr_paths=set().union(*(set(direct[i]) for i in train))
        te_paths=set().union(*(set(direct[i]) for i in test))
        emit("ILLUSTRATIVE_GROUP_HOLDOUT",{"seed":seed,"train_operations":len(train),"test_operations":len(test),"test_groups":chosen,
          "quarantine_unknown":n-len(known),"shared_publication_IDs":0,
          "test_role_state_features_seen_train":len(tr_features&te_features),"test_role_state_features_total":len(te_features),
          "test_direct_pathways_seen_train":len(tr_paths&te_paths),"test_direct_pathways_total":len(te_paths),
          "all_known_direct_pathways":len(all_direct),"test_roots":[[a,paths[a],len(os&test),len(os&train),len(os)] for a,os in sorted(rootops.items())],
          "test_operations_with_no_role_state_feature_in_train":sum(not(byop[i]&tr_features) for i in test)})
    emit("EVIDENCE_PILOT_LIMIT","whole-operation split; not known-node prediction; observed publication grouping does not exclude unobserved dependence or molecular near-duplicates")

def keyfun(level,sp,rc,escope,root,direct):
    def key(k,e,i):
        z=(k,)
        if level>=1:z+=(sp[e],rc[i],escope[i])
        if level>=2:z+=(root[i],)
        if level>=3:z+=(direct[i],)
        return z
    return key

def margins(edge,key):
    return (Counter((key(k,e,i),e) for k,e,i in edge),
            Counter((key(k,e,i),i) for k,e,i in edge))

def diagnostics(edge,base,rootmask,direct):
    out=defaultdict(set);inn=defaultdict(set)
    per=Counter(k for k,e,i in edge);overlap=Counter(k for k,e,i in edge&base)
    for k,e,i in edge:
        if k=="material_output":out[e].add(i)
        elif k=="material_input":inn[e].add(i)
    count=same=dsame=0
    for e in out.keys()&inn.keys():
        for i in out[e]:
            for j in inn[e]:
                if i==j:continue
                count+=1;same+=bool(rootmask[i]&rootmask[j]);dsame+=bool(direct[i]&direct[j])
    return {"original_edge_retention":len(edge&base)/len(base),
            "per_role_retention":{k:overlap[k]/v for k,v in sorted(per.items())},
            "flow_witnesses":count,"flow_shared_root_witnesses":same,"flow_shared_direct_witnesses":dsame,
            "flow_cross_root_fraction":(count-same)/count if count else None}

def chain(base,key,seed,sweeps,rootmask,direct):
    rng=random.Random(seed);slots=sorted(base);current=set(base);buckets=defaultdict(list)
    for j,(k,e,i) in enumerate(slots):buckets[key(k,e,i)].append(j)
    keys=[key(*z) for z in slots];M=len(slots);original_margins=margins(base,key)
    count=Counter();checkpoints=[];t0=time.monotonic()
    for s in range(1,sweeps[-1]+1):
        for _ in range(M):
            a=rng.randrange(M);pool=buckets[keys[a]]
            if len(pool)==1:count["singleton"]+=1;continue
            b=rng.choice(pool)
            k,e,i=slots[a];l,f,j=slots[b]
            if e==f or i==j:count["same_endpoint"]+=1;continue
            x=(k,e,j);y=(l,f,i)
            if x in current or y in current:count["existing_cross"]+=1;continue
            current.remove(slots[a]);current.remove(slots[b]);current.update((x,y));slots[a]=x;slots[b]=y
            count["accepted"]+=1
        if s in sweeps:
            need(len(current)==M and len(set(slots))==M and all(key(*z)==q for z,q in zip(slots,keys)),"swap constraints")
            need(margins(current,key)==original_margins,"stratified degree invariants")
            info=diagnostics(current,base,rootmask,direct)
            info.update({"sweep":s,"attempts":s*M,"accepted":count["accepted"],"seconds":round(time.monotonic()-t0,3)})
            checkpoints.append(info)
    return current,checkpoints

def main():
    R,edge,sp,rc,escope,root,direct=load()
    evidence(R,edge,root,direct)
    roots=sorted(set().union(*(set(v) for v in root.values())))
    bits={a:1<<i for i,a in enumerate(roots)}
    rootmask={i:sum(bits[a] for a in v) for i,v in root.items()}
    ds={i:set(v) for i,v in direct.items()}
    emit("OBSERVED_NATIVE",diagnostics(edge,edge,rootmask,ds))
    for level in range(4):
        key=keyfun(level,sp,rc,escope,root,direct)
        buckets=Counter(key(*z) for z in edge)
        emit("CHAIN_LEVEL","N"+str(level),{"edges":len(edge),"strata":len(buckets),"singleton_edges":sum(v==1 for v in buckets.values())})
        final=[]
        for seed in (20261127,20261128,20261129):
            z,checkpoints=chain(edge,key,seed,(1,5,10,20),rootmask,ds)
            final.append(z)
            for x in checkpoints:emit("CHAIN_CHECKPOINT","N"+str(level),seed,x)
        emit("CHAIN_BETWEEN_FINAL_RETENTION","N"+str(level),[len(a&b)/len(edge) for a,b in itertools.combinations(final,2)])
    emit("P1D02_PILOT_COMPLETE_NO_INFERENCE_OR_SAVED_SPLITS")

if __name__=="__main__":main()
