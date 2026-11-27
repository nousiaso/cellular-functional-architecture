#!/usr/bin/env python3
"""Full-scale computational pilot on randomized native incidence.
No observed-network partition, holdout assignment, biological labels, saved
memberships, inference p-values, or certified null draws. N2/N3 are finite-sweep
workloads only. Composition is a separate fixed-topology workload.
"""
import argparse
from collections import Counter
import json
import math
import sys
import time
from p1d02_mixing_evidence_pilot import (
    ROOT, B01, C01, load, chain, keyfun, rows, sha, need, emit
)
RECEIPT_HASH="d6c476b68def958d8140d901ab896119630b9a00b401600f1a98b350d9d7eacc"
ROLES=("material_input","material_output","catalyst","positiveregulator",
       "negativeregulator","degradation_input","template","template_product",
       "component","memberPhysicalEntity")
OUTPUT={"material_output","template_product"}

def make_graph(gt,edge,n,comp):
    """Only active states in each workload, all operations. Distinct typed edges."""
    states=sorted({e for k,e,i in edge}|{x for k,a,b in comp for x in (a,b)})
    em={e:i for i,e in enumerate(states)}
    ne=len(states);g=gt.Graph(directed=True);g.add_vertex(ne+n)
    kind=g.new_vertex_property("int64_t");kind.a[:ne]=0;kind.a[ne:]=1
    ec=g.new_edge_property("int64_t")
    km={k:i for i,k in enumerate(ROLES)}
    data=[]
    for k,e,i in sorted(edge):
        a,b=(ne+i,em[e]) if k in OUTPUT else (em[e],ne+i)
        data.append((a,b,km[k]))
    for k,a,b in sorted(comp):
        data.append((em[a],em[b],km[k]))
    # Layer labels must be consecutive; template is empty in the native input.
    active=sorted({k for a,b,k in data})
    renum={k:i for i,k in enumerate(active)}
    data=[(a,b,renum[k]) for a,b,k in data]
    need(len(data)==len(set(data)) and all(a!=b for a,b,k in data),"typed graph duplicate or loop")
    g.add_edge_list(data,eprops=[ec])
    need(g.num_edges()==len(data),"graph edge count")
    expected=Counter(k for a,b,k in data);observed=Counter()
    for ed in g.edges():
        k=int(ec[ed]);a,b=int(kind[ed.source()]),int(kind[ed.target()])
        name=ROLES[active[k]]
        want=(0,0) if name in ("component","memberPhysicalEntity") else (1,0) if name in OUTPUT else (0,1)
        need((a,b)==want,"edge endpoint direction")
        observed[k]+=1
    need(expected==observed,"layer counts")
    return g,kind,ec,ne,{ROLES[active[k]]:v for k,v in sorted(observed.items())}

def fit_workload(gt,np,LayeredBlockState,minimize_blockmodel_dl,edge,n,comp,label,seed):
    gt.seed_rng(seed)
    t=time.monotonic()
    g,kind,ec,ne,counts=make_graph(gt,edge,n,comp)
    emit("WORKLOAD_GRAPH",label,{"vertices":g.num_vertices(),"active_states":ne,
        "operations":n,"edges":g.num_edges(),"layers":counts,"build_seconds":round(time.monotonic()-t,3)})
    args=dict(ec=ec,independent=True,deg_corr=True,clabel=kind,pclabel=kind)
    b0=g.new_vertex_property("int64_t");b0.a[:]=kind.a
    start=LayeredBlockState(g,b=b0,**args)
    initial=float(start.entropy());need(math.isfinite(initial),"initial description length")
    del start
    emit("WORKLOAD_OPTIMIZER_START",label,initial)
    t=time.monotonic()
    state=minimize_blockmodel_dl(g,state=LayeredBlockState,state_args=args)
    elapsed=time.monotonic()-t
    dl=float(state.entropy());b=np.asarray(state.get_blocks().a).copy()
    need(math.isfinite(dl) and len(b)==ne+n and not(set(b[:ne])&set(b[ne:])),"fitted state validity")
    report={"type_only_dl":initial,"fitted_dl":dl,"improves_type_only":dl<=initial+1e-7,
            "fit_seconds":round(elapsed,3),"operation_blocks":len(set(b[ne:])),
            "state_blocks":len(set(b[:ne])),"type_constraint_pass":True}
    emit("WORKLOAD_FIT_COMPLETE",label,report)
    # Do not persist or report individual memberships or use them to choose
    # a biological model. Scores across workloads have different observations.
    return report

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--workload",choices=("N2_roles","N3_roles","N2_composition","all"),default="all")
    opt=a.parse_args()
    receipt=ROOT/".conda/cfa_inference_v1.explicit.txt"
    need(sha(receipt)==RECEIPT_HASH,"environment receipt fingerprint")
    need((ROOT/".conda/cfa_inference_v1").resolve().as_posix()==__import__("pathlib").Path(sys.prefix).resolve().as_posix(),"use dedicated Python")
    emit("BACKEND_IMPORT_START")
    t=time.monotonic()
    import numpy as np
    import scipy
    import graph_tool as gt
    from graph_tool.inference import LayeredBlockState, minimize_blockmodel_dl
    need((gt.__version__,np.__version__,scipy.__version__)==("3.9","2.5.5","1.18.1"),"backend versions")
    emit("BACKEND_IMPORT_PASS",gt.__version__,np.__version__,scipy.__version__,round(time.monotonic()-t,3))
    R,edge,sp,rc,escope,root,direct=load()
    n=len(R);need(len(edge)==66440,"native incidence count")
    roots=sorted(set().union(*(set(v) for v in root.values())))
    bits={a:1<<i for i,a in enumerate(roots)}
    rootmask={i:sum(bits[a] for a in v) for i,v in root.items()}
    ds={i:set(v) for i,v in direct.items()}
    requested=("N2_roles","N3_roles","N2_composition") if opt.workload=="all" else (opt.workload,)
    random_edges={}
    for level in sorted({int(label[1]) for label in requested}):
        z,cp=chain(edge,keyfun(level,sp,rc,escope,root,direct),20261130,(20,),rootmask,ds)
        need(z!=edge and len(z)==len(edge),"randomized workload")
        random_edges[level]=z
        emit("RANDOMIZED_WORKLOAD",f"N{level}",cp[-1])
    comp=set()
    if "N2_composition" in requested:
        E={r["state_id"]:int(r["state_index"]) for r in rows(C01/"state_index.tsv")}
        nr=0
        for r in rows(B01/"composition.tsv"):
            k=r["relation"];a=E[r["owner_id"]];b=E[r["child_id"]]
            need(k in ("component","memberPhysicalEntity") and a!=b,"composition")
            comp.add((k,a,b));nr+=1
        need(nr==76196,"composition source rows")
        emit("COMPOSITION_FIXED_WORKLOAD",{"source_rows":nr,"distinct_typed_edges":len(comp)})
    reports={}
    for label in requested:
        reports[label]=fit_workload(gt,np,LayeredBlockState,minimize_blockmodel_dl,
            random_edges[int(label[1])],n,comp if label.endswith("composition") else set(),label,20261130)
    emit("P1D05_WORKLOAD_SUMMARY",reports)
    emit("P1D05_SCALE_PILOT_COMPLETE_NO_OBSERVED_FIT_OR_SAVED_PARTITIONS")
if __name__=="__main__":main()
