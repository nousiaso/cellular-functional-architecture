#!/usr/bin/env python3
"""Bounded single-vertex optimization pilot, full objective checked every sweep.
No observed fit, saved partitions, final null draws or biological claims.
Pilot schedules/start counts are not confirmatory settings.
"""
import argparse
import math
import random
import time
import numpy as np
import graph_tool as gt
from graph_tool.inference import LayeredBlockState
from p1d05_native_scale_pilot import make_graph, RECEIPT_HASH
from p1d05_entropy_audit import canonical
from p1d05_update_recovery import close
from p1d02_mixing_evidence_pilot import ROOT, B01, C01, load, chain, keyfun, rows, sha, need, emit

SEEDS=(20261130,20261201)
STARTS=(1,2,4,8,16)
SCHEDULE=((0.2,3),(0.5,3),(1.0,4),(np.inf,8))

def ari(a,b):
    a=canonical(a);b=canonical(b)
    _,cells=np.unique(np.stack((a,b),axis=1),axis=0,return_counts=True)
    ca=np.unique(a,return_counts=True)[1];cb=np.unique(b,return_counts=True)[1]
    pairs=lambda x:float(np.sum(x*(x-1)/2))
    t=pairs(cells);x=pairs(ca);y=pairs(cb);n=len(a)*(len(a)-1)/2
    if n==0:return 1.0
    expected=x*y/n;den=(x+y)/2-expected
    return (t-expected)/den if den else float(np.array_equal(a,b))

def synthetic(m,p,q):
    rng=random.Random(20261127);n=2*m
    g=gt.Graph(directed=True);g.add_vertex(2*n)
    ec=g.new_edge_property("int64_t");kind=g.new_vertex_property("int64_t")
    kind.a[:n]=0;kind.a[n:]=1
    data=[]
    for e in range(n):
        for i in range(n):
            for k in range(3):
                match=(e//m==i//m) if k!=1 else (e//m!=i//m)
                if rng.random()<(p if match else q):
                    a,b=(n+i,e) if k==1 else (e,n+i)
                    data.append((a,b,k))
    g.add_edge_list(data,eprops=[ec])
    need(g.num_edges()==len(data)==len(set(data)),"synthetic edges")
    truth=np.repeat(np.arange(4,dtype="int64"),m)
    return g,kind,ec,n,truth

def search(label,g,kind,ec,ne,truth=None):
    t=time.monotonic();orig=kind.a.copy();oe=ec.a.copy();ge=g.get_edges().copy()
    n=g.num_vertices();nr=n-ne
    need(ne>0 and nr>0,"nonempty vertex types")
    def vp(a):
        p=g.new_vertex_property("int64_t");p.a[:]=a;return p
    def make(a):
        cl=vp(orig)
        return LayeredBlockState(g,b=vp(canonical(a)),ec=ec,independent=True,
                                 deg_corr=True,clabel=cl,pclabel=cl)
    def score(a):
        s=make(a);v=float(s.entropy())
        need(math.isfinite(v),"nonfinite fresh score")
        return v
    def describe(a):
        out={"state_blocks":len(set(a[:ne])),"operation_blocks":len(set(a[ne:]))}
        if truth is not None:
            out.update(state_ARI=ari(a[:ne],truth[:ne]),operation_ARI=ari(a[ne:],truth[ne:]))
        return out
    baseline=score(orig);best_score=baseline;best=orig.copy()
    emit("SEARCH_START",label,{"vertices":n,"edges":g.num_edges(),"baseline":baseline,
        "planted_score":score(truth) if truth is not None else None,
        "start_blocks_per_type":STARTS,"seeds":SEEDS,"schedule":[(str(b),c) for b,c in SCHEDULE]})
    checks=0;total_moves=0
    for k in STARTS:
        for seed in SEEDS:
            gt.seed_rng(seed);rng=np.random.default_rng(seed)
            a=np.concatenate((rng.integers(0,k,ne),rng.integers(k,2*k,nr))).astype("int64")
            s=make(a);initial=float(s.entropy());local_best=initial
            if initial<best_score:
                best_score=initial;best=canonical(s.get_blocks().a.copy())
            moves=0;sweeps=0;t1=time.monotonic()
            for beta,ns in SCHEDULE:
                for i in range(ns):
                    before=float(s.entropy());b0=canonical(s.get_blocks().a.copy())
                    ret=s.mcmc_sweep(beta=beta,niter=1)
                    after=float(s.entropy());b=canonical(s.get_blocks().a.copy())
                    fresh=score(b);delta=after-before
                    ok=(math.isfinite(after) and close(after,fresh) and close(delta,ret[0])
                        and (not np.array_equal(b0,b) or close(delta,0.0))
                        and not(set(b[:ne])&set(b[ne:]))
                        and np.array_equal(kind.a,orig) and np.array_equal(ec.a,oe)
                        and np.array_equal(g.get_edges(),ge))
                    if math.isinf(beta) and delta>1e-6:ok=False
                    if not ok:
                        emit("SEARCH_CONSISTENCY_FAIL",label,k,seed,str(beta),i,
                             {"before":before,"after":after,"fresh":fresh,"ret":ret})
                        raise RuntimeError("single-vertex objective consistency gate")
                    checks+=1;sweeps+=1;moves+=int(ret[2]);local_best=min(local_best,fresh)
                    if fresh<best_score:
                        best_score=fresh;best=b.copy()
                emit("SEARCH_PHASE",label,k,seed,str(beta),{"sweeps":sweeps,
                     "DL":after,"local_best":local_best,"moves":moves,**describe(b)})
            total_moves+=moves
            emit("SEARCH_START_COMPLETE",label,k,seed,{"initial":initial,"final":after,
                 "local_best":local_best,"seconds":round(time.monotonic()-t1,3)})
            del s
    need(close(score(best),best_score) and best_score<=baseline+1e-6,"best candidate validity")
    result={"baseline":baseline,"best_DL":best_score,"gain":baseline-best_score,
            "checks":checks,"reported_moves":total_moves,"seconds":round(time.monotonic()-t,3),
            **describe(best)}
    if truth is not None:
        result["planted_DL"]=score(truth)
    emit("SEARCH_RESULT",label,result)
    return result

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--workload",choices=("all","synthetic","N2_roles","N3_roles","N2_composition"),default="all")
    opt=a.parse_args()
    need(sha(ROOT/".conda/cfa_inference_v1.explicit.txt")==RECEIPT_HASH,"receipt")
    need(gt.__version__=="3.9" and np.__version__=="2.5.5","backend versions")
    results={}
    if opt.workload in ("all","synthetic"):
        g,kind,ec,ne,truth=synthetic(60,.23,.012)
        need(g.num_edges()==5248,"planted smoke graph fingerprint count")
        r=search("planted_dense",g,kind,ec,ne,truth)
        results["planted_dense"]=r
        # Favorable known positive control; require both recovery and improvement.
        need(r["gain"]>0 and min(r["state_ARI"],r["operation_ARI"])>=.95,
             "dense planted recovery gate")
        del g,kind,ec
        g,kind,ec,ne,truth=synthetic(300,.012,.0006)
        results["planted_sparser"]=search("planted_sparser",g,kind,ec,ne,truth)
        del g,kind,ec
    if opt.workload!="synthetic":
        R,edge,sp,rc,escope,root,direct=load()
        roots=sorted(set().union(*(set(v) for v in root.values())))
        bits={a:1<<i for i,a in enumerate(roots)}
        rootmask={i:sum(bits[a] for a in v) for i,v in root.items()}
        ds={i:set(v) for i,v in direct.items()}
        labels=("N2_roles","N3_roles","N2_composition") if opt.workload=="all" else (opt.workload,)
        for label in labels:
            level=int(label[1])
            z,cp=chain(edge,keyfun(level,sp,rc,escope,root,direct),20261130,(20,),rootmask,ds)
            emit("RANDOMIZED_WORKLOAD",label,cp[-1])
            comp=set()
            if label.endswith("composition"):
                E={r["state_id"]:int(r["state_index"]) for r in rows(C01/"state_index.tsv")}
                for r in rows(B01/"composition.tsv"):
                    comp.add((r["relation"],E[r["owner_id"]],E[r["child_id"]]))
                need(len(comp)==76196,"composition count")
            g,kind,ec,ne,counts=make_graph(gt,z,len(R),comp)
            results[label]=search(label,g,kind,ec,ne)
            del g,kind,ec
    emit("P1D05_SEARCH_PILOT_SUMMARY",results)
    emit("P1D05_SINGLE_VERTEX_PILOT_COMPLETE_NO_OBSERVED_FIT")

if __name__=="__main__":main()
