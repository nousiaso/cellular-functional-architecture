#!/usr/bin/env python3
"""Synthetic-only initialization/recovery diagnostic for the unchanged full objective.
Truth-based starts are oracle diagnostics, never usable for observed inference.
Two groups per type are known control design, not an inferred biological block count.
"""
import math
import time
import numpy as np
import graph_tool as gt
from graph_tool.inference import LayeredBlockState
from p1d05_single_vertex_search import synthetic, ari, SEEDS, SCHEDULE
from p1d05_native_scale_pilot import RECEIPT_HASH
from p1d05_entropy_audit import canonical
from p1d05_update_recovery import close
from p1d02_mixing_evidence_pilot import ROOT, sha, need, emit

def two_means(x):
    # Exact 1D two-means cut, no forced equal sizes.
    order=np.argsort(x,kind="stable");v=x[order];n=len(v)
    s=np.cumsum(v);q=np.cumsum(v*v);k=np.arange(1,n)
    cost=q[:-1]-s[:-1]**2/k+q[-1]-q[:-1]-(s[-1]-s[:-1])**2/(n-k)
    cost[v[:-1]==v[1:]]=np.inf
    need(np.isfinite(cost).any(),"nonconstant spectral coordinate")
    cut=int(np.argmin(cost))+1
    b=np.ones(n,dtype="int64");b[order[:cut]]=0
    return b

def spectral(g,ec,ne):
    # Synthetic-only bipartite role matrices; do not collapse role channels.
    nr=g.num_vertices()-ne
    mats=np.zeros((3,ne,nr),dtype=float)
    for e in g.edges():
        u,v,k=int(e.source()),int(e.target()),int(ec[e])
        need(0<=k<3,"synthetic role")
        if k==1:
            need(u>=ne and v<ne,"reverse role orientation");u,v=v,u
        else:
            need(u<ne and v>=ne,"forward role orientation")
        mats[k,u,v-ne]+=1
    need(mats.sum()==g.num_edges(),"spectral incidence count")
    # Center each role separately; leading eigenvector of concatenated
    # centered features' Gram matrix for each vertex type independently.
    gs=np.zeros((ne,ne));go=np.zeros((nr,nr))
    for a in mats:
        x=a-a.mean(axis=0,keepdims=True)-a.mean(axis=1,keepdims=True)+a.mean()
        gs+=x@x.T;go+=x.T@x
    def split(gram):
        eig,vec=np.linalg.eigh(gram)
        need(np.isfinite(eig).all() and eig[-1]>0,"spectral eigenvalues")
        return two_means(vec[:,-1]),eig[-3:].tolist()
    bs,es=split(gs);bo,eo=split(go)
    return np.concatenate((bs,bo+2)),{"state_top_eigenvalues":es,"operation_top_eigenvalues":eo}

def run(label,g,kind,ec,ne,truth):
    t=time.monotonic();orig=kind.a.copy();oe=ec.a.copy();ge=g.get_edges().copy()
    def vp(a):
        p=g.new_vertex_property("int64_t");p.a[:]=a;return p
    def make(a):
        cl=vp(orig)
        return LayeredBlockState(g,b=vp(canonical(a)),ec=ec,independent=True,
                                 deg_corr=True,clabel=cl,pclabel=cl)
    def score(a):
        val=float(make(a).entropy());need(math.isfinite(val),"fresh score")
        return val
    def describe(a):
        return {"state_ARI":ari(a[:ne],truth[:ne]),"operation_ARI":ari(a[ne:],truth[ne:]),
                "state_blocks":len(set(a[:ne])),"operation_blocks":len(set(a[ne:]))}
    base=score(orig);planted=score(truth)
    spec,eigs=spectral(g,ec,ne)
    emit("INIT_CONTROL",label,{"baseline_DL":base,"planted_DL":planted,
        "spectral_DL":score(spec),**describe(spec),**eigs})
    checks=0
    for seed in SEEDS:
        rng=np.random.default_rng(seed)
        starts={"spectral":spec.copy(),"oracle_exact":truth.copy()}
        for frac in (.10,.25,.50):
            a=truth.copy()
            for lo,hi in ((0,ne),(ne,len(a))):
                idx=rng.choice(np.arange(lo,hi),size=round(frac*(hi-lo)),replace=False)
                a[idx]^=1
            starts[f"oracle_flip_{frac:.2f}"]=a
        for name,a in starts.items():
            for schedule_name,schedule in (("greedy",((np.inf,18),)),("anneal",SCHEDULE)):
                gt.seed_rng(seed);s=make(a);initial=score(a);best=initial;best_b=canonical(a)
                moves=0
                emit("INIT_START",label,name,seed,schedule_name,{"DL":initial,**describe(a)})
                for beta,ns in schedule:
                    for _ in range(ns):
                        before=float(s.entropy());b0=canonical(s.get_blocks().a.copy())
                        ret=s.mcmc_sweep(beta=beta,niter=1)
                        after=float(s.entropy());b=canonical(s.get_blocks().a.copy())
                        fresh=score(b);delta=after-before
                        need(math.isfinite(after) and close(after,fresh) and close(delta,ret[0])
                            and (not np.array_equal(b0,b) or close(delta,0.0))
                            and not(set(b[:ne])&set(b[ne:]))
                            and np.array_equal(kind.a,orig) and np.array_equal(ec.a,oe)
                            and np.array_equal(g.get_edges(),ge)
                            and (not math.isinf(beta) or delta<=1e-6),
                            "single-vertex consistency gate")
                        checks+=1;moves+=int(ret[2])
                        if fresh<best:best=fresh;best_b=b.copy()
                need(close(score(best_b),best),"best score")
                emit("INIT_RESULT",label,name,seed,schedule_name,
                     {"initial_DL":initial,"best_DL":best,"final_DL":fresh,
                      "gap_to_planted":best-planted,"moves":moves,
                      "best":describe(best_b),"final":describe(b)})
                del s
    emit("INIT_CONTROL_COMPLETE",label,{"checks":checks,"seconds":round(time.monotonic()-t,3)})

def main():
    need(sha(ROOT/".conda/cfa_inference_v1.explicit.txt")==RECEIPT_HASH,"receipt")
    need(gt.__version__=="3.9" and np.__version__=="2.5.5","backend versions")
    for label,m,p,q,edges in (("dense",60,.23,.012,5248),("sparser",300,.012,.0006,6756)):
        g,kind,ec,ne,truth=synthetic(m,p,q)
        need(g.num_edges()==edges,"control edge count")
        run(label,g,kind,ec,ne,truth)
    emit("P1D05_INITIALIZATION_DIAGNOSTIC_COMPLETE_NO_OBSERVED_FIT")

if __name__=="__main__":main()
