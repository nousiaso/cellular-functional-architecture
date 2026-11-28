#!/usr/bin/env python3
"""Unknown-block-count synthetic benchmark. No observed fit or final settings.
Sparse implicit centered role Gram operators, no dense n-by-n projection.
Truth is used only for simulation and evaluation, never initialization/selection.
"""
import math
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import LinearOperator, eigsh
import graph_tool as gt
from graph_tool.inference import LayeredBlockState
from p1d05_single_vertex_search import synthetic, ari, SCHEDULE
from p1d05_native_scale_pilot import RECEIPT_HASH
from p1d05_entropy_audit import canonical
from p1d05_update_recovery import close
from p1d02_mixing_evidence_pilot import ROOT, sha, need, emit

# Provisional common search grid, not truth-dependent or final.
PAIRS=((1,1),(2,2),(3,3),(4,4),(6,6),(2,3),(3,2),(3,4),(4,3))
GRAPH_SEEDS=(20261202,20261203)
CASES=(("two_degree12",2,2,12,False,True),
       ("two_degree3",2,2,3,False,True),
       ("three_degree6",3,3,6,False,True),
       ("asymmetric_degree12",3,4,12,True,True),
       ("asymmetric_degree3",3,4,3,True,True),
       ("heterogeneous_two_degree3",2,2,3,True,True),
       ("no_signal_degree3",2,2,3,True,False))

def generate(ks,ko,degree,heterogeneous,signal,seed):
    rng=np.random.default_rng(seed);ne=600;no=400
    def groups(n,k):
        a=np.arange(n)%k
        if heterogeneous:
            a=rng.choice(k,n,p=np.arange(1,k+1)/sum(range(1,k+1)))
        rng.shuffle(a)
        need(len(set(a))==k,"planted group coverage")
        return a
    bs,bo=groups(ne,ks),groups(no,ko)
    ws=rng.lognormal(0,.7,ne).clip(.2,4) if heterogeneous else np.ones(ne)
    wo=rng.lognormal(0,.7,no).clip(.2,4) if heterogeneous else np.ones(no)
    g=gt.Graph(directed=True);g.add_vertex(ne+no)
    ec=g.new_edge_property("int64_t");kind=g.new_vertex_property("int64_t")
    kind.a[:ne]=0;kind.a[ne:]=1
    edges=[]
    for r in range(3):
        affinity=np.where((bs[:,None]+r)%ko==bo[None,:],20.,1.) if signal else np.ones((ne,no))
        # Additional role-specific mapping lets unequal groups have distinct profiles.
        if signal and ks!=ko:
            affinity+=8.*(((bo[None,:]+r)%ks)==bs[:,None])
        w=affinity*ws[:,None]*wo[None,:]
        p=w*(degree*(ne+no)/6)/w.sum()
        need(p.max()<=1,"invalid Bernoulli probability")
        u,v=np.nonzero(rng.random(p.shape)<p)
        edges.extend((int(ne+j),int(i),r) if r==1 else (int(i),int(ne+j),r)
                     for i,j in zip(u,v))
    g.add_edge_list(edges,eprops=[ec])
    need(g.num_edges()==len(edges)==len(set(edges)),"synthetic incidence")
    truth=np.concatenate((bs,bo+ks)).astype("int64") if signal else None
    return g,kind,ec,ne,truth

def embedding(g,ec,ne):
    no=g.num_vertices()-ne;items=[([],[]) for _ in range(3)]
    for e in g.edges():
        u,v,r=int(e.source()),int(e.target()),int(ec[e])
        need(0<=r<3,"role index")
        if r==1:u,v=v,u
        need(u<ne and v>=ne,"role direction")
        items[r][0].append(u);items[r][1].append(v-ne)
    mats=[sp.csr_matrix((np.ones(len(u)),(u,v)),shape=(ne,no)) for u,v in items]
    need(sum(a.sum() for a in mats)==g.num_edges(),"role matrix count")
    def center(x):return x-x.mean()
    def eig(matrices,n):
        def mv(v):
            v=center(v);out=np.zeros(n)
            for a in matrices:
                out+=center(a@center(a.T@v))
            return out
        op=LinearOperator((n,n),matvec=mv,dtype=np.float64)
        rng=np.random.default_rng(20261204)
        vals,vec=eigsh(op,k=max(max(x) for x in PAIRS)-1,which="LA",
                       v0=rng.normal(size=n),tol=1e-7)
        order=np.argsort(vals)[::-1];vals,vec=vals[order],vec[:,order]
        need(np.isfinite(vec).all() and np.isfinite(vals).all(),"spectral finite")
        residual=max(np.linalg.norm(mv(vec[:,j])-vals[j]*vec[:,j])/max(1,abs(vals[j]))
                     for j in range(len(vals)))
        need(residual<1e-5,"spectral residual")
        return vec,vals.tolist(),float(residual)
    a,ea,ra=eig(mats,ne);b,eb,rb=eig([x.T.tocsr() for x in mats],no)
    return a,b,{"state_eigenvalues":ea,"operation_eigenvalues":eb,
               "state_residual":ra,"operation_residual":rb}

def cluster(x,k,rng):
    n=len(x)
    if k==1:return np.zeros(n,dtype="int64")
    x=x[:,:k-1]
    centers=[x[rng.integers(n)]]
    for _ in range(1,k):
        d=np.min(((x[:,None,:]-np.array(centers)[None,:,:])**2).sum(axis=2),axis=1)
        need(d.sum()>0,"insufficient distinct embedding rows")
        centers.append(x[rng.choice(n,p=d/d.sum())])
    centers=np.array(centers);prev=None
    for _ in range(100):
        d=((x[:,None,:]-centers[None,:,:])**2).sum(axis=2)
        b=d.argmin(axis=1);counts=np.bincount(b,minlength=k)
        for j in np.flatnonzero(counts==0):
            candidates=np.flatnonzero(counts[b]>1)
            i=candidates[np.argmax(d[candidates,b[candidates]])]
            counts[b[i]]-=1;b[i]=j;counts[j]+=1
        if prev is not None and np.array_equal(b,prev):break
        prev=b.copy();centers=np.array([x[b==j].mean(axis=0) for j in range(k)])
    need(len(set(b))==k,"cluster occupancy")
    return b.astype("int64")

def run(label,seed,g,kind,ec,ne,truth):
    t=time.monotonic();orig=kind.a.copy();oe=ec.a.copy();ge=g.get_edges().copy()
    n=g.num_vertices();no=n-ne
    active=np.array([v.in_degree()+v.out_degree()>0 for v in g.vertices()])
    def vp(a):
        p=g.new_vertex_property("int64_t");p.a[:]=a;return p
    def make(a):
        cl=vp(orig)
        return LayeredBlockState(g,b=vp(canonical(a)),ec=ec,independent=True,
                                 deg_corr=True,clabel=cl,pclabel=cl)
    def score(a):
        val=float(make(a).entropy());need(math.isfinite(val),"fresh DL");return val
    def desc(a):
        d={"state_blocks":len(set(a[:ne])),"operation_blocks":len(set(a[ne:]))}
        if truth is not None:
            d.update(state_ARI=ari(a[:ne],truth[:ne]),operation_ARI=ari(a[ne:],truth[ne:]),
                     active_state_ARI=ari(a[:ne][active[:ne]],truth[:ne][active[:ne]]),
                     active_operation_ARI=ari(a[ne:][active[ne:]],truth[ne:][active[ne:]]))
        return d
    base=score(orig);planted=score(truth) if truth is not None else None
    xs,xo,meta=embedding(g,ec,ne)
    emit("GRID_CONTROL",label,seed,{"vertices":n,"edges":g.num_edges(),
        "mean_degree":2*g.num_edges()/n,"isolates_state":int((~active[:ne]).sum()),
        "isolates_operation":int((~active[ne:]).sum()),"baseline_DL":base,
        "planted_DL":planted,"pairs":PAIRS,**meta})
    checks=0;overall={};schedules=(("greedy",((np.inf,18),)),("anneal",SCHEDULE))
    for method in ("random","spectral"):
        for schedule_name,schedule in schedules:
            overall[(method,schedule_name)]=(base,orig.copy(),(1,1))
        for ks,ko in PAIRS:
            rng=np.random.default_rng(seed+100*ks+ko)
            if method=="random":
                a=np.concatenate((rng.integers(ks,size=ne),rng.integers(ko,size=no)+ks))
            else:a=np.concatenate((cluster(xs,ks,rng),cluster(xo,ko,rng)+ks))
            for schedule_name,schedule in schedules:
                gt.seed_rng(seed);s=make(a);initial=score(a);best=initial;best_b=canonical(a)
                moves=0
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
                            and (not math.isinf(beta) or delta<=1e-6),"objective gate")
                        checks+=1;moves+=int(ret[2])
                        if fresh<best:best=fresh;best_b=b.copy()
                need(close(score(best_b),best),"best DL")
                key=(method,schedule_name)
                if best<overall[key][0]:overall[key]=(best,best_b.copy(),(ks,ko))
                emit("GRID_START_RESULT",label,seed,method,schedule_name,ks,ko,
                     {"initial_DL":initial,"best_DL":best,"final_DL":fresh,
                      "moves":moves,"initial":desc(a),"best":desc(best_b)})
                del s
    for (method,schedule_name),(dl,b,pair) in overall.items():
        emit("GRID_BEST",label,seed,method,schedule_name,{"start_pair":pair,
            "best_DL":dl,"gain":base-dl,"gap_to_planted":None if planted is None else dl-planted,
            **desc(b)})
    emit("GRID_CONTROL_COMPLETE",label,seed,{"checks":checks,"seconds":round(time.monotonic()-t,3)})

def main():
    need(sha(ROOT/".conda/cfa_inference_v1.explicit.txt")==RECEIPT_HASH,"receipt")
    need(gt.__version__=="3.9" and np.__version__=="2.5.5","backend")
    g,kind,ec,ne,truth=synthetic(300,.012,.0006)
    need(g.num_edges()==6756,"prior control count")
    run("previous_sparse",20261202,g,kind,ec,ne,truth)
    for label,ks,ko,degree,heterogeneous,signal in CASES:
        for seed in GRAPH_SEEDS:
            g,kind,ec,ne,truth=generate(ks,ko,degree,heterogeneous,signal,seed)
            run(label,seed,g,kind,ec,ne,truth)
    emit("P1D05_BLOCK_GRID_BENCHMARK_COMPLETE_NO_OBSERVED_FIT")

if __name__=="__main__":main()
