#!/usr/bin/env python3
"""Read-only diagnosis of layered objective drift on previously audited N2 workload."""
import inspect
import math
import time
import hashlib
import numpy as np
import graph_tool as gt
from graph_tool.inference import LayeredBlockState
from p1d05_native_scale_pilot import make_graph, RECEIPT_HASH
from p1d02_mixing_evidence_pilot import ROOT, load, chain, keyfun, sha, need, emit

def canonical(a):
    # Labels are immaterial; preserve exact membership equivalence.
    _,first,inv=np.unique(a,return_index=True,return_inverse=True)
    return np.argsort(np.argsort(first))[inv.reshape(-1)]

def main():
    need(sha(ROOT/".conda/cfa_inference_v1.explicit.txt")==RECEIPT_HASH,"receipt")
    need(gt.__version__=="3.9","backend")
    R,edge,sp,rc,escope,root,direct=load()
    roots=sorted(set().union(*(set(v) for v in root.values())))
    bits={a:1<<i for i,a in enumerate(roots)}
    rootmask={i:sum(bits[a] for a in v) for i,v in root.items()}
    ds={i:set(v) for i,v in direct.items()}
    z,cp=chain(edge,keyfun(2,sp,rc,escope,root,direct),20261130,(20,),rootmask,ds)
    emit("N2_WORKLOAD",cp[-1])
    g,kind,ec,ne,counts=make_graph(gt,z,len(R),set())
    original_kind=kind.a.copy()
    original_ec=ec.a.copy()
    original_edges=g.get_edges().copy()
    args=dict(ec=ec,independent=True,deg_corr=True,clabel=kind,pclabel=kind)
    def make(a,**kw):
        b=g.new_vertex_property("int64_t");b.a[:]=a
        return LayeredBlockState(g,b=b,**args,**kw)
    def properties(s):
        out={"repr":str(s),"occupied_labels":len(np.unique(s.get_blocks().a)),
             "min_label":int(np.min(s.get_blocks().a)),"max_label":int(np.max(s.get_blocks().a))}
        for k in ("B","get_B","get_N","get_E","get_entropy_args"):
            try:
                x=getattr(s,k);x=x() if callable(x) else x
                out[k]=str(x)
            except Exception as e:out[k]=str(e)
        try:out["block_graph_vertices"]=s.bg.num_vertices()
        except Exception:pass
        return out
    options={"full":{},"likelihood_only":{"dl":False},"without_partition_dl":{"partition_dl":False},
             "without_degree_dl":{"degree_dl":False},"without_edges_dl":{"edges_dl":False},
             "without_layer_membership":{"lmemb":False},"without_adjacency":{"adjacency":False}}
    def scores(s):
        out={}
        for k,kw in options.items():
            try:out[k]=float(s.entropy(**kw))
            except Exception as e:out[k]=str(e)
        return out
    gt.seed_rng(20261130)
    s=make(original_kind)
    emit("BASELINE_PROPERTIES",properties(s))
    emit("BASELINE_ENTROPY",scores(s))
    emit("ENTROPY_SIGNATURE",str(inspect.signature(s.entropy)))
    b=canonical(s.get_blocks().a.copy())
    for i in range(5):
        before=float(s.entropy())
        ret=s.multiflip_mcmc_sweep(beta=np.inf,niter=1)
        after=float(s.entropy())
        a=s.get_blocks().a.copy()
        need(np.array_equal(kind.a,original_kind) and np.array_equal(ec.a,original_ec)
             and np.array_equal(g.get_edges(),original_edges),"graph or input properties changed")
        fresh=make(a)
        compact=make(canonical(a))
        emit("SWEEP",i+1,{"reported":ret,"measured_delta":after-before,
            "same_partition":bool(np.array_equal(canonical(a),b)),
            "fresh_dl":float(fresh.entropy()),"compact_dl":float(compact.entropy()),
            "live_dl":after})
        if not math.isclose(after,float(fresh.entropy()),rel_tol=1e-10,abs_tol=1e-6):
            emit("LIVE_PROPERTIES",properties(s))
            emit("FRESH_PROPERTIES",properties(fresh))
            emit("COMPACT_PROPERTIES",properties(compact))
            emit("LIVE_ENTROPY_TERMS",scores(s))
            emit("FRESH_ENTROPY_TERMS",scores(fresh))
            emit("COMPACT_ENTROPY_TERMS",scores(compact))
            # Distinguish nominal block capacity from occupied membership.
            try:
                B=int(s.B)
                if 2<=B<=g.num_vertices()+100:
                    matched=make(a,B=B)
                    emit("MATCHED_NOMINAL_B_PROPERTIES",properties(matched))
                    emit("MATCHED_NOMINAL_B_ENTROPY",scores(matched))
                else:emit("NOMINAL_B_RECONSTRUCTION_SKIPPED",B)
            except Exception as e:emit("NOMINAL_B_RECONSTRUCTION_UNAVAILABLE",str(e))
            emit("OBJECTIVE_MISMATCH_REPRODUCED")
            break
        b=canonical(a)
    emit("P1D05_ENTROPY_AUDIT_COMPLETE_NO_OBSERVED_FIT")

if __name__=="__main__":main()
