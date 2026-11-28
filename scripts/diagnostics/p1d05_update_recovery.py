#!/usr/bin/env python3
"""Bounded recovery gate: equivalent full-objective construction and update methods.
Randomized N2 workload only; no biological partition, no saved memberships.
"""
import math
import time
import numpy as np
import graph_tool as gt
from graph_tool.inference import LayeredBlockState
from p1d05_native_scale_pilot import make_graph, RECEIPT_HASH
from p1d05_entropy_audit import canonical
from p1d02_mixing_evidence_pilot import ROOT, load, chain, keyfun, sha, need, emit

def close(a,b):
    return math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-6)

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
    orig=kind.a.copy();oe=ec.a.copy();ge=g.get_edges().copy()
    def vp(a):
        p=g.new_vertex_property("int64_t");p.a[:]=a;return p
    def make(a,mode):
        if mode=="original":cl=pc=vp(orig)
        else:cl,pc=vp(orig),vp(orig)
        args=dict(ec=ec,independent=True,deg_corr=True,clabel=cl,pclabel=pc)
        if mode=="explicit_unit_weights":
            ew=g.new_edge_property("int64_t");ew.a[:]=1
            vw=g.new_vertex_property("int64_t");vw.a[:]=1
            args.update(eweight=ew,vweight=vw)
        return LayeredBlockState(g,b=vp(canonical(a)),**args)
    rng=np.random.default_rng(20261130)
    split=np.concatenate((rng.integers(0,4,ne),rng.integers(4,8,len(R)))).astype("int64")
    starts={"type_only":orig,"random_8_blocks":split}
    ref={name:float(make(a,"original").entropy()) for name,a in starts.items()}
    emit("REFERENCE_FULL_DL",ref)
    results=[]
    for mode in ("original","separate_label_maps","explicit_unit_weights"):
        scores={name:float(make(a,mode).entropy()) for name,a in starts.items()}
        equivalent=all(close(scores[k],ref[k]) for k in ref)
        emit("CONSTRUCTOR_EQUIVALENCE",mode,scores,equivalent)
        if not equivalent:
            emit("CONSTRUCTOR_REJECTED_DIFFERENT_OBJECTIVE",mode);continue
        for method in ("mcmc_sweep","multiflip_mcmc_sweep"):
            passed=True;total_moves=0;nonzero_delta=0
            for name,a in starts.items():
                if not passed:break
                for seed in (20261130,20261201):
                    if not passed:break
                    gt.seed_rng(seed);s=make(a,mode)
                    for i in range(3):
                        before=float(s.entropy());b0=canonical(s.get_blocks().a.copy())
                        # Finite temperature exercises accepted updates; this is
                        # correctness testing, not a fitted partition or sampler.
                        ret=getattr(s,method)(beta=1.0,niter=1)
                        after=float(s.entropy());b=s.get_blocks().a.copy()
                        fresh=float(make(b,mode).entropy())
                        delta=after-before
                        invariant=(np.array_equal(kind.a,orig) and np.array_equal(ec.a,oe)
                                   and np.array_equal(g.get_edges(),ge))
                        # Verify against the original constructor as well as
                        # the candidate, not just against itself.
                        original_score=float(make(b,"original").entropy())
                        typeok=not bool(set(b[:ne])&set(b[ne:]))
                        ok=(invariant and typeok and math.isfinite(after)
                            and close(after,fresh) and close(fresh,original_score)
                            and close(delta,ret[0])
                            and (not np.array_equal(b0,canonical(b)) or close(delta,0.0)))
                        total_moves+=int(ret[2]);nonzero_delta+=int(not close(delta,0.0))
                        emit("UPDATE_CHECK",mode,method,name,seed,i+1,{
                            "before":before,"live":after,"fresh":fresh,"original_constructor":original_score,
                            "reported_delta":float(ret[0]),"actual_delta":delta,
                            "moves":int(ret[2]),"type_ok":typeok,"invariant":invariant,"pass":ok})
                        if not ok:passed=False;break
            result={"mode":mode,"method":method,"all_checks_pass":passed,
                    "reported_moves":total_moves,"nonzero_delta_checks":nonzero_delta,
                    "candidate_for_further_validation":passed and total_moves>0 and nonzero_delta>0}
            results.append(result);emit("RECOVERY_CANDIDATE",result)
    emit("P1D05_RECOVERY_SUMMARY",results)
    emit("P1D05_RECOVERY_GATE_COMPLETE_NO_OBSERVED_FIT")
if __name__=="__main__":main()
