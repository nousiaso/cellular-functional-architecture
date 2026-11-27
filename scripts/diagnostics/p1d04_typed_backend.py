#!/usr/bin/env python3
"""Check existing environments and run a synthetic typed-SBM smoke test.
Optional --create is explicit opt-in to a dedicated project Conda prefix.
No biological data, inference outputs or changes to base/other environments.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
PREFIX=ROOT/".conda/cfa_inference_v1"
CACHE=ROOT/".conda/pkgs"

def smoke():
    import graph_tool as gt
    from graph_tool.inference import LayeredBlockState, minimize_blockmodel_dl
    import numpy as np
    import scipy
    gt.seed_rng(20261127)
    rng=random.Random(20261127)
    # Two operation blocks, two molecular-state blocks, three directed roles.
    # Planted labels are used only for evaluating this synthetic test.
    m=60;n=2*m
    g=gt.Graph(directed=True);g.add_vertex(2*n)
    role=g.new_edge_property("int")
    kind=g.new_vertex_property("int")
    kind.a[:n]=0;kind.a[n:]=1
    data=[]
    for e in range(n):
        for i in range(n):
            for k in range(3):
                match=(e//m==i//m) if k!=1 else (e//m!=i//m)
                if rng.random() < (.23 if match else .012):
                    a,b=(n+i,e) if k==1 else (e,n+i)
                    data.append((a,b,k))
    if len(data)!=len(set(data)):raise RuntimeError("duplicate role incidence")
    g.add_edge_list(data,eprops=[role])
    if g.num_edges()!=len(data):raise RuntimeError("edge count")
    if any(kind[e.source()]==kind[e.target()] for e in g.edges()):
        raise RuntimeError("non-bipartite incidence")
    if Counter((int(role[e])) for e in g.edges())!=Counter(k for a,b,k in data):
        raise RuntimeError("role map mismatch")
    print("BACKEND_VERSIONS",gt.__version__,np.__version__,scipy.__version__,flush=True)
    print("SYNTHETIC_GRAPH",g.num_vertices(),g.num_edges(),dict(Counter(k for a,b,k in data)),flush=True)
    args=dict(ec=role,independent=True,deg_corr=True,clabel=kind,pclabel=kind)
    start=LayeredBlockState(g,b=kind,**args)
    null_dl=float(start.entropy())
    t=time.monotonic()
    state=minimize_blockmodel_dl(g,state=LayeredBlockState,state_args=args)
    dl=float(state.entropy());b=np.asarray(state.get_blocks().a).copy()
    if not math.isfinite(dl) or not math.isfinite(null_dl):raise RuntimeError("nonfinite description length")
    if len(b)!=2*n or set(b[:n])&set(b[n:]):raise RuntimeError("type constraint failed")
    # Planted recovery is a diagnostic, not a guarantee of global optimum.
    def purity(labels):
        groups=defaultdict(Counter)
        for i,x in enumerate(labels):groups[int(x)][i//m]+=1
        return sum(max(z.values()) for z in groups.values())/len(labels)
    print("TYPED_SBM_SMOKE",json.dumps(dict(
        degree_corrected=True,independent_layers=True,operation_blocks=len(set(b[n:])),
        state_blocks=len(set(b[:n])),type_only_dl=null_dl,fitted_dl=dl,
        synthetic_operation_purity=purity(b[n:]),synthetic_state_purity=purity(b[:n]),
        seconds=round(time.monotonic()-t,3),type_constraint_pass=True)),flush=True)
    print("P1D04_SYNTHETIC_BACKEND_PASS",flush=True)

def run(prefix,conda):
    print("TEST_ENVIRONMENT",str(prefix),flush=True)
    x=subprocess.run([conda,"run","--no-capture-output","-p",str(prefix),"python","-u",str(Path(__file__).resolve()),"--smoke"],timeout=300)
    return x.returncode==0

def main():
    p=argparse.ArgumentParser();p.add_argument("--smoke",action="store_true");p.add_argument("--create",action="store_true")
    a=p.parse_args()
    if a.smoke:smoke();return
    conda=shutil.which("conda")
    if not conda:raise RuntimeError("conda not found")
    info=json.loads(subprocess.check_output([conda,"info","--json"],text=True,timeout=60))
    prefixes=sorted({Path(e) for e in info.get("envs",[]) if Path(e).is_dir()})
    if PREFIX.is_dir() and PREFIX not in prefixes:prefixes.append(PREFIX)
    candidates=[]
    for e in prefixes:
        if (e/"bin/python").is_file() and any(e.glob("lib/python*/site-packages/graph_tool/__init__.py")):
            candidates.append(e)
    print("EXISTING_GRAPH_TOOL_CANDIDATES",len(candidates),[x.name for x in candidates],flush=True)
    for e in candidates:
        if run(e,conda):
            print("P1D04_EXISTING_BACKEND_PASS",str(e),flush=True);return
    if not a.create:
        print("P1D04_BACKEND_UNAVAILABLE_NO_INSTALLATION_PERFORMED",flush=True)
        return
    if PREFIX.exists():
        raise RuntimeError("project prefix already exists; inspect before changing; no overwrite")
    CACHE.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env["CONDA_PKGS_DIRS"]=str(CACHE)
    # Pin candidate graph-tool major/minor and Python minor; exact build lock
    # recorded after successful solve. No channels inherited from base.
    subprocess.run([conda,"create","--yes","--prefix",str(PREFIX),
                    "--override-channels","--channel","conda-forge","--strict-channel-priority",
                    "python=3.12","graph-tool-base=3.9","numpy","scipy"],env=env,check=True)
    if not run(PREFIX,conda):raise RuntimeError("new backend smoke failed")
    lock=PREFIX.parent/"cfa_inference_v1.explicit.txt"
    text=subprocess.check_output([conda,"list","--prefix",str(PREFIX),"--explicit"],text=True)
    with lock.open("x") as f:f.write(text)
    print("P1D04_CREATED_BACKEND_PASS",str(PREFIX),"EXPLICIT_LOCK",str(lock),flush=True)
if __name__=="__main__":main()
