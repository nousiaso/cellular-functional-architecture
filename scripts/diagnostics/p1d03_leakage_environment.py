#!/usr/bin/env python3
"""Read-only leakage and inference-environment feasibility; no split or fit."""
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from p1d02_mixing_evidence_pilot import load, B01, B03, C01, rows, arr, need, emit, DSU

def main():
    R,edge,sp,rc,escope,root,direct=load()
    E={r["state_id"]:int(r["state_index"]) for r in rows(C01/"state_index.tsv")}
    n=len(R);pub=defaultdict(set)
    for r in rows(B03/"bp_publication_links.tsv"):
        for a in arr(r,"namespaces_json"):
            for b in arr(r,"identifiers_json"):
                if isinstance(a,str) and a and isinstance(b,str) and b:pub[a.lower(),b].add(R[r["owner_id"]])
    known=set().union(*pub.values());d=DSU(n)
    for es in pub.values():
        a=min(es)
        for b in es:d.union(a,b)
    need(d.report(known)["largest_10"][0]==4962,"publication baseline")
    original={i:d.find(i) for i in range(n)}
    byop=defaultdict(set)
    for k,e,i in edge:byop[i].add((k,e))
    need(len(byop)==n,"role coverage")
    # Same operation class and exact complete observed binary role profile.
    # Candidate duplicates only; identity and source records are never merged.
    fp={}
    for r in rows(B01/"nodes.tsv"):
        if r["node_id"] in E:
            x=json.loads(r["payload_json"]);e=E[r["node_id"]]
            f=x.get("fingerprint_candidate")
            if isinstance(f,str) and f:fp[e]=f
    emit("STATE_FINGERPRINT_COVERAGE",len(fp),len(E))
    for mode in ("exact_state_roles","candidate_fingerprint_roles"):
        if mode=="candidate_fingerprint_roles" and not fp:
            emit("FINGERPRINT_PROFILE_UNAVAILABLE");continue
        groups=defaultdict(list)
        for i,z in sorted(byop.items()):
            if mode=="exact_state_roles":sig=tuple(sorted(z))
            else:sig=tuple(sorted({(k,"FP:"+fp[e] if e in fp else "ID:"+str(e)) for k,e in z}))
            groups[rc[i],sig].append(i)
        dup=[g for g in groups.values() if len(g)>1]
        bridge=[g for g in dup if len({original[i] for i in g})>1]
        for g in dup:
            for i in g[1:]:d.union(g[0],i)
        emit("MOLECULAR_DUPLICATE_GROUP_PILOT",mode,{"distinct_profiles":len(groups),
             "multi_operation_groups":len(dup),"operations_in_groups":sum(map(len,dup)),
             "largest_groups":sorted(map(len,dup),reverse=True)[:10],
             "groups_crossing_original_publication_components":len(bridge),
             "operations_in_crossing_groups":sum(map(len,bridge)),
             "combined_components_all":d.report(range(n)),"combined_components_known":d.report(known)})
    # The publication-connected component includes 4962 operations; report
    # exact profile linkage to unknown evidence rather than treat them as independent.
    comps=defaultdict(set)
    for i in range(n):comps[d.find(i)].add(i)
    emit("COMBINED_EVIDENCE_UNKNOWN",{"unknown":n-len(known),
        "unknown_in_components_with_known":sum(len(z-known) for z in comps.values() if z&known)})
    # Role-masked new-operation assignment feasibility, using disjoint role
    # views only as a proposed design. No held-out labels are used for fitting.
    for mask in ({"material_output"},{"material_input"},{"catalyst","positiveregulator","negativeregulator"}):
        seen=defaultdict(set)
        for i,z in byop.items():
            for k,e in z:
                if k not in mask:seen[k,e].add(d.find(i))
        eligible=0;has_target=0;has_assignment=0;othergroup=0
        for i,z in byop.items():
            a={x for x in z if x[0] not in mask};b={x for x in z if x[0] in mask}
            has_target+=bool(b);has_assignment+=bool(a)
            eligible+=bool(a and b)
            othergroup+=bool(b and any(seen[x]-{d.find(i)} for x in a))
        emit("DISJOINT_ROLE_VIEW_CAPACITY",sorted(mask),{"has_target":has_target,"has_assignment":has_assignment,
             "both":eligible,"target_and_assignment_feature_seen_outside_component":othergroup})
    # No environment installation or network calls. Existing package caches
    # and available tools determine a safe next installation path.
    emit("ENVIRONMENT",{"python":sys.version.split()[0],"conda":bool(shutil.which("conda")),
        "mamba":bool(shutil.which("mamba")),"micromamba":bool(shutil.which("micromamba"))})
    for mod in ("numpy","scipy","networkx","sklearn","graph_tool","igraph","leidenalg","infomap"):
        try:
            from importlib.metadata import version
            v=version({"graph_tool":"graph-tool","sklearn":"scikit-learn","igraph":"igraph"}.get(mod,mod))
        except Exception:v=None
        emit("PACKAGE_VERSION",mod,v)
    # Report only environment names and channel names, not private paths or URLs.
    if shutil.which("conda"):
        x=subprocess.run(["conda","info","--json"],capture_output=True,text=True,timeout=60)
        if x.returncode==0:
            z=json.loads(x.stdout)
            emit("CONDA_CAPACITY",{"version":z.get("conda_version"),"platform":z.get("platform"),
              "environment_names":[Path(e).name for e in z.get("envs",[])],
              "package_cache_count":len(z.get("pkgs_dirs",[]))})
        else:emit("CONDA_INFO_FAILED",x.returncode)
    emit("P1D03_READONLY_COMPLETE_NO_SPLIT_OR_FIT")
if __name__=="__main__":main()
