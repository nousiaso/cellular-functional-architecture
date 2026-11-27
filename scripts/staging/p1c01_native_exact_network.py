#!/usr/bin/env python3
"""P1C01: native binary role incidence and separate exact-state projections.
All source operations retained. No graph fitting or eligibility selection.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile

from p1b04_harmonization_certificate import ROOT, B01, B02, B03, OUT as B04, verify_chain, rows, sha, need, emit

OUT=ROOT/"results/P1C01_NATIVE_EXACT_NETWORK_V1"
LOCK=ROOT/"manifests/P1B_HARMONIZATION_LOCK_V1.json"
B04_HASH="82fc4ac1213edcfe4ea68c6b447870ea5e2027da9379f9276b52ef1be1bf6c7a"
csv.field_size_limit(64*1024*1024)
ROLES=("material_input","material_output","catalyst","positiveregulator",
       "negativeregulator","degradation_input","template","template_product")
LAYERS=(("flow","material_output","material_input"),
        ("catalytic_dependency","material_output","catalyst"),
        ("positive_regulation_candidate","material_output","positiveregulator"),
        ("negative_regulation_candidate","material_output","negativeregulator"),
        ("shared_catalyst_baseline","catalyst","catalyst"))

def verify_frozen():
    _,lh=verify_chain()
    lock=json.loads(LOCK.read_text())
    need(lock["schema_version"]=="P1B_LOCK_v1" and lock["status"]=="FROZEN_APPROVED","P1B freeze")
    need(len(lock["manifests"])==4,"P1B lock inventory")
    for r,d in zip(lock["manifests"],(B01,B02,B03,B04)):
        p=d/"manifest.json"
        need(r["path"]==str(p.relative_to(ROOT)) and sha(p)==r["sha256"],"P1B manifest "+str(d))
    p=B04/"manifest.json";need(sha(p)==B04_HASH,"P1B04 fingerprint")
    m=json.loads(p.read_text())
    need(m["schema_version"]=="P1B04_v1" and m["status"]=="TECHNICAL_PASS_REVIEW_REQUIRED"
         and m["source_lock_sha256"]==lh and m["git_head"]==lock["certificate_execution_commit"],"P1B04 lineage")
    for n,x in m["files"].items():
        need(Path(n).name==n,"unsafe certificate file")
        f=B04/n;need(f.stat().st_size==x["bytes"] and sha(f)==x["sha256"],"P1B04 file "+n)
    emit("P1B_FROZEN_CHAIN_VERIFIED",len(m["files"]))
    return lh,sha(LOCK)

def projected_rows(left,right,n):
    """Stream one source row at a time. All witnesses reconstructible from incidence."""
    by_op=defaultdict(list);stats={}
    for e in sorted(left.keys()&right.keys()):
        l,r=left[e],right[e]
        d=len(l|r);b=len(l)*len(r)-len(l&r)
        need(0<=d<=n and b>=0,"degree arithmetic")
        if not b:continue
        stats[e]=(d,b,math.log((n+1)/(d+1)),1/b)
        for i in sorted(l):by_op[i].append(e)
    def generate():
        for i in sorted(by_op):
            rr={}
            for e in by_op[i]:
                _,_,w,z=stats[e]
                for j in sorted(right[e]):
                    if i==j:continue
                    a=rr.setdefault(j,[0,0.0,0.0])
                    a[0]+=1;a[1]+=w;a[2]+=z
            for j,a in sorted(rr.items()):yield i,j,*a
    return stats,generate()

def self_test():
    stats,rr=projected_rows({0:{0,1},1:{0}},{0:{1,2},1:{2}},3)
    x={(i,j):(c,w,z) for i,j,c,w,z in rr}
    need(set(x)=={(0,1),(0,2),(1,2)} and x[0,2][0]==2,"projection fixture")
    need(stats[0][1]==3 and stats[1][1]==1,"bridge denominator")
    need(math.isclose(sum(v[2] for v in x.values()),2),"bridge normalized mass")
    need(math.isclose(x[0,2][1],math.log(4/3)),"frequency mass")
    s,r=projected_rows({0:{1}},{0:{1}},3)
    need(not s and not list(r),"self-only state")
    emit("P1C01_SELF_TEST_PASS")

def build(out):
    states={};ops={}
    for r in rows(B01/"nodes.tsv"):
        if r["source"]!="BP":continue
        if r["kind"]=="state":states[r["node_id"]]=r["class"]
        elif r["kind"]=="operation":ops[r["node_id"]]=r["class"]
    need((len(states),len(ops))==(61418,16187),"node universe")
    eids=sorted(states);rids=sorted(ops);n=len(rids)
    ei={e:i for i,e in enumerate(eids)};ri={e:i for i,e in enumerate(rids)}
    inc={k:Counter() for k in ROLES}
    for r in rows(B01/"assertions.tsv"):
        if r["source"]!="BP":continue
        k=r["layer"];a,b=r["from_id"],r["to_id"]
        if k not in ("material_input","material_output","degradation_input","template","template_product"):continue
        e,op=(a,b) if k in ("material_input","degradation_input","template") else (b,a)
        need(e in ei and op in ri,"material endpoint")
        inc[k][ei[e],ri[op]]+=1
    ctrl=0
    for r in rows(B01/"direct_control_candidates.tsv"):
        k=r["role"];e=r["controller_id"];op=r["event_id"]
        need(k in ("catalyst","positiveregulator","negativeregulator") and e in ei and op in ri,"control endpoint")
        inc[k][ei[e],ri[op]]+=1;ctrl+=1
    need(ctrl==10145,"control occurrences")
    raw={k:sum(x.values()) for k,x in inc.items()};unique={k:len(x) for k,x in inc.items()}
    need(all(raw[k]==v for k,v in {"material_input":31282,"material_output":24986,"degradation_input":20,"template":0,"template_product":8}.items()),"material counts")
    need(all(unique[k]==v for k,v in {"catalyst":6871,"positiveregulator":2162,"negativeregulator":1111}.items()),"control counts")
    def writer(f,header):
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(header.split("\t"));return w
    for fn,ids,classes,header in (("state_index.tsv",eids,states,"state_index\tstate_id\tclass"),("operation_index.tsv",rids,ops,"operation_index\toperation_id\tclass")):
        with (out/fn).open("w",newline="") as f:
            w=writer(f,header)
            for i,e in enumerate(ids):w.writerow([i,e,classes[e]])
    by_state={};active=set()
    with (out/"typed_incidence.tsv").open("w",newline="") as f:
        w=writer(f,"role\tstate_index\toperation_index\tbinary_value\tsource_record_occurrences")
        for k in ROLES:
            z=defaultdict(set)
            for (e,op),m in sorted(inc[k].items()):
                w.writerow([k,e,op,1,m]);z[e].add(op);active.add(op)
            by_state[k]=z
    emit("NATIVE_INCIDENCE_PASS",unique,"raw",raw,"incident_operations",len(active))
    summaries={}
    with (out/"exact_projection_candidates.tsv").open("w",newline="") as f,(out/"projection_state_weights.tsv").open("w",newline="") as g:
        w=writer(f,"layer\tfrom_operation_index\tto_operation_index\tunweighted_count\tinverse_frequency_weight\tbridge_normalized_weight")
        q=writer(g,"layer\tstate_index\tleft_degree\tright_degree\tunion_degree\tordered_nonself_bridge_occurrences\tinverse_frequency_weight\tbridge_normalized_weight")
        for layer,l,r in LAYERS:
            left,right=by_state[l],by_state[r]
            stats,rr=projected_rows(left,right,n)
            for e,(d,b,a,z) in sorted(stats.items()):
                q.writerow([layer,e,len(left[e]),len(right[e]),d,b,format(a,".17g"),format(z,".17g")])
            edges=count=0;mass=idf_mass=0.0;src=set();dst=set()
            for i,j,c,a,z in rr:
                w.writerow([layer,i,j,c,format(a,".17g"),format(z,".17g")])
                edges+=1;count+=c;mass+=z;idf_mass+=a;src.add(i);dst.add(j)
            need(count==sum(x[1] for x in stats.values()),"witness mass "+layer)
            need(math.isclose(mass,len(stats),rel_tol=1e-9,abs_tol=1e-8),"normalized mass "+layer)
            need(math.isclose(idf_mass,sum(x[1]*x[2] for x in stats.values()),rel_tol=1e-9,abs_tol=1e-8),"frequency mass "+layer)
            if layer=="flow":need(count==345409,"frozen exact material bridge count")
            s={"unique_directed_nonself_pairs":edges,"ordered_nonself_witness_occurrences":count,
               "bridging_states":len(stats),"source_operations":len(src),"target_operations":len(dst),
               "incident_operations":len(src|dst),"isolated_operations":n-len(src|dst),
               "bridge_normalized_total_mass":mass,
               "top_5_state_occurrences":[[eids[e],x[1]] for e,x in sorted(stats.items(),key=lambda t:(-t[1][1],t[0]))[:5]]}
            summaries[layer]=s;emit("EXACT_LAYER",layer,s)
    return {"states":len(eids),"operations":n,"native_unique_role_pairs":unique,"native_raw_role_occurrences":raw,
            "operations_with_any_native_role":len(active),"exact_layers":summaries,"scope_filters":0,"primary_model_selected":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--execute",action="store_true");p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    lh,bh=verify_frozen()
    if not a.execute:emit("PLAN_ONLY",OUT.relative_to(ROOT));return
    need(not OUT.exists(),"output exists; no overwrite")
    with tempfile.TemporaryDirectory(prefix=".p1c01_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir()
        s=build(tmp)
        m={"schema_version":"P1C01_v1","status":"CANDIDATE_NETWORKS_REVIEW_REQUIRED",
           "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
           "source_lock_sha256":lh,"harmonization_lock_sha256":bh,"P1B04_manifest_sha256":B04_HASH,
           "index_base":0,"projection_definition":{"roles":LAYERS,"inverse_frequency":"ln((N+1)/(d+1)), d=union of left/right role operation sets, N=all 16187 operations","bridge_normalized":"1/b, b=ordered nonself left/right occurrences","self_pairs":"excluded only from projections","witness_reconstruction":"intersection of left-role states at source and right-role states at target in typed_incidence.tsv"},
           "summary":s,"files":{f.name:{"bytes":f.stat().st_size,"sha256":sha(f)} for f in sorted(tmp.iterdir())},
           "limits":["Full source universe diagnostic, not final eligibility or confirmatory training graph","No pathway, precedence, GO-CAM, gene or publication labels used for edges","Native composition and original assertion provenance remain in frozen 1B; no OR/AND expansion yet","Degradation/template roles separate from ordinary material flow","Signs represented by separate nonnegative layers; shared catalyst baseline symmetric","Projected links are neither flux nor independent observations","No weights selected or graph models fitted"]}
        (tmp/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared");tmp.rename(OUT)
    emit("P1C01_MANIFEST_SHA256",sha(OUT/"manifest.json"));emit("P1C01_COMPLETE_REVIEW_REQUIRED")
if __name__=="__main__":main()
