#!/usr/bin/env python3
"""P1B01: lossless source-scoped harmonization, no equivalence merging or graph fitting."""
import argparse
from collections import Counter,defaultdict
from contextlib import ExitStack
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"diagnostics"))
from p1a04k_scope_feasibility_audit import ROOT,TYPED,STATE,rows,sha,need
from p1a04t_full_sbml_export_comparison import control_role
from p1a04f_typed_source_layer import RO,MOL,EXPECTED
OUT=ROOT/"results/P1B01_SOURCE_SCOPED_HARMONIZATION_V1"
START=time.monotonic()
def js(x):return json.dumps(x,sort_keys=True,ensure_ascii=False,separators=(",",":"))
def emit(*x):print(*x,"ELAPSED_S",round(time.monotonic()-START,2),flush=True)
def digest(b):return hashlib.sha256(b).hexdigest()
def key(*xs):
    need(all(isinstance(x,str) and x and "|" not in x for x in xs),"invalid scoped ID "+repr(xs))
    return "|".join(xs)
def bp(x):return key("BP",x)
def go(m,k,x):return key("GO",m,k,x)
def lock_inputs():
    lockpath=ROOT/"manifests/P1A_SOURCE_LOCK_V1.json"
    lock=json.loads(lockpath.read_text())
    need(lock["status"]=="FROZEN_APPROVED","source lock not approved")
    candidate=ROOT/lock["candidate_path"];manifest=ROOT/lock["manifest_path"]
    need(sha(candidate)==lock["candidate_sha256"] and sha(manifest)==lock["manifest_sha256"],"candidate/manifest hash mismatch")
    c=json.loads(candidate.read_text())
    need(c["git_head"]==lock["candidate_git_commit"] and c["source_manifest_sha256"]==lock["manifest_sha256"] and c["policy_and_registry_sha256"]==lock["policy_and_registry_sha256"],"historical candidate mismatch")
    for name,h in lock["policy_and_registry_sha256"].items():
        b=subprocess.check_output(["git","show",lock["candidate_git_commit"]+":"+name],cwd=ROOT)
        need(digest(b)==h,"historical policy mismatch "+name)
    need(not subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip(),"tracked working tree modified")
    rr=list(rows(manifest));need(len(rr)==len({x["id"] for x in rr})==lock["manifest_rows"]==2085,"source manifest count")
    idx={x["id"]:x for x in rr}
    need(len({x["local_relative_path"] for x in rr})==len(rr),"duplicate manifest paths")
    emit("SOURCE_LOCK_VERIFIED",lock["candidate_git_commit"],lock["manifest_sha256"])
    return lock,idx,sha(lockpath)
def file_for(idx,id,verify=True):
    x=idx[id];p=(ROOT/x["local_relative_path"]).resolve()
    need(p.is_relative_to(ROOT.resolve()) and p.is_file() and p.stat().st_size==int(x["bytes"]),"missing/changed input "+id)
    if verify:need(sha(p)==x["sha256"],"frozen input changed "+id)
    return p
def load_tables(idx):
    for base,v in ((TYPED,"P1A04F_v1"),(STATE,"P1A04H_v1")):
        p=file_for(idx,"DERIVED_MANIFEST:"+v)
        need(p==base/"manifest.json","derived manifest location")
        m=json.loads(p.read_text());need(m["schema_version"]==v,"manifest version")
        for n,x in m["files"].items():
            need(Path(n).name==n,"unsafe table name")
            p=file_for(idx,v+":"+n)
            need(p==base/n and idx[v+":"+n]["sha256"]==x["sha256"] and int(idx[v+":"+n]["bytes"])==x["bytes"],"derived manifest mismatch")
    emit("FROZEN_DERIVED_TABLES_VERIFIED",13)
def self_test():
    need(go("m1","molecule","x")!=go("m2","molecule","x"),"GO model scope")
    need(go("m1","activity","x")!=go("m1","molecule","x"),"GO kind scope")
    need(bp("x")!=go("m1","activity","x"),"source scope")
    need(control_role({"class":"Control","control_type":'["INHIBITION"]'})=="negativeregulator","control role")
    need(RO["RO:0002233"][1]=="reverse" and RO["RO:0002234"][1]=="forward","role direction")
    emit("P1B01_SELF_TEST_PASS")
def build(out,idx):
    nodes={};nc=Counter();ac=Counter();roles=Counter()
    def node(id,source,kind,cl,mid,local,ordinal,payload):
        v=(source,kind,cl,mid,local,str(ordinal),js(payload))
        need(id not in nodes,"duplicate node "+id)
        nodes[id]=v;nc[source,kind]+=1
    def writer(stack,name,cols):
        f=stack.enter_context((out/name).open("w",encoding="utf-8",newline=""))
        w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(cols.split("\t"));return w
    def bp_kind(cl):
        if cl in ("BiochemicalReaction","Degradation","TemplateReaction"):return "operation"
        if cl in ("Protein","Complex","PhysicalEntity","SmallMolecule","Dna","Rna"):return "state"
        if cl in ("Control","Catalysis","TemplateReactionRegulation"):return "control"
        return "auxiliary"
    ops={x["rdf_id"]:x for x in rows(TYPED/"biopax_operations.tsv")}
    states={x["rdf_id"]:x for x in rows(STATE/"biopax_states.tsv")}
    controls={x["rdf_id"]:x for x in rows(TYPED/"biopax_control_records.tsv")}
    need((len(ops),len(states),len(controls))==(16187,61418,10145),"BP catalogue counts")
    for records in (ops,states,controls):
        for id,x in sorted(records.items()):
            node(bp(id),"BP",bp_kind(x["class"]),x["class"],"",id,"",x)
    # GO model-scoped IDs; original activity and molecule objects carry detailed evidence.
    go_compact={(x["model_id"],x["activity_id"]):x for x in rows(TYPED/"gocam_activities.tsv")}
    enablers={(x["model_id"],x["activity_id"]):x for x in rows(STATE/"gocam_enablers.tsv")}
    need(len(go_compact)==len(enablers)==7385,"GO compact activity inventory")
    go_expected={}
    for x in rows(TYPED/"gocam_incidence.tsv"):
        molecular=x["original_predicate"] in MOL
        k=(x["model_id"],int(x["activity_ordinal"]),"molecular_associations" if molecular else "causal_associations",int(x["association_ordinal"]))
        need(k not in go_expected,"duplicate GO association key")
        go_expected[k]=x
    need(len(go_expected)==15491,"GO compact assertion inventory")
    with ExitStack() as stack:
        assertions=writer(stack,"assertions.tsv","assertion_id\tsource\tfrom_id\tto_id\tlayer\tqualitative_effect\toriginal_predicate\tprovenance_json\tevidence_json")
        ctrw=writer(stack,"direct_control_candidates.tsv","event_id\tcontroller_id\trole\tcontrol_id\tcontroller_property_json\tcontrolled_property_json")
        partw=writer(stack,"participant_coefficients.tsv","event_id\tentity_id\toriginal_role\tcoefficient_status\towner_stoichiometry_json")
        compw=writer(stack,"composition.tsv","owner_id\tchild_id\trelation\tstoichiometry_id\tcoefficient")
        modelw=writer(stack,"gocam_model_context.tsv","model_id\tprimary_taxon_json\tindexed_taxa_json\tpayload_sha256\tlocal_relative_path")
        def assertion(id,s,a,b,layer,effect,pred,prov,evidence):
            need(a in nodes and b in nodes,"unresolved assertion endpoints "+id)
            assertions.writerow([id,s,a,b,layer,effect,pred,js(prov),js(evidence)])
            ac[s]+=1;roles[s,layer]+=1
        controllers=defaultdict(list);targets=defaultdict(list);n=0
        for x in rows(TYPED/"biopax_incidence.tsv"):
            for id,cl in ((x["source_id"],x["source_class"]),(x["target_id"],x["target_class"])):
                k=bp(id)
                if k not in nodes:
                    need(bp_kind(cl)=="auxiliary","missing primary BP node "+id)
                    node(k,"BP","auxiliary",cl,"",id,"",{})
                else:need(nodes[k][2]==cl,"BP class conflict")
            prov={"record_ordinal":int(x["record_ordinal"]),"property_ordinal":int(x["property_ordinal"])}
            assertion(key("BP_ASSERT",x["record_ordinal"],x["property_ordinal"]),"BP",bp(x["source_id"]),bp(x["target_id"]),x["layer"],"",x["original_predicate"],prov,None)
            if x["layer"]=="controller_of_record":controllers[x["target_id"]].append((x["source_id"],prov))
            if x["layer"]=="controlled_by_record":targets[x["source_id"]].append((x["target_id"],prov))
            n+=1
        need(n==441156,"BP assertion count")
        pairs=defaultdict(set);control_n=0
        for c,x in sorted(controls.items()):
            need(len(controllers[c])==len(targets[c])==1,"control target cardinality "+c)
            e,ep=targets[c][0];a,ap=controllers[c][0]
            need(e in ops and a in states,"control endpoints")
            z=control_role(x);pairs[z].add((e,a))
            ctrw.writerow([bp(e),bp(a),z,bp(c),js(ap),js(ep)]);control_n+=1
        need(control_n==10145 and {k:len(v) for k,v in pairs.items()}=={"catalyst":6871,"positiveregulator":2162,"negativeregulator":1111},"direct control pair counts")
        n=0
        for x in rows(STATE/"biopax_event_participants.tsv"):
            need(x["event_rdf_id"] in ops and x["entity_rdf_id"] in states,"participant endpoints")
            partw.writerow([bp(x["event_rdf_id"]),bp(x["entity_rdf_id"]),x["role"],x["coefficient_status"],x["owner_stoichiometry_json"]]);n+=1
        need(n==sum(v for (s,z),v in roles.items() if s=="BP" and z in ("material_input","material_output","degradation_input","template","template_product")),"participant count versus typed incidence")
        n=0
        for x in rows(STATE/"biopax_composition.tsv"):
            need(x["owner_id"] in states and x["child_id"] in states,"composition endpoints")
            compw.writerow([bp(x["owner_id"]),bp(x["child_id"]),x["relation"],bp(x["stoichiometry_id"]) if x["stoichiometry_id"] else "",x["coefficient"]]);n+=1
        need(n==76196,"composition count")
        model_count=0;go_preds=Counter();go_seen=set();seen_activities=set()
        inventory=json.loads(file_for(idx,"GO_INVENTORY").read_text())
        need(len(inventory)==1520,"GO inventory size")
        for rec in sorted(inventory,key=lambda v:v["model"]["model_id"]):
            mid=rec["model"]["model_id"];id="GO_MODEL:"+mid;p=file_for(idx,id,False)
            raw=p.read_bytes();need(digest(raw)==idx[id]["sha256"]==rec["sha256"],"GO raw hash")
            obj=json.loads(raw);need(obj["id"]=="gomodel:"+mid,"GO model ID")
            modelw.writerow([mid,js(obj.get("taxon")),js(rec["indexed_taxa"]),rec["sha256"],idx[id]["local_relative_path"]])
            for field,kind in (("activities","activity"),("molecules","molecule")):
                for i,x in enumerate(obj[field]):
                    sid=x["id"];nid=go(mid,kind,sid)
                    node(nid,"GO",kind,kind,mid,sid,i,x)
                    if kind=="activity":
                        k=(mid,sid);compact=go_compact[k];en=enablers[k]
                        need(int(compact["activity_ordinal"])==int(en["activity_ordinal"])==i,"GO activity ordinal")
                        eb=x.get("enabled_by") or {}
                        need(eb==json.loads(en["enabler_json"]),"GO enabler mismatch")
                        need((eb.get("term") or "")==compact["enabler_term"] and (x.get("molecular_function",{}).get("term") or "")==compact["molecular_function"],"GO term mismatch")
                        seen_activities.add(k)
            for i,x in enumerate(obj["activities"]):
                aid=x["id"];nid=go(mid,"activity",aid);eb=x.get("enabled_by") or {}
                # Enabler observation is scoped to this activity. Not an inferred
                # global protein/complex node and not an extra evidence assertion.
                enid=go(mid,"enabler_observation",aid)
                node(enid,"GO","enabler_observation",eb.get("type","unspecified"),mid,aid,i,eb)
                assertion(key("GO_ENABLER",mid,str(i)),"GO_ENABLER",enid,nid,"enabled_by_observation","","enabled_by",{"model_id":mid,"activity_ordinal":i,"json_path":"enabled_by"},eb.get("evidence"))
                for field,targetfield,kind in (("causal_associations","downstream_activity","activity"),("molecular_associations","molecule","molecule")):
                    for j,y in enumerate(x.get(field,[])):
                        pred=y["predicate"];layer,direction,effect=RO[pred];tid=y[targetfield]
                        need((kind=="molecule")== (pred in MOL),"GO role kind")
                        a,b=(go(mid,kind,tid),nid) if direction=="reverse" else (nid,go(mid,kind,tid))
                        k=(mid,i,field,j);v=go_expected[k]
                        need(k not in go_seen,"duplicate GO relation");go_seen.add(k)
                        left,right,lk,rk=(tid,aid,kind,"activity") if direction=="reverse" else (aid,tid,"activity",kind)
                        need((v["source_id"],v["target_id"],v["source_kind"],v["target_kind"],v["layer"],v["qualitative_effect"],v["original_predicate"])==(left,right,lk,rk,layer,effect,pred),"GO typed assertion mismatch")
                        ev=y.get("evidence");need(v["evidence_n"]==(str(len(ev)) if isinstance(ev,list) else ""),"GO evidence count mismatch")
                        assertion(key("GO_ASSERT",mid,str(i),field,str(j)),"GO",a,b,layer,effect,pred,{"model_id":mid,"activity_ordinal":i,"field":field,"association_ordinal":j},ev)
                        go_preds[pred]+=1
            model_count+=1
            if model_count%250==0:emit("GO_MODELS_HARMONIZED",model_count)
        need(model_count==1520 and seen_activities==set(go_compact)==set(enablers) and go_seen==set(go_expected),"GO coverage")
        need(dict(go_preds)==EXPECTED and nc["GO","molecule"]==7322 and nc["GO","activity"]==7385 and nc["GO","enabler_observation"]==7385,"GO inventory")
        need(ac=={"BP":441156,"GO":15491,"GO_ENABLER":7385},"assertion counts")
        nw=writer(stack,"nodes.tsv","node_id\tsource\tkind\tclass\tmodel_id\tlocal_id\tsource_ordinal\tpayload_json")
        for id,v in sorted(nodes.items()):nw.writerow([id,*v])
    summary={"nodes":len(nodes),"node_counts":[[list(k),v] for k,v in sorted(nc.items())],"assertions":dict(ac),"role_counts":[[list(k),v] for k,v in sorted(roles.items())],"direct_control_records":control_n,"direct_control_unique_pairs":{k:len(v) for k,v in sorted(pairs.items())},"go_models":model_count,"scope_filters":0,"merges":0}
    emit("P1B01_COUNTS",summary)
    return summary
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--self-test",action="store_true");ap.add_argument("--execute",action="store_true");a=ap.parse_args()
    self_test()
    if a.self_test:return
    lock,idx,lh=lock_inputs();load_tables(idx)
    if not a.execute:
        emit("P1B01_PLAN_ONLY",OUT.relative_to(ROOT));return
    need(not OUT.exists(),"output exists; no overwrite")
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".p1b01_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir()
        summary=build(tmp,idx)
        m={"schema_version":"P1B01_v1","source_lock_sha256":lh,"source_manifest_sha256":lock["manifest_sha256"],"git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),"summary":summary,"files":{p.name:{"bytes":p.stat().st_size,"sha256":sha(p)} for p in sorted(tmp.iterdir())},"limitations":["BP evidence remains linked by property ordinals and frozen staged records; explicit publication materialization follows","GO enabler observation not canonical molecule","No equivalence merges, scope filters, graph fitting, or inferred coefficients"]}
        (tmp/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared");tmp.rename(OUT)
    emit("P1B01_MANIFEST_SHA256",sha(OUT/"manifest.json"))
    emit("P1B01_SOURCE_SCOPED_HARMONIZATION_COMPLETE; NO_MERGES_OR_FILTERS")
if __name__=="__main__":main()
