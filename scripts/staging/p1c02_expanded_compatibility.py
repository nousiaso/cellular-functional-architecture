#!/usr/bin/env python3
"""P1C02: deduplicated OR and direct AND candidate layers, sparse context.
Retains molecular witnesses; never equates states or adds implied chemistry.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile
from p1c01_native_exact_network import ROOT, B01, B02, OUT as C01, verify_frozen, rows, sha, need, emit
from p1b02_molecular_compatibility import STRICT

OUT=ROOT/"results/P1C02_EXPANDED_COMPATIBILITY_V1"
C01_HASH="c19922ead412401c862b11faad63a4b9c1de7cc683e2d91c8695ac753938b09c"

def js(x):return json.dumps(x,separators=(",",":"),sort_keys=True)
def writer(f,header):
    w=csv.writer(f,delimiter="\t",lineterminator="\n");w.writerow(header.split("\t"));return w

def verify():
    lh,bh=verify_frozen()
    p=C01/"manifest.json";need(sha(p)==C01_HASH,"C01 manifest")
    m=json.loads(p.read_text())
    need(m["schema_version"]=="P1C01_v1" and m["source_lock_sha256"]==lh and m["harmonization_lock_sha256"]==bh,"C01 lineage")
    for n,x in m["files"].items():
        need(Path(n).name==n,"unsafe filename")
        f=C01/n;need(f.stat().st_size==x["bytes"] and sha(f)==x["sha256"],"C01 file "+n)
    need(len(m["files"])==5,"C01 file count")
    emit("P1C01_MANIFEST_VERIFIED",len(m["files"]))
    return lh,bh

def or_pairs(real,supported,O,I):
    """Index by terminal. Each state pair appears once, with all shared leaves."""
    pairs=defaultdict(set);good=defaultdict(set);n=Counter()
    for leaf,es in real.items():
        for a in sorted(e for e in es if O.get(e)):
            for b in sorted(e for e in es if I.get(e)):
                if a==b:continue
                pairs[a,b].add(leaf);n["all_pairs"]+=1
                mass=len(O[a])*len(I[b])-len(O[a]&I[b]);n["all_occurrences"]+=mass
                if a in supported.get(leaf,set()) and b in supported.get(leaf,set()):
                    good[a,b].add(leaf);n["supported_pairs"]+=1;n["supported_occurrences"]+=mass
    return pairs,good,n

def self_test():
    O={0:{0},1:{0,1}};I={1:{1,2},2:{2}}
    pairs,good,n=or_pairs({4:{0,1,2},5:{0,2}},{4:{0,1},5:{0,2}},O,I)
    need(pairs[0,2]=={4,5} and good[0,2]=={5},"OR multiple leaves")
    need(n["all_pairs"]==4 and n["supported_pairs"]==2,"OR pair counts")
    need((1,1) not in pairs,"OR exact state excluded")
    emit("P1C02_SELF_TEST_PASS")

def build(out):
    idx={r["state_id"]:int(r["state_index"]) for r in rows(C01/"state_index.tsv")}
    ri={r["operation_id"]:int(r["operation_index"]) for r in rows(C01/"operation_index.tsv")}
    need(set(idx.values())==set(range(61418)) and set(ri.values())==set(range(16187)),"indices")
    roles=defaultdict(lambda:defaultdict(set))
    for r in rows(C01/"typed_incidence.tsv"):
        need(r["binary_value"]=="1","binary incidence")
        roles[r["role"]][int(r["state_index"])].add(int(r["operation_index"]))
    O=roles["material_output"];I=roles["material_input"];N=len(ri)
    def b(a,c):return len(O.get(a,set()))*len(I.get(c,set()))-len(O.get(a,set())&I.get(c,set()))
    common=O.keys()&I.keys()
    need(len(common)==8098 and sum(b(e,e)==0 for e in common)==12,"exact shared state definition")
    exact={(int(r["from_operation_index"]),int(r["to_operation_index"])) for r in rows(C01/"exact_projection_candidates.tsv") if r["layer"]=="flow"}
    need(len(exact)==339542,"exact flow pair count")
    real=defaultdict(set);good=defaultdict(set);term=defaultdict(set);goodterm=defaultdict(set)
    for r in rows(B02/"or_realizations.tsv"):
        a=idx[r["set_id"]];c=idx[r["terminal_id"]]
        term[a].add(c);real[c].add(a)
        if r["location_path_class"]=="location_supported":
            good[c].add(a);goodterm[a].add(c)
        else:need(r["location_path_class"]=="location_uncertain","OR context class")
    owners=set(term)
    for e in idx.values():
        if e not in owners:
            term[e].add(e);goodterm[e].add(e);real[e].add(e);good[e].add(e)
    pairs,supported,raw=or_pairs(real,good,O,I)
    need(raw=={"all_pairs":7011,"all_occurrences":22890,"supported_pairs":6556,"supported_occurrences":22352},"OR preflight reconciliation "+repr(raw))
    need(max(map(len,real.values()))==28,"OR inverse membership")
    # Keep all direct AND requirements, including nonbridging siblings.
    components=defaultdict(set);androws={};andraw=Counter()
    for r in rows(B02/"composition_interfaces.tsv"):
        if r["relation"]!="component":continue
        a=idx[r["child_id"]];c=idx[r["owner_id"]]
        components[c].add(a)
        if b(a,c):
            need((a,c) not in androws,"duplicate direct AND")
            androws[a,c]=r;andraw[r["location_status"]]+=b(a,c)
    need(andraw=={"exact_source_location":41003,"different_recorded_location":27304},"AND preflight")
    assembly={(idx[r["complex_id"]],ri[r["operation_id"]]):r["status"] for r in rows(B02/"complex_production_candidates.tsv")}
    # Each (a,c) counted once regardless of terminal multiplicity. Distinct
    # possible vs recorded inclusion and location support are independent flags.
    # We never infer that the recorded alternatives exhaust biology.
    layers=defaultdict(lambda:defaultdict(list))
    paircounts=Counter();occ=Counter()
    with (out/"or_state_pair_witnesses.tsv").open("w",newline="") as f:
        w=writer(f,"output_state_index\tinput_state_index\tshared_terminal_indices_json\tlocation_supported_terminal_indices_json\toutput_recorded_alternatives_subset_input\toutput_all_recorded_alternatives_location_supported\tordered_nonself_operation_occurrences")
        for (a,c),leaves in sorted(pairs.items()):
            gs=supported.get((a,c),set());mass=b(a,c)
            inclusion=term[a]<=term[c]
            allgood=term[a]<=gs
            w.writerow([a,c,js(sorted(leaves)),js(sorted(gs)),int(inclusion),int(allgood),mass])
            if not mass:continue
            tag="OR_location_supported" if gs else "OR_uncertain_only"
            paircounts[tag]+=1;occ[tag]+=mass
            for i in O[a]:
                for j in I[c]:
                    if i!=j:layers[tag][i,j].append((a,c))
    with (out/"and_state_pair_witnesses.tsv").open("w",newline="") as f:
        w=writer(f,"output_component_state_index\tinput_complex_state_index\tlocation_status\tfeature_id_status\tstoichiometry_id\tcoefficient\tother_direct_required_component_indices_json\tordered_nonself_operation_occurrences")
        for (a,c),r in sorted(androws.items()):
            mass=b(a,c);tag="AND_location_supported" if r["location_status"] in STRICT else "AND_uncertain_only"
            w.writerow([a,c,r["location_status"],r["feature_id_status"],r["stoichiometry_id"],r["coefficient"],js(sorted(components[c]-{a})),mass])
            paircounts[tag]+=1;occ[tag]+=mass
            for i in O[a]:
                for j in I[c]:
                    if i!=j:layers[tag][i,j].append((a,c))
    # The uncertain-only suffix concerns state-pair witnesses, not operation
    # pairs: an operation pair can have supported AND/OR witnesses too.
    # Keep exact-overlap and inter-layer overlap rather than summing layers.
    summary={}
    with (out/"expanded_operation_pairs.tsv").open("w",newline="") as f:
        w=writer(f,"layer\tfrom_operation_index\tto_operation_index\tstate_pair_witnesses_json\tn_unique_state_pairs\texact_flow_pair\tother_layers_json\tbridge_normalized_state_pair_weight")
        for tag,edges in sorted(layers.items()):
            total=0;mass=0.0;nodes=set();added=0
            for (i,j),ws in sorted(edges.items()):
                need(len(ws)==len(set(ws)),"duplicate state pair witness")
                weight=sum(1/b(a,c) for a,c in ws);mass+=weight
                ex=(i,j) in exact;added+=not ex;nodes.update((i,j));total+=len(ws)
                others=sorted(k for k,z in layers.items() if k!=tag and (i,j) in z)
                w.writerow([tag,i,j,js(sorted(ws)),len(ws),int(ex),js(others),format(weight,".17g")])
            need(total==occ[tag] and math.isclose(mass,paircounts[tag],rel_tol=1e-9),"expanded witness mass")
            summary[tag]={"unique_operation_pairs":len(edges),"new_relative_to_exact_flow":added,"overlaps_exact_flow":len(edges)-added,
                          "incident_operations":len(nodes),"unique_state_pairs_with_nonself_witnesses":paircounts[tag],
                          "operation_state_pair_witness_occurrences":total,"bridge_normalized_state_pair_mass":mass}
            emit("EXPANDED_LAYER",tag,summary[tag])
    all_or=set().union(*(set(z) for k,z in layers.items() if k.startswith("OR_")))
    all_and=set().union(*(set(z) for k,z in layers.items() if k.startswith("AND_")))
    # Explicit complex production candidates remain separate evidence.
    # Count whether the upstream operation already records output of complex.
    and_explicit=sum(any((c,i) in assembly for a,c in ws) for k,z in layers.items() if k.startswith("AND_") for (i,j),ws in z.items())
    # U is binary location x operation with typed role witnesses; no U^T U expansion.
    context=defaultdict(set);witness=defaultdict(set);name={}
    for r in rows(B01/"nodes.tsv"):
        if r["node_id"] not in idx:continue
        e=idx[r["node_id"]];p=json.loads(r["payload_json"])
        for loc in json.loads(p["location_ids_json"]):
            context[loc] # retain empty groups
            for role,z in roles.items():
                for op in z.get(e,set()):
                    context[loc].add(op);witness[loc,op].add((e,role))
        if e in common:name[e]=json.loads(p["name_json"])
    upper=sum(len(v)*(len(v)-1) for v in context.values())
    need(len(context)==138 and upper==115947160,"context preflight")
    with (out/"context_incidence.tsv").open("w",newline="") as f:
        w=writer(f,"location_id\toperation_index\tstate_role_witnesses_json")
        for (loc,op),ws in sorted(witness.items()):
            w.writerow([loc,op,js(sorted(ws))])
    with (out/"context_groups.tsv").open("w",newline="") as f:
        w=writer(f,"location_id\toperation_count\tordered_nonself_pair_upper_bound")
        for loc,es in sorted(context.items()):w.writerow([loc,len(es),len(es)*(len(es)-1)])
    inv={v:k for k,v in idx.items()}
    s={"OR_preflight":dict(raw),"AND_preflight":dict(andraw),"expanded_layers":summary,
       "OR_union_operation_pairs":len(all_or),"OR_union_added_only":len(all_or-exact),
       "AND_union_operation_pairs":len(all_and),"AND_union_added_only":len(all_and-exact),
       "OR_AND_overlap":len(all_or&all_and),"expanded_union_added_only":len((all_or|all_and)-exact),
       "AND_layer_pair_occurrences_with_upstream_explicit_complex_output":and_explicit,
       "context_groups":len(context),"context_incidence_nnz":len(witness),"context_pair_occurrence_upper_bound":upper,
       "context_pairs_materialized":0,"exact_shared_states":len(common),"exact_self_only_states":12,
       "top_exact_flow_states":[[inv[e],name[e],b(e,e)] for e in sorted(common,key=lambda e:(-b(e,e),e))[:10]],
       "scope_filters":0,"primary_model_selected":False}
    emit("P1C02_COUNTS",s);return s

def main():
    p=argparse.ArgumentParser();p.add_argument("--execute",action="store_true");p.add_argument("--self-test",action="store_true");a=p.parse_args()
    self_test()
    if a.self_test:return
    lh,bh=verify()
    if not a.execute:emit("PLAN_ONLY",OUT.relative_to(ROOT));return
    need(not OUT.exists(),"output exists; no overwrite")
    with tempfile.TemporaryDirectory(prefix=".p1c02_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir();s=build(tmp)
        m={"schema_version":"P1C02_v1","status":"CANDIDATE_LAYERS_REVIEW_REQUIRED",
           "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
           "source_lock_sha256":lh,"harmonization_lock_sha256":bh,"P1C01_manifest_sha256":C01_HASH,
           "summary":s,"files":{f.name:{"bytes":f.stat().st_size,"sha256":sha(f)} for f in sorted(tmp.iterdir())},
           "limits":["OR witnesses are recorded alternatives, not complete biological realizations or equivalent biochemical states. Inclusion is of recorded sets only.","OR location support means both recorded paths supported; feature differences remain in frozen B02 witnesses; no compatibility claim beyond these flags.","AND direct required components only; no reverse transfer, nested inferred assembly or proof of sufficient inputs. Other direct requirements and B02 composition retained.","Context is factored incidence only, not causal connectivity.","Uncertain-only refers to individual state-pair witnesses; operation pairs may occur in multiple layers.","State-pair normalization is a comparison measure, not independent evidence or selected weight.","Exact overlap and layer overlap recorded; no layer aggregation or fitted model."]}
        (tmp/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared");tmp.rename(OUT)
    emit("P1C02_MANIFEST_SHA256",sha(OUT/"manifest.json"));emit("P1C02_COMPLETE_REVIEW_REQUIRED")
if __name__=="__main__":main()
