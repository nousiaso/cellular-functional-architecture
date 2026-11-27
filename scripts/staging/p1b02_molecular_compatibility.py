#!/usr/bin/env python3
"""P1B02: factored molecular compatibility and gene candidates, not reaction edges."""
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from functools import lru_cache
import csv
import json
from pathlib import Path
import re
import subprocess
import tempfile
from p1b01_source_scoped_harmonization import (
    ROOT, OUT as B01, lock_inputs, file_for, rows, sha, need, js, emit, bp
)
from p1a04p_chemical_reference_audit import chebi_id
from p1a04m_uniprot_identity_pilot import parent

OUT = ROOT / "results/P1B02_MOLECULAR_COMPATIBILITY_V1"
B01_HASH = "bb532cdca6d16e6357d06d3e76de6054275c44346c638680a1cfdbba6bcec714"
csv.field_size_limit(64 * 1024 * 1024)

def loc_profile(a,b,ga,gb):
    if not a or not b:
        return "missing_location"
    if a == b:
        return "exact_source_location"
    if len(ga) == len(gb) == 1 and ga == gb:
        return "exact_GO_location"
    if a & b or ga & gb:
        return "overlapping_location"
    return "different_recorded_location"

STRICT = {"exact_source_location", "exact_GO_location"}

def realization_index(members, edge_context):
    # One deterministic witness per terminal and path class; all paths can be
    # reconstructed from the direct composition table. Strict means every
    # step has exact recorded location agreement, not biochemical equivalence.
    @lru_cache(None)
    def walk(e):
        if not members.get(e):
            return {(e, True): (e,)}
        result = {}
        for c in sorted(members[e]):
            for (leaf, good), path in walk(c).items():
                k = (leaf, good and edge_context[e,c] in STRICT)
                candidate = (e,) + path
                if k not in result or (len(candidate),candidate) < (len(result[k]),result[k]):
                    result[k] = candidate
        return result
    return walk

def self_test():
    need(loc_profile({"x"},{"y"},{"GO:1"},{"GO:1"}) == "exact_GO_location","ontology location")
    need(loc_profile({"x"},{"y"},set(),set()) == "different_recorded_location","different location")
    need(loc_profile(set(),{"y"},set(),set()) == "missing_location","missing location")
    m={"s":{"a","t"},"t":{"a","b"}}
    ctx={("s","a"):"exact_source_location",("s","t"):"different_recorded_location",
         ("t","a"):"exact_source_location",("t","b"):"exact_source_location"}
    got=realization_index(m,ctx)("s")
    need(got=={("a",True):("s","a"),("a",False):("s","t","a"),("b",False):("s","t","b")},"OR witnesses")
    # AND is never traversed by the OR index.
    need(realization_index({}, {})("complex")=={("complex",True):("complex",)},"AND not expanded")
    emit("P1B02_SELF_TEST_PASS")

def verify_inputs():
    lock,idx,lh=lock_inputs()
    p=B01/"manifest.json"
    need(sha(p)==B01_HASH,"P1B01 manifest mismatch")
    m=json.loads(p.read_text())
    need(m["schema_version"]=="P1B01_v1" and m["source_lock_sha256"]==lh and
         m["source_manifest_sha256"]==lock["manifest_sha256"],"P1B01 lineage")
    need(m["summary"]["nodes"]==235140 and m["summary"]["merges"]==m["summary"]["scope_filters"]==0,"P1B01 summary")
    for name,x in m["files"].items():
        need(Path(name).name==name,"unsafe file")
        f=B01/name
        need(f.stat().st_size==x["bytes"] and sha(f)==x["sha256"],"P1B01 file "+name)
    emit("P1B01_MANIFEST_VERIFIED",len(m["files"]))
    return idx,lh

def build(out,idx):
    states={}; go_nodes={}; ops=set(); n=0
    for r in rows(B01/"nodes.tsv"):
        n+=1
        if r["source"]=="BP" and r["kind"]=="state":
            states[r["node_id"]]=json.loads(r["payload_json"])
        elif r["source"]=="BP" and r["kind"]=="operation":
            ops.add(r["node_id"])
        elif r["source"]=="GO" and r["kind"] in ("molecule","enabler_observation"):
            go_nodes[r["node_id"]]=(r["kind"],json.loads(r["payload_json"]))
    need((n,len(states),len(ops),len(go_nodes))==(235140,61418,16187,14707),"node counts")
    def dec(e,k): return json.loads(states[e][k])
    loc={e:set(dec(e,"location_ids_json")) for e in states}
    gloc={e:{t for _,xs in dec(e,"location_xrefs_json") for _,t in xs if re.fullmatch(r"GO:[0-9]+",t)} for e in states}
    feat={e:(set(dec(e,"feature_ids_json")),set(dec(e,"not_feature_ids_json"))) for e in states}
    members=defaultdict(set); components=defaultdict(set); direct=[]; context={}
    for r in rows(B01/"composition.tsv"):
        a,b,rel=r["owner_id"],r["child_id"],r["relation"]
        need(a in states and b in states and a!=b,"composition endpoints")
        need(rel in ("component","memberPhysicalEntity"),"composition relation")
        dest=members if rel=="memberPhysicalEntity" else components
        need(b not in dest[a],"duplicate composition")
        dest[a].add(b)
        direct.append(r)
        context[a,b]=loc_profile(loc[a],loc[b],gloc[a],gloc[b])
    need(len(direct)==76196 and len(members)==6022 and not set(members)&set(components),"composition inventory")
    need(sum(map(len,components.values()))==40083 and sum(map(len,members.values()))==36113,"AND OR inventory")
    # Confirm acyclicity over both relations, including nested complexes.
    visiting=set();done=set()
    def visit(e):
        if e in done:return
        need(e not in visiting,"composition cycle")
        visiting.add(e)
        for c in members.get(e,set())|components.get(e,set()):visit(c)
        visiting.remove(e);done.add(e)
    for e in states:visit(e)
    inputs=defaultdict(set);outputs=defaultdict(set);event_inputs=defaultdict(set)
    for r in rows(B01/"assertions.tsv"):
        if r["source"]!="BP":continue
        a,b,k=r["from_id"],r["to_id"],r["layer"]
        if k=="material_input":
            need(a in states and b in ops,"input endpoint")
            inputs[a].add(b);event_inputs[b].add(a)
        elif k=="material_output":
            need(a in ops and b in states,"output endpoint")
            outputs[b].add(a)
    def bridges(a,b):return len(outputs[a])*len(inputs[b])-len(outputs[a]&inputs[b])
    exact=sum(bridges(e,e) for e in states)
    need(exact==345409,"exact material bridge inventory")
    h=defaultdict(set)
    for r in rows(file_for(idx,"H01")):
        if r["status"]=="Approved":
            for a in r["uniprot_ids"].split("|"):
                if a.strip():h[a.strip()].add(r["hgnc_id"])
    target={}
    for label,n_expected in (("enabler",520),("participant",45)):
        inv=json.loads(file_for(idx,"UNIPROT_INVENTORY:"+label).read_text())
        rr=inv["token_records"]
        need(len(rr)==n_expected,"targeted token inventory")
        for r in rr:
            need(r["token"] not in target,"overlapping targeted token")
            target[r["token"]]=r
    need(len(target)==565,"targeted tokens")
    refs=defaultdict(set); tax=defaultdict(set); unresolved=set()
    acc_tax=defaultdict(set);acc_unknown=set()
    term_bp=defaultdict(set); bp_identifiers=defaultdict(set)
    for e,s in states.items():
        for ref,xs in dec(e,"reference_xrefs_json"):
            if s["class"]=="Protein":
                refs[ref].update(tuple(x) for x in xs)
                tax[ref].update(dec(e,"taxon_union_json"))
                if s["unresolved_taxon"]=="1":unresolved.add(ref)
            for ns,a in xs:
                t=None
                if s["class"]=="Protein" and ns in ("UniProt","UniProt Isoform"):
                    t="UniProtKB:"+a
                    acc_tax[a].update(dec(e,"taxon_union_json"))
                    if s["unresolved_taxon"]=="1":acc_unknown.add(a)
                elif s["class"]=="SmallMolecule" and ns=="ChEBI":
                    t=chebi_id(a)
                    need(t is not None,"ChEBI identifier")
                if t:
                    term_bp[t].add(e);bp_identifiers[e].add((ref,ns,a,t))
    def gene_sets(tokens):
        exact=set().union(*(h.get(x,set()) for x in tokens))
        base=set().union(*(h.get(parent(x),set()) for x in tokens))
        return exact,base
    def gene_status(exact,base):
        return ("exact_unique" if len(exact)==1 else "exact_ambiguous" if exact else
                "parent_candidate_unique" if len(base)==1 else "parent_candidate_ambiguous" if base else "no_candidate")
    with ExitStack() as stack:
        def w(name,cols):
            f=stack.enter_context((out/name).open("w",newline="",encoding="utf-8"))
            z=csv.writer(f,delimiter="\t",lineterminator="\n");z.writerow(cols.split("\t"));return z
        cw=w("composition_interfaces.tsv","owner_id\tchild_id\trelation\tstoichiometry_id\tcoefficient\tlocation_status\tfeature_id_status\tchild_to_owner_material_occurrences\towner_to_child_material_occurrences\tinterpretation")
        oc=Counter();bridge_counts=Counter();feature_counts=Counter()
        for r in sorted(direct,key=lambda x:(x["owner_id"],x["child_id"],x["relation"])):
            a,b,rel=r["owner_id"],r["child_id"],r["relation"]
            x,y=bridges(b,a),bridges(a,b)
            fs="equal_recorded_feature_ids" if feat[a]==feat[b] else "different_recorded_feature_ids"
            cw.writerow([a,b,rel,r["stoichiometry_id"],r["coefficient"],context[a,b],fs,x,y,
                         "alternative_realization" if rel=="memberPhysicalEntity" else "required_component_not_material_transfer"])
            oc[rel,context[a,b]]+=1;feature_counts[rel,fs]+=1
            bridge_counts[rel,"child_to_owner"]+=x;bridge_counts[rel,"owner_to_child"]+=y
        need(dict(bridge_counts)=={("memberPhysicalEntity","child_to_owner"):10790,
             ("memberPhysicalEntity","owner_to_child"):5392,("component","child_to_owner"):68307,
             ("component","owner_to_child"):154802},"composition bridge counts")
        rw=w("or_realizations.tsv","set_id\tterminal_id\tlocation_path_class\twitness_path_json\tdepth\tfeature_ids_differ_along_witness")
        walk=realization_index(members,context); total=0;witness_counts=Counter();max_size=0
        for e in sorted(members):
            rr=walk(e);leaves={k[0] for k in rr}
            total+=len(leaves);max_size=max(max_size,len(leaves))
            for (leaf,good),path in sorted(rr.items()):
                status="location_supported" if good else "location_uncertain"
                rw.writerow([e,leaf,status,js(path),len(path)-1,int(any(feat[a]!=feat[b] for a,b in zip(path,path[1:])))])
                witness_counts[status]+=1
        need(total==40324 and max_size==1301 and sum(any(x in members for x in v) for v in members.values())==354,"OR realization inventory")
        # Record explicit material reactions that produce a complex. A complete
        # direct component-input match is a candidate assembly witness, not
        # stoichiometric or biochemical proof. No reverse dissociation inferred.
        aw=w("complex_production_candidates.tsv","complex_id\toperation_id\tmatched_direct_components_json\tunmatched_direct_components_json\tstatus")
        assembly=Counter()
        for e,cs in sorted(components.items()):
            for r in sorted(outputs[e]):
                found=cs&event_inputs[r];missing=cs-found
                tag="all_direct_components_recorded_as_inputs" if not missing else "partial_direct_component_inputs" if found else "no_direct_component_input_match"
                aw.writerow([e,r,js(sorted(found)),js(sorted(missing)),tag]);assembly[tag]+=1
        # Factored identifier memberships avoid cross-products between sources.
        iw=w("identifier_memberships.tsv","node_id\tkind\toriginal_reference_or_path\toriginal_namespace\toriginal_identifier\tshared_identifier")
        for e,rr in sorted(bp_identifiers.items()):
            for ref,ns,a,t in sorted(rr):iw.writerow([e,"BP_state",ref,ns,a,t])
        gw=w("gene_candidates.tsv","source_kind\treference_or_token\taccessions_json\ttaxonomy_status\texact_HGNC_json\tparent_HGNC_json\ttargeted_UniProt_HGNC_json\tstatus")
        gc=Counter()
        for ref,xs in sorted(refs.items()):
            tokens={a for ns,a in xs if ns in ("UniProt","UniProt Isoform")}
            tx="human" if tax[ref]=={"9606"} and ref not in unresolved else "nonhuman_or_uncertain"
            exact_g,base=gene_sets(tokens) if tx=="human" else (set(),set())
            tag=gene_status(exact_g,base) if tx=="human" else "not_assigned"
            gw.writerow(["BP_reference",bp(ref),js(sorted(tokens)),tx,js(sorted(exact_g)),js(sorted(base)),"[]",tag])
            if tx=="human":gc[tag]+=1
        need(dict(gc)=={"exact_unique":11308,"parent_candidate_unique":310,"exact_ambiguous":30,"no_candidate":34},"BP human gene mapping")
        def terms(obj,path="$"):
            if isinstance(obj,dict):
                t=obj.get("term")
                if isinstance(t,str):yield path,t
                yield from terms(obj.get("has_part",[]),path+".has_part")
            elif isinstance(obj,list):
                for j,x in enumerate(obj):yield from terms(x,path+"["+str(j)+"]")
        go_tokens=set();mol_counts=Counter();chem_profile=Counter()
        cross=w("cross_source_molecule_candidates.tsv","go_node_id\tidentifier\tgo_direct_location\tBP_exact_location_ids_json\tBP_missing_location_ids_json\tBP_other_location_ids_json\tstatus")
        for e,(kind,obj) in sorted(go_nodes.items()):
            tt=list(terms(obj))
            for path,t in tt:
                if t.startswith("UniProtKB:"):go_tokens.add(t.split(":",1)[1])
                ns,_,a=t.partition(":")
                iw.writerow([e,"GO_"+kind,path,ns,a,t])
            if kind!="molecule":continue
            t=obj["term"];ns=t.split(":",1)[0];mol_counts[ns]+=1
            if ns not in ("CHEBI","UniProtKB"):continue
            q=obj.get("located_in") or {}
            lt=q.get("term") if isinstance(q,dict) else None
            es=term_bp.get(t,set())
            same=sorted(x for x in es if lt and lt in gloc[x])
            missing=sorted(x for x in es if not lt or not gloc[x])
            other=sorted(es-set(same)-set(missing))
            tag=("no_exact_BP_identifier" if not es else "GO_location_unrecorded" if not lt else
                 "exact_direct_GO_location_match" if same else "no_exact_location_match")
            cross.writerow([e,t,lt or "",js(same),js(missing),js(other),tag])
            if ns=="CHEBI":chem_profile[tag]+=1
        need(mol_counts["CHEBI"]==3346 and mol_counts["UniProtKB"]==3864 and sum(mol_counts.values())==7322,"GO molecular inventory")
        need(dict(chem_profile)=={"no_exact_BP_identifier":503,"GO_location_unrecorded":2550,"exact_direct_GO_location_match":149,"no_exact_location_match":144},"chemical correspondence inventory")
        go_gene=Counter()
        for a in sorted(go_tokens):
            rec=target.get(a)
            if rec is not None:
                tx="human" if rec["taxon_id"]==9606 else "nonhuman" if rec["taxon_id"] is not None else "unresolved"
                source="targeted_"+rec["identity"]+"_"+rec["substate"]
                verified=set(rec.get("HGNC_xrefs",[])) if tx=="human" else set()
                # Use returned parent only as a gene candidate, not a state identity.
                exact_g=h.get(a,set()) if tx=="human" else set()
                base=h.get(rec["primary"],set()) if tx=="human" else set()
            else:
                ts=acc_tax.get(a,set())
                tx="human" if ts=={"9606"} and a not in acc_unknown else "nonhuman" if ts and "9606" not in ts and a not in acc_unknown else "unresolved"
                source="BP_exact_taxonomy_observation";verified=set()
                exact_g,base=gene_sets({a}) if tx=="human" else (set(),set())
            tag=gene_status(exact_g,base) if tx=="human" else "not_assigned"
            if rec is not None and tx=="human" and base!=verified:tag="HGNC_parent_disagreement_review"
            gw.writerow(["GO_token",a,js([a]),source+":"+tx,js(sorted(exact_g)),js(sorted(base)),js(sorted(verified)),tag])
            go_gene[tag]+=1
    summary={"exact_material_occurrences":exact,"direct_composition":len(direct),
        "composition_location_profiles":[[list(k),n] for k,n in sorted(oc.items())],
        "composition_feature_profiles":[[list(k),n] for k,n in sorted(feature_counts.items())],
        "direct_bridge_occurrences":[[list(k),n] for k,n in sorted(bridge_counts.items())],
        "OR_owners":len(members),"OR_unique_terminal_pairs":total,"OR_max_realizations":max_size,
        "OR_witness_rows":dict(witness_counts),"complex_production_profiles":dict(assembly),
        "BP_human_reference_gene_candidates":dict(gc),"GO_token_gene_candidates":dict(go_gene),
        "GO_unique_UniProt_tokens":len(go_tokens),"GO_chemical_correspondence":dict(chem_profile),
        "graph_edges_created":0,"merges":0,"filters":0}
    emit("P1B02_COUNTS",summary)
    return summary

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true");ap.add_argument("--execute",action="store_true")
    a=ap.parse_args();self_test()
    if a.self_test:return
    idx,lh=verify_inputs()
    if not a.execute:
        emit("PLAN_ONLY",OUT.relative_to(ROOT));return
    need(not OUT.exists(),"output exists; no overwrite")
    with tempfile.TemporaryDirectory(prefix=".p1b02_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir()
        summary=build(tmp,idx)
        m={"schema_version":"P1B02_v1","git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
           "source_lock_sha256":lh,"P1B01_manifest_sha256":B01_HASH,"summary":summary,
           "files":{p.name:{"bytes":p.stat().st_size,"sha256":sha(p)} for p in sorted(tmp.iterdir())},
           "limits":["OR witness location agreement not biochemical equivalence","AND prerequisites not material edges","Exact identifier/location candidates not molecular state merges","Protein gene candidates not activity or complex sufficiency","Complex production candidates not verified assembly or dissociation"]}
        (tmp/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared");tmp.rename(OUT)
    emit("P1B02_MANIFEST_SHA256",sha(OUT/"manifest.json"))
    emit("P1B02_COMPLETE; NO_GRAPH_EDGES_OR_MERGES")
if __name__=="__main__":main()
