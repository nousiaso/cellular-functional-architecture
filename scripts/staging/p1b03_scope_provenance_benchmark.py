#!/usr/bin/env python3
"""P1B03: scope, native pathway benchmarks and provenance; no filtering or inference."""
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
import xml.etree.ElementTree as ET
import zipfile

from p1b01_source_scoped_harmonization import ROOT, OUT as B01, rows, sha, need, js, emit, bp, file_for
from p1b02_molecular_compatibility import OUT as B02, B01_HASH, verify_inputs
from p1a04k_scope_feasibility_audit import all_of, any_of
from p1a04e_ontology_references import parse_obo
from p1a04j_provenance_audit import norm_pmid

OUT=ROOT/"results/P1B03_SCOPE_PROVENANCE_BENCHMARK_V1"
B02_HASH="ba2ee290a9aeac833721b72a321273f8be64dd54dd0c47fbbcd7ce4b6d0cbbec"
RDF="{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
XML="{http://www.w3.org/XML/1998/namespace}"
BP="{http://www.biopax.org/release/biopax-level3.owl#}"
OWL_ONTOLOGY="{http://www.w3.org/2002/07/owl#}Ontology"
STABLE=re.compile(r"R-[A-Z]+-[0-9]+")
EVENT={"BiochemicalReaction","Degradation","TemplateReaction"}
CONTROL={"Catalysis","Control","TemplateReactionRegulation"}

def references(x):
    if isinstance(x,dict):
        for k,v in x.items():
            if k=="reference" and isinstance(v,str):yield v
            yield from references(v)
    elif isinstance(x,list):
        for v in x:yield from references(v)

def evidence_terms(x):
    if isinstance(x,list):
        return {str(v.get("term","missing")) for v in x if isinstance(v,dict)}
    return set()

def validate_id(a,tag,ordinal,rid,about,classes):
    # The frozen source contains one empty rdf:about on its initial OWL
    # ontology declaration. Preserve this record and its ordinal.
    metadata=(a=="" and ordinal==1 and tag==OWL_ONTOLOGY and rid is None and about=="")
    need((bool(a) or metadata) and a not in classes,
         "duplicate/missing BP ID ordinal="+str(ordinal)+" tag="+tag+" id="+repr(a))
    return metadata

def self_test():
    need(validate_id("",OWL_ONTOLOGY,1,None,"",{}),"empty ontology metadata")
    need(not validate_id("Protein1",BP+"Protein",2,"Protein1",None,{}),"ordinary BP ID")
    for args in (("",BP+"Protein",2,None,"",{}),
                 (None,BP+"Protein",2,None,None,{}),
                 ("Protein1",BP+"Protein",3,"Protein1",None,{"Protein1":"Protein"}),
                 ("",OWL_ONTOLOGY,2,None,"",{})):
        try:validate_id(*args)
        except RuntimeError:pass
        else:raise RuntimeError("invalid BP ID accepted in fixture")
    need(all_of(["H","F"])=="F" and any_of(["F","H"])=="H","composition logic")
    need(list(references({"a":[{"reference":"PMID: 1"}]}))==["PMID: 1"],"nested references")
    need(norm_pmid("PMID: 1")=="1","PMID")
    need(evidence_terms([{"term":"ECO:0000313"}])=={"ECO:0000313"},"ECO")
    emit("P1B03_SELF_TEST_PASS")

def verify():
    idx,lh=verify_inputs()
    p=B02/"manifest.json";need(sha(p)==B02_HASH,"B02 manifest")
    m=json.loads(p.read_text())
    need(m["schema_version"]=="P1B02_v1" and m["P1B01_manifest_sha256"]==B01_HASH and m["source_lock_sha256"]==lh,"B02 lineage")
    for name,x in m["files"].items():
        need(Path(name).name==name,"unsafe filename")
        f=B02/name;need(f.stat().st_size==x["bytes"] and sha(f)==x["sha256"],"B02 file "+name)
    emit("P1B02_MANIFEST_VERIFIED",len(m["files"]))
    return idx,lh

def scan(path):
    classes={}; links=[]; xrefs={}; owners=defaultdict(list); names={}
    stack=[];base=None;ordinal=0
    def norm(s):
        if s is None:return None
        if s.startswith("#"):return s[1:]
        if base and s.startswith(base.rstrip("#")+"#"):return s[len(base.rstrip("#"))+1:]
        return s
    with zipfile.ZipFile(path) as z,z.open("Homo_sapiens.owl") as f:
        for ev,e in ET.iterparse(f,events=("start","end")):
            if ev=="start":
                stack.append(e)
                if len(stack)==1:
                    need(e.tag==RDF+"RDF","XML root");base=e.get(XML+"base")
                continue
            if len(stack)==2:
                ordinal+=1
                a=norm(e.get(RDF+"ID") or e.get(RDF+"about"));cl=e.tag.split("}")[-1]
                if validate_id(a,e.tag,ordinal,e.get(RDF+"ID"),e.get(RDF+"about"),classes):
                    emit("BIOPAX_EMPTY_ONTOLOGY_METADATA_RETAINED",ordinal)
                classes[a]=cl
                if cl in ("UnificationXref","PublicationXref"):
                    xrefs[a]=(cl,[c.text for c in e if c.tag==BP+"db"],[c.text for c in e if c.tag==BP+"id"])
                if cl in EVENT|CONTROL|{"Pathway"}:
                    for j,c in enumerate(e):
                        if c.tag==BP+"xref":
                            owners[a].append((norm(c.get(RDF+"resource")),ordinal,j))
                    if cl=="Pathway":names[a]=[c.text for c in e if c.tag==BP+"displayName"]
                for j,c in enumerate(e):
                    if c.tag in (BP+"pathwayComponent",BP+"pathwayOrder",BP+"stepProcess",BP+"nextStep"):
                        links.append((a,norm(c.get(RDF+"resource")),c.tag.split("}")[-1],ordinal,j))
                stack[-2].remove(e);e.clear()
                if ordinal%100000==0:emit("BIOPAX_SEQUENTIAL_RECORDS",ordinal)
            stack.pop()
    need(len(classes)==548688 and classes.get("")=="Ontology","BP record count/ontology metadata")
    return classes,links,xrefs,owners,names

def build(out,idx):
    states={};ops={};goact={}
    for r in rows(B01/"nodes.tsv"):
        if r["source"]=="BP" and r["kind"] in ("state","operation"):
            (states if r["kind"]=="state" else ops)[r["local_id"]]=json.loads(r["payload_json"])
        elif r["source"]=="GO" and r["kind"]=="activity":
            goact[r["node_id"]]=(r["model_id"],r["local_id"],json.loads(r["payload_json"]))
    need((len(states),len(ops),len(goact))==(61418,16187,7385),"catalogue counts")
    stable={x["stable_id"]:e for e,x in ops.items()}
    need(len(stable)==16187,"stable event uniqueness")
    classes,links,xrefs,owners,names=scan(file_for(idx,"R01"))
    pathways=set(names)
    need(len(pathways)==2889,"pathway count")
    def ids(e):
        return {a for x,_,_ in owners[e] if x in xrefs and xrefs[x][0]=="UnificationXref" and "Reactome" in xrefs[x][1] for a in xrefs[x][2] if a and STABLE.fullmatch(a)}
    need(all(ids(e)=={x["stable_id"]} for e,x in ops.items()),"native event stable IDs")
    pid={e:ids(e) for e in pathways}
    children=defaultdict(set);parents=defaultdict(set);direct=defaultdict(set);stepowners=defaultdict(set);stepops=defaultdict(set)
    path_counts=Counter()
    for a,b,k,_,_ in links:
        need(a in classes and b in classes,"unresolved native path endpoint")
        if k in ("pathwayComponent","pathwayOrder"):
            path_counts[k,classes[a],classes[b]]+=1
        if k=="pathwayComponent" and a in pathways:
            if b in pathways:children[a].add(b);parents[b].add(a)
            elif b in ops:direct[b].add(a)
        elif k=="pathwayOrder" and a in pathways:stepowners[b].add(a)
        elif k=="stepProcess" and b in ops:stepops[a].add(b)
    need(dict(path_counts)=={("pathwayComponent","Pathway","BiochemicalReaction"):17063,
        ("pathwayComponent","Pathway","Degradation"):20,("pathwayComponent","Pathway","Pathway"):2905,
        ("pathwayComponent","Pathway","TemplateReaction"):8,("pathwayOrder","Pathway","PathwayStep"):19996},"pathway predicate counts")
    visiting=set()
    @lru_cache(None)
    def ancestors(e):
        need(e not in visiting,"pathway hierarchy cycle")
        visiting.add(e)
        try:return frozenset(parents[e]).union(*(ancestors(a) for a in parents[e]))
        finally:visiting.remove(e)
    for e in pathways:ancestors(e)
    step_members=defaultdict(set)
    for s,ee in stepops.items():
        for e in ee:step_members[e].update(stepowners[s])
    r04=defaultdict(set)
    for x in rows(file_for(idx,"R04")):
        if x["reaction_id"] in stable:r04[stable[x["reaction_id"]]].add(x["pathway_id"])
    need(len(r04)==15771,"R04 event coverage")
    def pairfile(id):
        with file_for(idx,id).open(newline="",encoding="utf-8") as f:
            for r in csv.reader(f,delimiter="\t"):
                if r:
                    need(len(r)==2,"pair export schema "+id)
                    yield tuple(r)
    exported_hierarchy={(a,b) for a,b in pairfile("R09") if a.startswith("R-HSA-") and b.startswith("R-HSA-")}
    need(len(exported_hierarchy)==2905,"R09 human hierarchy")
    native_hierarchy={(next(iter(pid[a])),next(iter(pid[b]))) for a in children for b in children[a] if len(pid[a])==len(pid[b])==1}
    disease={a for a,b in pairfile("R11") if STABLE.fullmatch(a)}
    # Only assign disease labels using exact stable IDs; incomplete native
    # identity is explicitly flagged, not assumed to be nondisease.
    r06=defaultdict(set)
    for a,b in pairfile("R06"):
        if a in stable and re.fullmatch(r"[0-9]+",b):r06[stable[a]].add(b)
    need(len(r06)==14821,"R06 numeric coverage")
    comp=defaultdict(list);members=defaultdict(list)
    for r in rows(B01/"composition.tsv"):
        a,b=r["owner_id"][3:],r["child_id"][3:]
        (comp if r["relation"]=="component" else members)[a].append(b)
    @lru_cache(None)
    def scope(e):
        s=states[e]
        if comp.get(e) and members.get(e):return "U"
        if comp.get(e):return all_of(scope(x) for x in comp[e])
        if members.get(e):return any_of(scope(x) for x in members[e])
        t=set(json.loads(s["taxon_union_json"]))
        if s["class"]=="SmallMolecule":return "H"
        if s["unresolved_taxon"]=="1" or not t or ("9606" in t and len(t)>1):return "U"
        return "H" if t=={"9606"} else "F"
    @lru_cache(None)
    def human_protein_leaf(e):
        cs=comp.get(e,[])+members.get(e,[])
        return any(human_protein_leaf(x) for x in cs) if cs else states[e]["class"]=="Protein" and scope(e)=="H"
    state_counts=Counter(scope(e) for e in states)
    need(dict(state_counts)=={"H":55616,"F":2527,"U":3275},"state scope")
    roles=defaultdict(lambda:defaultdict(set))
    for r in rows(B01/"participant_coefficients.tsv"):
        roles[r["event_id"][3:]][r["original_role"]].add(r["entity_id"][3:])
    for r in rows(B01/"direct_control_candidates.tsv"):
        roles[r["event_id"][3:]][r["role"]].add(r["controller_id"][3:])
    def profile(es):return "+".join(sorted({scope(e) for e in es})) if es else "none_recorded"
    _,eco=parse_obo(file_for(idx,"ECO").read_bytes())
    @lru_cache(None)
    def eco_import(e):
        # Follow pinned is_a ancestry only; no free-text classification of citations.
        if e not in eco:return False
        if e=="ECO:0000313" or "imported information" in " ".join(eco[e].get("name",[])).lower():return True
        return any(eco_import(x.split()[0]) for x in eco[e].get("is_a",[]))
    bp_pmids=defaultdict(set);pub_counts=Counter();event_pub_pairs=set();event_pub_ns=Counter()
    with ExitStack() as stack:
        def w(name,cols):
            f=stack.enter_context((out/name).open("w",encoding="utf-8",newline=""))
            a=csv.writer(f,delimiter="\t",lineterminator="\n");a.writerow(cols.split("\t"));return a
        pw=w("bp_publication_links.tsv","owner_id\towner_class\txref_id\trecord_ordinal\tproperty_ordinal\tnamespaces_json\tidentifiers_json\tpublication_identity_status")
        for e in sorted(owners):
            if classes[e] not in EVENT|CONTROL:continue
            for x,ord_,j in owners[e]:
                if x not in xrefs or xrefs[x][0]!="PublicationXref":continue
                _,ns,ii=xrefs[x]
                kind="event" if e in ops else "control"
                status=("namespace_and_identifier_recorded" if ns and ii else
                        "partial_identity_fields" if ns or ii else "missing_namespace_and_identifier")
                pw.writerow([bp(e),classes[e],bp(x),ord_,j,js(ns),js(ii),status])
                pub_counts[kind]+=1
                pub_counts[kind+"_"+status]+=1
                if e in ops:
                    event_pub_pairs.add((e,x))
                    event_pub_ns[tuple(sorted(map(str,ns)))]+=1
                if any(isinstance(a,str) and a.lower()=="pubmed" for a in ns):
                    bp_pmids[e].update(a for a in ii if isinstance(a,str) and re.fullmatch(r"[0-9]+",a))
        emit("BP_PUBLICATION_COUNTS",dict(pub_counts),"unique_event_xref_pairs",len(event_pub_pairs),
             "event_namespaces",sorted(event_pub_ns.items()))
        need(pub_counts["event"]==len(event_pub_pairs)==60015,"native event publication occurrence/pair count")
        need(pub_counts["event_namespace_and_identifier_recorded"]==59951 and
             pub_counts["event_missing_namespace_and_identifier"]==64 and
             pub_counts["event_partial_identity_fields"]==0,"publication identity coverage")
        need(dict(event_pub_ns)=={("Pubmed",):59811,("ISBN",):140,():64},"event publication namespace profile")
        bp_pmids={e:s for e,s in bp_pmids.items() if s}
        evpub={e:s for e,s in bp_pmids.items() if e in ops}
        need(len(evpub)==15701 and all(r06[e]<=evpub.get(e,set()) for e in r06),"BP R06 coverage")
        all_bp=set().union(*bp_pmids.values())
        wpath=w("pathways.tsv","pathway_id\tstable_id_candidates_json\tname_json\tparent_ids_json\troot_ids_json")
        for e in sorted(pathways):
            anc=ancestors(e);roots={a for a in anc|{e} if not parents[a]}
            wpath.writerow([bp(e),js(sorted(pid[e])),js(names[e]),js(sorted(map(bp,parents[e]))),js(sorted(map(bp,roots)))])
        wl=w("native_pathway_assertions.tsv","from_id\tto_id\tpredicate\trecord_ordinal\tproperty_ordinal")
        for a,b,k,o,j in links:wl.writerow([bp(a),bp(b),k,o,j])
        wm=w("event_pathway_benchmark.tsv","operation_id\tstable_id\tnative_direct_json\tnative_ancestor_json\tnative_root_json\tstep_process_candidate_json\tR04_stable_pathways_json\tnative_direct_stable_json\tR04_comparison_status\tnative_path_identity_unresolved\tR04_disease_context")
        cov=Counter();dis=Counter()
        for e in sorted(ops):
            ds=direct[e];anc=set().union(*(ancestors(a) for a in ds))
            roots={a for a in ds|anc if not parents[a]}
            ss=set().union(*(pid[a] for a in ds if len(pid[a])==1))
            rr=r04[e];unresolved=any(len(pid[a])!=1 for a in ds|anc)
            tag="equal" if ss==rr else "R04_subset" if rr<ss else "native_subset" if ss<rr else "partial_overlap" if ss&rr else "disjoint"
            d="no_R04_mapping" if not rr else "disease_only" if rr<=disease else "mixed" if rr&disease else "no_listed_disease_context"
            wm.writerow([bp(e),ops[e]["stable_id"],js(sorted(map(bp,ds))),js(sorted(map(bp,anc))),js(sorted(map(bp,roots))),js(sorted(map(bp,step_members[e]))),js(sorted(rr)),js(sorted(ss)),tag,int(unresolved),d])
            cov["native_mapped"]+=bool(ds);cov["step_mapped"]+=bool(step_members[e])
            cov["native_rescues_R04_missing"]+=bool(ds and not rr);cov["comparison_"+tag]+=1;dis[d]+=1
        need(dict(dis)=={"no_R04_mapping":416,"disease_only":2052,"mixed":95,"no_listed_disease_context":13624},"R04 disease context")
        ws=w("state_scope.tsv","state_id\tcomposition_feasibility\trecorded_human_protein_leaf\tclass\ttaxon_union_json\tunresolved_taxon")
        for e,s in sorted(states.items()):
            ws.writerow([bp(e),scope(e),int(human_protein_leaf(e)),s["class"],s["taxon_union_json"],s["unresolved_taxon"]])
        we=w("event_scope.tsv","operation_id\tcatalyst_profile\tpositive_regulator_profile\tnegative_regulator_profile\tinput_profile\toutput_profile\tcatalyst_human_protein_leaf\tforeign_required_catalyst_with_human_leaf\tcandidate_scope_stratum\tno_material_output_recorded")
        scope_counts=Counter();cat_counts=Counter();foreign_human=0
        for e in sorted(ops):
            rr=roles[e];cat=rr["catalyst"];inp=rr["left"]|rr["template"];op=rr["right"]|rr["product"]
            cp=profile(cat); hp=any(human_protein_leaf(x) for x in cat)
            fh=any(scope(x)=="F" and human_protein_leaf(x) for x in cat)
            if not cat:tag="no_catalyst_recorded"
            elif cp=="H":tag="human_compatible_catalyst"
            elif cp=="F":tag="foreign_required_catalyst_with_human_leaf" if fh else "foreign_required_catalyst_no_recorded_human_leaf"
            else:tag="mixed_or_unresolved_catalyst"
            we.writerow([bp(e),cp,profile(rr["positiveregulator"]),profile(rr["negativeregulator"]),profile(inp),profile(op),int(hp),int(fh),tag,int(not op)])
            scope_counts[tag]+=1;cat_counts[cp]+=1;foreign_human+=int(fh and cp=="F")
        need(dict(cat_counts)=={"H":6092,"F":256,"U":516,"none_recorded":9323},"catalyst scope")
        need(foreign_human==41,"foreign-required with human leaf")
        # Factored evidence memberships, not assumed independent blocks.
        # Shared PMIDs may generate huge connected components; assess later.
        wg=w("go_activity_provenance.tsv","activity_id\tmodel_id\tstable_event_candidate\tReactome_event_citations_json\tall_references_json\tPMIDs_json\tBP_publication_overlap_json\tdependency_flag")
        act_flags={};act_pmids={};ac=Counter()
        for e,(mid,sid,obj) in sorted(goact.items()):
            refs=set(references(obj));ps={x for s in refs if (x:=norm_pmid(s))}
            rc={s.strip().split(":",1)[1] for s in refs if s.strip().startswith("Reactome:")}
            raw=sid.removeprefix("gomodel:");match=raw if raw in stable else ""
            flag=bool(mid.startswith("R-HSA-") or rc or match)
            act_flags[e]=flag;act_pmids[e]=ps
            wg.writerow([e,mid,match,js(sorted(rc)),js(sorted(refs)),js(sorted(ps)),js(sorted(ps&all_bp)),int(flag)])
            ac["bare_match"]+=bool(match);ac["dependency_flag"]+=flag
        need(ac["bare_match"]==682,"GO exact activity crosswalk")
        wa=w("go_assertion_provenance.tsv","assertion_id\tmodel_id\tlayer\tECO_terms_json\treferences_json\tPMIDs_json\tBP_publication_overlap_json\tendpoint_activity_PMIDs_overlap_json\timport_marker\tautomatic_inference_marker\tendpoint_or_model_Reactome_dependency\tcorroboration_stratum")
        ee=Counter();et=Counter()
        for r in rows(B01/"assertions.tsv"):
            if r["source"]!="GO":continue
            pr=json.loads(r["provenance_json"]);mid=pr["model_id"]
            ev=json.loads(r["evidence_json"]);terms=evidence_terms(ev);et.update(terms)
            refs=set(references(ev));ps={x for s in refs if (x:=norm_pmid(s))}
            imp=any(eco_import(t) for t in terms)
            auto="ECO:0000312" in terms or any("automatic assertion" in " ".join(eco.get(t,{}).get("name",[])).lower() for t in terms)
            endpoints=[a for a in (r["from_id"],r["to_id"]) if a in goact]
            dep=mid.startswith("R-HSA-") or any(act_flags[a] for a in endpoints) or any(s.strip().startswith("Reactome:") for s in refs)
            shared=ps&all_bp;ap=set().union(*(act_pmids[a] for a in endpoints))&all_bp
            if dep or imp:tag="known_or_flagged_Reactome_or_import_dependence"
            elif not terms or not refs:tag="insufficient_evidence_metadata"
            elif auto:tag="automatic_inference_review"
            elif not ps:tag="no_recognized_PMID_review"
            elif shared or ap:tag="publication_overlap_review"
            else:tag="candidate_corroboration_no_observed_dependence_marker"
            wa.writerow([r["assertion_id"],mid,r["layer"],js(sorted(terms)),js(sorted(refs)),js(sorted(ps)),js(sorted(shared)),js(sorted(ap)),int(imp),int(auto),int(dep),tag]);ee[tag]+=1
        need(sum(ee.values())==15491 and len(et)==22,"GO assertion inventory")
        wt=w("go_evidence_terms.tsv","ECO_id\tname_json\tis_a_json\timport_ancestry_flag\tassertion_occurrences")
        for t,n in sorted(et.items()):
            need(t in eco,"unresolved ECO term "+t)
            wt.writerow([t,js(eco[t].get("name",[])),js(eco[t].get("is_a",[])),int(eco_import(t)),n])
    summary={"native_pathway_counts":[[list(k),n] for k,n in sorted(path_counts.items())],
        "pathway_stable_identity_cardinality":dict(Counter(str(len(x)) for x in pid.values())),
        "native_R09_hierarchy_equal_pairs":len(native_hierarchy&exported_hierarchy),
        "native_hierarchy_only":len(native_hierarchy-exported_hierarchy),"R09_only":len(exported_hierarchy-native_hierarchy),
        "pathway_coverage":dict(cov),"R04_disease_context":dict(dis),"BP_publication_links":dict(pub_counts),"BP_event_unique_publication_pairs":len(event_pub_pairs),
        "BP_events_with_PMIDs":len(evpub),"R06_numeric_events":len(r06),"state_scope":dict(state_counts),
        "event_catalyst_profiles":dict(cat_counts),"event_scope_strata":dict(scope_counts),
        "foreign_required_catalyst_with_human_leaf":foreign_human,"GO_activity_profiles":dict(ac),
        "GO_assertion_strata":dict(ee),"GO_ECO_counts":dict(et),"scope_filters":0,"graph_edges_created":0}
    emit("P1B03_COUNTS",summary);return summary

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--execute",action="store_true");ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args();self_test()
    if a.self_test:return
    idx,lh=verify()
    if not a.execute:emit("PLAN_ONLY",OUT.relative_to(ROOT));return
    need(not OUT.exists(),"output exists; no overwrite")
    with tempfile.TemporaryDirectory(prefix=".p1b03_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir()
        s=build(tmp,idx)
        m={"schema_version":"P1B03_v1","git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
           "source_lock_sha256":lh,"P1B01_manifest_sha256":B01_HASH,"P1B02_manifest_sha256":B02_HASH,"summary":s,
           "files":{p.name:{"bytes":p.stat().st_size,"sha256":sha(p)} for p in sorted(tmp.iterdir())},
           "limits":["Native pathway, ancestor and step memberships are separate benchmark records; no primary fitting","H/F/U is recorded compositional feasibility, not activity; no event exclusion","GO absence of dependence markers does not prove nonimport or independent experiments","Evidence memberships factored; final leakage-safe splits and domain coverage tests require graph feasibility","No inferred pathway memberships from shared genes"]}
        (tmp/"manifest.json").write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared");tmp.rename(OUT)
    emit("P1B03_MANIFEST_SHA256",sha(OUT/"manifest.json"));emit("P1B03_COMPLETE; NO_FILTERS_OR_GRAPH_EDGES")
if __name__=="__main__":main()
