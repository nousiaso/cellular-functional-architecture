#!/usr/bin/env python3
"""Create a review candidate manifest from actual verified local receipts, not a freeze."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from p1a04s_remaining_identity_and_format_preflight import core,ROOT,STATE,DB,rows,sha,need,PART_HASH
from p1a04k_scope_feasibility_audit import TYPED,check_manifest
from p1a04q_state_correspondence_audit import CHEBI_HASH,UNI_HASH
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"acquisition"))
import p1a04n_acquire_targeted_uniprot as u
import p1a04r_remaining_participant_identities as r
RAW=ROOT/"data/raw"
OUT=ROOT/"results/P1A05_SOURCE_LOCK_CANDIDATE_V1"
GO_HASH="03418cc71f0195d49f89f3bcbad09ca1e6416867c58afa0d0fe48304a4034d09"
DB_HASH="9fc75c327a76de4ceb82c9095e9f5c82e34c6a202b1468f0836a53a498cdd63b"
def emit(*x):print(*x,flush=True)
def main():
    need(not OUT.exists(),"candidate output already exists; no overwrite")
    check_manifest(TYPED,"P1A04F_v1");check_manifest(STATE,"P1A04H_v1")
    result=[];seen=set()
    def add(id,resource,release,role,path,url,utc,license,status,expected=None,receipt=None,payload_hash=None,payload_bytes=None):
        path=path.resolve()
        need(path.is_relative_to(ROOT.resolve()) and path.is_file(),"missing or unsafe path "+id)
        rel=str(path.relative_to(ROOT.resolve()))
        need(rel not in seen,"duplicate path "+rel);seen.add(rel)
        h=sha(path);n=path.stat().st_size
        if expected:need(h==expected,"hash mismatch "+id)
        rec_rel="";rec_hash=""
        if receipt is not None:
            receipt=receipt.resolve()
            need(receipt.is_relative_to(ROOT.resolve()),"unsafe receipt")
            rec_rel=str(receipt.relative_to(ROOT.resolve()));rec_hash=sha(receipt)
        result.append(dict(id=id,resource=resource,release=release,role=role,url=url,local_relative_path=rel,bytes=n,sha256=h,acquired_utc=utc or "not_recorded",license=license,status=status,receipt_path=rec_rel,receipt_sha256=rec_hash,original_payload_sha256=payload_hash or h,original_payload_bytes=n if payload_bytes is None else payload_bytes))
    cp=RAW/"P1A03_CORE_RECEIPTS.json";cr=json.loads(cp.read_text())
    need(len(cr)==19 and len({x["id"] for x in cr})==19,"core receipt count")
    pins={x["id"]:x for x in rows(ROOT/"manifests/P1A03_CORE_REPORTED_SHA256.tsv")}
    need(set(pins)=={x["id"] for x in cr},"core receipt IDs")
    for x in cr:
        id=x["id"];p=(ROOT/x["path"]).resolve();side=p.with_name(p.name+".receipt.json")
        need(p.is_relative_to(RAW.resolve()),"core path")
        need(json.loads(side.read_text())==x and p.stat().st_size==x["bytes"]==int(pins[id]["bytes"]) and x["sha256"]==pins[id]["sha256"],"core receipt mismatch "+id)
        lic="CC-BY-4.0" if x["resource"]=="GO" or id=="R02" else "CC0-1.0"
        status="DEFERRED_SEMANTICS_BYTE_VERIFIED" if id=="R02" else "VERIFIED_WITHIN_AUDITED_ROLE"
        add(id,x["resource"],x["release_candidate"],x["role"],p,x["url"],x["acquired_utc"],lic,status,x["sha256"],side)
    add("CORE_RECEIPTS","CFA","P1A03","acquisition_inventory",cp,"","", "internal_metadata","RECEIPT_INVENTORY")
    emit("CORE_MANIFEST_ROWS",19)
    gp=RAW/"go_2026-08-05/human_indexed_models/MODEL_INVENTORY.json"
    need(sha(gp)==GO_HASH,"GO inventory hash")
    gs=json.loads(gp.read_text());need(len(gs)==1520,"GO model count")
    ids=set()
    for x in gs:
        id=x["model"]["model_id"]
        need(id not in ids and Path(id).name==id,"GO model ID");ids.add(id)
        path=gp.parent/(id+".json");side=gp.parent/(id+".receipt.json")
        need(json.loads(side.read_text())==x and path.stat().st_size==x["bytes"],"GO receipt mismatch "+id)
        add("GO_MODEL:"+id,"GO-CAM",x["release_archive"],"source_scoped_primary_candidate",path,x["url"],x["acquired_utc"],"CC-BY-4.0","VERIFIED_HUMAN_INDEXED_NOT_SCOPE_FILTERED",x["sha256"],side)
    add("GO_INVENTORY","GO-CAM","2026-08-05","acquisition_inventory",gp,"","","CC-BY-4.0","PINNED",GO_HASH)
    emit("GO_MANIFEST_ROWS",len(ids))
    for name,release,f,h in (("RO","2026-09-04","ro.obo","f2145eaecd649b395f5a908906f429fd6ac6dca2bc6816045b5c8a80da9a6b10"),("ECO","2026-07-10","eco.obo","c8b52a2da2b0f92224b50a1d8773fd7e3bca0977d920d4e31ef9fde960c84fd4")):
        path=RAW/"ontology_interpretation"/name/release/f;side=path.with_name(path.name+".receipt.json")
        x=json.loads(side.read_text());need(x["sha256"]==h and path.stat().st_size==x["bytes"] and x["release"]==release,"ontology receipt")
        add(name,name,release,"interpretation_reference",path,x["url"],x.get("acquired_utc"),"CC0-1.0","PINNED",h,side)
    path=RAW/"ontology_interpretation/ChEBI/monthly_snapshot/chebi.obo.gz";side=path.parent/"receipt.json";x=json.loads(side.read_text())
    need(x["sha256"]==CHEBI_HASH and x["data_version"]=="data-version: 257" and x["bytes"]==path.stat().st_size,"ChEBI receipt")
    add("ChEBI","ChEBI","257","chemical_reference",path,x["url"],x["acquired_utc"],"CC-BY-4.0","PINNED",CHEBI_HASH,side)
    parents=set()
    for label,base,h in (("enabler",u.OUT,UNI_HASH),("participant",r.OUT,PART_HASH)):
        path=base/"IDENTITY_INVENTORY.json";need(sha(path)==h,"UniProt inventory hash")
        inv=json.loads(path.read_text());need(inv["release"]==u.RELEASE,"UniProt release")
        add("UNIPROT_INVENTORY:"+label,"UniProtKB",u.RELEASE,"identity_inventory",path,"","","CC-BY-4.0","PINNED",h)
        plan=base/"PLAN.json";add("UNIPROT_PLAN:"+label,"CFA",u.RELEASE,"acquisition_plan",plan,"","","internal_metadata","RECORDED")
        for acc,x in sorted(inv["parent_receipts"].items()):
            origin=x.get("origin","enabler_snapshot")
            source=u.OUT if origin=="enabler_snapshot" else r.OUT
            fp=source/"records"/(acc+".json")
            need(sha(fp)==x["receipt_sha256"],"UniProt receipt hash "+acc)
            rec=json.loads(fp.read_text());obj=u.verify(acc,rec)
            need(obj is not None and rec["status"]==200 and rec["sha256"]==x["payload_sha256"],"UniProt payload "+acc)
            if acc in parents:continue
            parents.add(acc)
            add("UNIPROT_PARENT:"+acc,"UniProtKB",u.RELEASE,"identity_response_receipt",fp,rec["url"],rec["acquired_utc"],"CC-BY-4.0","VERIFIED_RESPONSE",x["receipt_sha256"],payload_hash=rec["sha256"],payload_bytes=rec["bytes"])
    need(len(parents)==521,"unique parent count")
    emit("UNIPROT_UNIQUE_PARENT_ROWS",len(parents))
    add("STAGING_DB","CFA","P1A04C_v1","staged_source_records",DB,"","","source_terms_apply","PINNED",DB_HASH)
    for base,version in ((TYPED,"P1A04F_v1"),(STATE,"P1A04H_v1")):
        path=base/"manifest.json";m=json.loads(path.read_text())
        add("DERIVED_MANIFEST:"+version,"CFA",version,"derived_manifest",path,"","","source_terms_apply","VERIFIED")
        for n,x in sorted(m["files"].items()):
            add(version+":"+n,"CFA",version,"derived_table",base/n,"","","source_terms_apply","VERIFIED",x["sha256"])
    head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
    need(not subprocess.check_output(["git","status","--porcelain","--untracked-files=no"],cwd=ROOT,text=True).strip(),"tracked files changed")
    # Bind protocol candidate and registries to this local manifest, without editing receipts.
    docs=["protocol/P1A_SOURCE_REPRESENTATION_POLICY.md","protocol/SOURCE_LICENSE_POLICY.md","protocol/PART1_MATHEMATICAL_PROTOCOL.md","config/part1.json","registries/DECISION_LOG.tsv","registries/AUDIT_REGISTRY.tsv","registries/LIMITATIONS_REGISTRY.tsv"]
    hashes={x:sha(ROOT/x) for x in docs}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".p1a_candidate_",dir=OUT.parent) as td:
        tmp=Path(td)/"payload";tmp.mkdir()
        manifest=tmp/"SOURCE_MANIFEST.tsv"
        with manifest.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(result)
        summary={"status":"REVIEW_REQUIRED_NOT_FROZEN","git_head":head,"source_manifest_sha256":sha(manifest),"manifest_rows":len(result),"core_sources":19,"go_models":1520,"uniprot_unique_parents":521,"deferred":["R02 native semantic verification","R14 full archive","seven unindexed GO-CAM models"],"policy_and_registry_sha256":hashes}
        (tmp/"CANDIDATE.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
        need(not OUT.exists(),"output appeared");tmp.rename(OUT)
    emit("SOURCE_MANIFEST_CANDIDATE",summary)
    emit("CANDIDATE_JSON_SHA256",sha(OUT/"CANDIDATE.json"))
    emit("P1A05_SOURCE_LOCK_CANDIDATE_COMPLETE; REVIEW_REQUIRED_NOT_FROZEN")
if __name__=="__main__":main()
