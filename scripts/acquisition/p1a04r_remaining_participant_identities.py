#!/usr/bin/env python3
"""Consolidated remaining participant identity and chemical xref candidate audit."""
import argparse
from collections import Counter,defaultdict
import gzip
import json
from pathlib import Path
import re
import sqlite3
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"diagnostics"))
import p1a04n_acquire_targeted_uniprot as u
from p1a04q_state_correspondence_audit import CHEBI_HASH,UNI_HASH
from p1a04p_chemical_reference_audit import DB,OUT as CHEBI_OUT,ROOT,STATE,rows,sha,need,check_manifest
OUT=ROOT/"data/raw/uniprot_2026_04/targeted_go_participants"
EXPECTED=set("B4DWJ1 O00187-PRO_0000027598 O15455-1 O43353-1 O60313-PRO_0000007397 O95467 O95825 O96017-10 P00751-PRO_0000027545 P01024-PRO_0000005909 P01031-PRO_0000005987 P02458-PRO_0000005731 P06681-PRO_0000027610 P09871-PRO_0000027586 P0C0L4-PRO_0000005968 P0CAP2 P16298-2 P17181-2 P19838-2 P31153-2 P40763-1 P41161 P48740-PRO_0000027592 P55017-3 P61960 P84996 P94416 Q12933-1 Q13323 Q13557-11 Q15437 Q15555 Q5GH70 Q5GH76 Q66K74 Q6SPF0 Q7L2J0 Q8NB37 Q8TB45-2 Q8TEA7 Q96QA5 Q9H6D3 Q9NNX6-12 Q9NQ89 Q9NX02".split())
def emit(*a): print(*a,flush=True)
def inventory():
    check_manifest(STATE,"P1A04H_v1")
    ip=u.OUT/"IDENTITY_INVENTORY.json"
    need(sha(ip)==UNI_HASH,"enabler inventory integrity")
    inv=json.loads(ip.read_text())
    ut={r["token"]:r for r in inv["token_records"]}
    bt=defaultdict(set);bu=set();chem=defaultdict(set)
    for r in rows(STATE/"biopax_states.tsv"):
        for _,xs in json.loads(r["reference_xrefs_json"]):
            for ns,acc in xs:
                if r["class"]=="Protein" and ns in ("UniProt","UniProt Isoform"):
                    bt[acc].update(json.loads(r["taxon_union_json"]))
                    if r["unresolved_taxon"]=="1":bu.add(acc)
                if r["class"]=="SmallMolecule" and ns!="ChEBI":
                    chem[ns].add(acc)
    need(DB.stat().st_size==1451589632,"staged DB size")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","staged schema")
    gaps=Counter()
    for payload, in db.execute("SELECT payload_json FROM go_record WHERE field='molecules'"):
        t=json.loads(payload)["term"]
        if not t.startswith("UniProtKB:"):continue
        x=t.split(":",1)[1]
        if x not in ut and (not bt.get(x) or x in bu):gaps[x]+=1
    db.close()
    need(set(gaps)==EXPECTED and sum(gaps.values())==91,"participant gap inventory changed")
    return inv,gaps,chem
def chemical_audit(chem):
    p=CHEBI_OUT/"chebi.obo.gz"
    need(sha(p)==CHEBI_HASH,"ChEBI hash")
    # Only observed, unambiguous namespace correspondences are considered.
    mapping={"pubchem.compound":"PubChem Compound","chemspider":"ChemSpider"}
    candidates=defaultdict(set); ident=None; active=True; xs=[]
    def finish():
        if not ident or not active:return
        for x in xs:
            # Actual release 257 xref syntax from Q; retain only exact numeric tokens.
            m=re.match(r"^(pubchem\.compound|chemspider):([0-9]+)(?:\s|$)",x)
            if m:
                ns,acc=mapping[m.group(1)],m.group(2)
                if acc in chem[ns]:candidates[(ns,acc)].add(ident)
    with gzip.open(p,"rt",encoding="utf-8") as f:
        for line in f:
            s=line.rstrip("\r\n")
            if s.startswith("[") and s.endswith("]"):
                finish();ident=None;active=True;xs=[]
            elif s.startswith("id: CHEBI:"):ident=s[4:]
            elif s=="is_obsolete: true":active=False
            elif s.startswith("xref: "):xs.append(s[6:])
        finish()
    for ns,ids in sorted(chem.items()):
        if ns not in mapping.values():
            emit("CHEMICAL_NAMESPACE_NOT_MAPPED",ns,len(ids))
            continue
        counts=Counter("unique_candidate" if len(candidates[(ns,x)])==1 else "multiple_candidates" if candidates[(ns,x)] else "no_candidate" for x in ids)
        emit("CHEMICAL_XREF_CANDIDATE_COUNTS",ns,dict(counts))
        emit("CHEMICAL_XREF_CANDIDATES",ns,[(x,sorted(candidates[(ns,x)])) for x in sorted(ids) if candidates[(ns,x)]])
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--execute",action="store_true");a=ap.parse_args()
    inv,gaps,chem=inventory();ts=sorted(gaps);parents=sorted({u.parent(x) for x in ts})
    need(all(u.TOKEN.fullmatch(x) for x in ts),"token syntax")
    old=inv["parent_receipts"];existing=[x for x in parents if x in old]
    emit("PARTICIPANT_ACQUISITION_PLAN","tokens",len(ts),"occurrences",sum(gaps.values()),"parents",len(parents),"existing_parent_receipts",len(existing),"new_parent_requests",len(parents)-len(existing),"release",u.RELEASE)
    chemical_audit(chem)
    if not a.execute:emit("PLAN_ONLY; NO_DOWNLOAD");return
    OUT.mkdir(parents=True,exist_ok=True)
    plan={"schema_version":"P1A04R_v1","release":u.RELEASE,"source_inventory_sha256":UNI_HASH,"tokens":ts,"parents":parents}
    pp=OUT/"PLAN.json"
    if pp.exists():need(json.loads(pp.read_text())==plan,"plan drift")
    else:u.atomic(pp,plan)
    rd=OUT/"records";rd.mkdir(exist_ok=True)
    need(not list(OUT.rglob("*.partial")),"partial files require review")
    need(set(x.name for x in rd.iterdir())<={x+".json" for x in parents if x not in old},"unexpected receipt")
    objs={};receipts={};counts=Counter()
    for acc in parents:
        if acc in old:
            p=u.OUT/"records"/(acc+".json")
            need(sha(p)==old[acc]["receipt_sha256"],"original receipt integrity "+acc)
            rec=json.loads(p.read_text());obj=u.verify(acc,rec);origin="enabler_snapshot"
        else:
            p=rd/(acc+".json")
            if p.exists():rec=json.loads(p.read_text());obj=u.verify(acc,rec);origin="participant_snapshot"
            else:
                try:rec=u.fetch(acc);obj=u.verify(acc,rec)
                except Exception:
                    emit("PARTICIPANT_ACQUISITION_FAILED",acc);raise
                u.atomic(p,rec);origin="participant_snapshot";time.sleep(.25)
        objs[acc]=obj
        receipts[acc]={"origin":origin,"receipt_sha256":sha(p),"payload_sha256":rec["sha256"],"status":rec["status"]}
        counts[origin]+=1
    result=[]
    for token in ts:
        obj=objs[u.parent(token)]
        r=u.identity(token,obj) if obj is not None else {"token":token,"queried_parent":u.parent(token),"identity":"lookup_unresolved","substate":"unresolved","taxon_id":None}
        result.append(r)
    out={"schema_version":"P1A04R_v1","release":u.RELEASE,"token_records":result,"parent_receipts":receipts}
    rp=OUT/"IDENTITY_INVENTORY.json"
    if rp.exists():need(json.loads(rp.read_text())==out,"inventory drift")
    else:u.atomic(rp,out)
    emit("PARTICIPANT_RECEIPT_ORIGINS",dict(counts))
    emit("PARTICIPANT_IDENTITY_STATUS",sorted(Counter((r["identity"],r["substate"]) for r in result).items()))
    def tax(r):
        return "human" if r["taxon_id"]==9606 else "unresolved" if r["taxon_id"] is None else "nonhuman"
    emit("PARTICIPANT_TAXONOMY_TOKENS",dict(Counter(tax(r) for r in result)))
    occ=Counter()
    for r in result:occ[tax(r)]+=gaps[r["token"]]
    emit("PARTICIPANT_TAXONOMY_OCCURRENCES",dict(occ))
    emit("PARTICIPANT_UNRESOLVED_EXACT",[r for r in result if r["identity"] not in ("primary","secondary") or r["substate"] not in ("parent_token","isoform_exact","feature_exact")])
    emit("PARTICIPANT_IDENTITY_INVENTORY_SHA256",sha(rp))
    emit("P1A04R_PARTICIPANT_AND_CHEMICAL_CANDIDATE_AUDIT_COMPLETE; NO_MERGES")
if __name__=="__main__":main()
