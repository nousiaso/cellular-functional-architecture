#!/usr/bin/env python3
"""Read-only chemical state correspondence and remaining identity inventory."""
from collections import Counter, defaultdict
import gzip
import json
import re
import sqlite3
from pathlib import Path
from p1a04p_chemical_reference_audit import ROOT,STATE,DB,OUT,chebi_id,rows,sha,need,check_manifest,emit
from p1a04k_scope_feasibility_audit import TYPED

CHEBI_HASH="fdc1d2d5c4d3e8227891a2eb955aa4d9433449a8f8830a6b68bc567f770bd216"
UNI_HASH="e9a48281b9b96e481902f09a9346a702aca19355623d81d2d97deecbd6d2b05d"
def main():
    check_manifest(STATE,"P1A04H_v1")
    check_manifest(TYPED,"P1A04F_v1")
    p=OUT/"chebi.obo.gz"; rec=json.loads((OUT/"receipt.json").read_text())
    need(p.stat().st_size==rec["bytes"]==47744176 and sha(p)==rec["sha256"]==CHEBI_HASH and rec["data_version"]=="data-version: 257","ChEBI integrity")
    ip=ROOT/"data/raw/uniprot_2026_04/targeted_go_enablers/IDENTITY_INVENTORY.json"
    need(sha(ip)==UNI_HASH,"UniProt inventory integrity")
    inv=json.loads(ip.read_text()); ut={r["token"]:r for r in inv["token_records"]}
    need(len(ut)==520,"UniProt token count")
    st={r["rdf_id"]:r for r in rows(STATE/"biopax_states.tsv")}
    need(len(st)==61418,"state inventory")
    bp=defaultdict(list); bp_loc=Counter(); bp_tax=defaultdict(set); bp_unknown=set(); bp_other=defaultdict(set)
    def locs(r):
        return {x for _,xs in json.loads(r["location_xrefs_json"]) for ns,x in xs if re.fullmatch(r"GO:[0-9]+",x)}
    for eid,r in st.items():
        refs=json.loads(r["reference_xrefs_json"])
        if r["class"]=="SmallMolecule":
            for _,xs in refs:
                for ns,acc in xs:
                    if ns=="ChEBI":bp[chebi_id(acc)].append(eid)
                    else:bp_other[ns].add(acc)
            bp_loc[tuple(sorted({ns for _,xs in json.loads(r["location_xrefs_json"]) for ns,_ in xs}))]+=1
        elif r["class"]=="Protein":
            for _,xs in refs:
                for ns,acc in xs:
                    if ns in ("UniProt","UniProt Isoform"):
                        bp_tax[acc].update(json.loads(r["taxon_union_json"]))
                        if r["unresolved_taxon"]=="1":bp_unknown.add(acc)
    need(len(bp)==2026,"BP ChEBI count")
    emit("BP_CHEMICAL_LOCATION_NAMESPACE_PROFILES",sorted(bp_loc.items()))
    emit("BP_CHEMICAL_LOCATION_XREF_EXAMPLES",[(eid,st[eid]["location_xrefs_json"]) for eid in sorted({e for es in bp.values() for e in es})[:5]])
    need(DB.stat().st_size==1451589632,"staged DB size")
    db=sqlite3.connect(DB.as_uri()+"?mode=ro",uri=True)
    db.execute("PRAGMA query_only=ON")
    need(dict(db.execute("SELECT key,value FROM metadata")).get("schema_version")=="P1A04C_v1","staged schema")
    # Exact shared identifiers, direct molecule location only. Activity location
    # is NOT imputed to substrates/products, and unmatched location is not proof of incompatibility.
    gm={}; go=defaultdict(list); proteins=Counter(); unknown=set(); other=defaultdict(set)
    for mid,sid,payload in db.execute("SELECT model_id,source_id,payload_json FROM go_record WHERE field='molecules'"):
        obj=json.loads(payload);t=obj["term"];gm[(mid,sid)]=obj
        if t.startswith("CHEBI:"):go[t].append((mid,sid))
        elif t.startswith("UniProtKB:"):
            acc=t.split(":",1)[1]
            if acc in ut:
                tax=ut[acc]["taxon_id"];k="targeted_human" if tax==9606 else "targeted_nonhuman" if tax is not None else "targeted_unresolved"
            elif bp_tax.get(acc) and acc not in bp_unknown:
                tax=bp_tax[acc];k="BP_human" if tax=={"9606"} else "BP_mixed" if "9606" in tax else "BP_nonhuman"
            else:k="not_resolved_by_existing_inventory";unknown.add(acc)
            proteins[k]+=1
        else:other[t.split(":",1)[0]].add(t)
    need(len(gm)==7322 and sum(map(len,go.values()))==3346 and sum(proteins.values())==3864,"GO molecule inventory")
    emit("GO_PROTEIN_PARTICIPANT_TAXONOMY",dict(proteins))
    emit("GO_PROTEIN_PARTICIPANT_ADDITIONAL_UNRESOLVED",len(unknown),sorted(unknown))
    emit("GO_OTHER_MOLECULE_TERMS",[(k,len(v),sorted(v)) for k,v in sorted(other.items())])
    profile=Counter(); examples=defaultdict(list); by_molecule={}
    for t,keys in go.items():
        for mid,sid in keys:
            obj=gm[(mid,sid)];loc=obj.get("located_in") or {};lt=loc.get("term") if isinstance(loc,dict) else None
            es=bp.get(t,[])
            same=[e for e in es if isinstance(lt,str) and lt in locs(st[e])]
            k="no_exact_BP_identifier" if not es else "GO_location_unrecorded" if not lt else "exact_direct_GO_location_match" if same else "no_exact_location_match"
            profile[k]+=1;by_molecule[(mid,sid)]=k
            if len(examples[k])<5:examples[k].append((t,lt,len(es),len(same)))
    emit("GO_CHEMICAL_STATE_CORRESPONDENCE",dict(profile))
    emit("GO_CHEMICAL_CORRESPONDENCE_EXAMPLES",dict(examples))
    roles=Counter()
    for mid,target,pred in db.execute("SELECT model_id,target_id,predicate FROM go_association WHERE field='molecular_associations'"):
        if (mid,target) in by_molecule:roles[(pred,by_molecule[(mid,target)])]+=1
    emit("GO_CHEMICAL_ROLE_CORRESPONDENCE",sorted(roles.items()))
    db.close()
    # Inspect all active ChEBI xref namespaces without guessing equivalence
    # from arbitrary numeric IDs. Report raw property keys for the used terms.
    wanted=set(bp)|set(go)
    xrefns=Counter();xref_examples=defaultdict(set);props=Counter();prop_examples=defaultdict(set);used=0
    ident=None;active=True;xs=[];ps=[]
    def finish():
        nonlocal used
        if not active or not ident:return
        for x in xs:
            ns=x.split(":",1)[0]
            xrefns[ns]+=1
            if len(xref_examples[ns])<3:xref_examples[ns].add(x)
        if ident in wanted:
            used+=1
            keys=set()
            for x in ps:
                key=x.split(None,1)[0];keys.add(key)
                if len(prop_examples[key])<2:prop_examples[key].add(x[:300])
            props.update(keys)
    with gzip.open(p,"rt",encoding="utf-8") as f:
        for line in f:
            s=line.rstrip("\r\n")
            if s.startswith("[") and s.endswith("]"):
                finish();ident=None;active=True;xs=[];ps=[]
            elif s.startswith("id: CHEBI:"):ident=s[4:]
            elif s=="is_obsolete: true":active=False
            elif s.startswith("xref: "):xs.append(s[6:])
            elif s.startswith("property_value: "):ps.append(s[16:])
        finish()
    need(used==len(wanted),"used ChEBI term inventory")
    emit("CHEBI_USED_TERM_PROPERTY_COVERAGE",dict(props))
    emit("CHEBI_USED_TERM_PROPERTY_EXAMPLES",[(k,sorted(v)) for k,v in sorted(prop_examples.items())])
    emit("CHEBI_XREF_NAMESPACE_INVENTORY",[(k,n,sorted(xref_examples[k])) for k,n in xrefns.most_common()])
    emit("P1A04Q_STATE_CORRESPONDENCE_AUDIT_COMPLETE; NO_MERGES_OR_EXCLUSIONS")
if __name__=="__main__":main()
