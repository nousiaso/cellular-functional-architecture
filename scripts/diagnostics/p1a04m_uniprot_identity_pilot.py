#!/usr/bin/env python3
"""Small read-only UniProtKB endpoint and identity schema pilot; no downloads retained."""
import hashlib
import json
import re
import time
from urllib.request import Request,urlopen
from urllib.error import HTTPError
from urllib.parse import urlsplit
from collections import Counter

# Actual unresolved examples from P1A.04L, plus two previously verified controls.
PILOT=("O00160","A0A024B7W1","O00468-6","O00187-PRO_0000027600",
       "P00797-1","P00973","P29728")
BASE="https://rest.uniprot.org/uniprotkb/"
def parent(token):return re.sub(r"-(?:[0-9]+|PRO_[0-9]+)$","",token)
def emit(*a):print(*a,flush=True)
def main():
    releases=set(); results=Counter()
    for token in PILOT:
        acc=parent(token)
        url=BASE+acc+".json"
        req=Request(url,headers={"Accept":"application/json","User-Agent":"CFA-identity-pilot/1.0"})
        try:
            with urlopen(req,timeout=45) as r:
                body=r.read()
                h=r.headers
                final=r.geturl()
                if urlsplit(final).hostname!="rest.uniprot.org":
                    raise RuntimeError("unexpected final host")
                release=h.get("X-UniProt-Release")
                release_date=h.get("X-UniProt-Release-Date")
                if not release or not release_date:raise RuntimeError("release headers absent")
                releases.add((release,release_date))
                obj=json.loads(body)
        except HTTPError as e:
            if e.code not in (404,410):raise
            emit("PILOT_NOT_FOUND",token,acc,e.code)
            results["not_found"]+=1
            continue
        primary=obj.get("primaryAccession")
        secondary=obj.get("secondaryAccessions",[])
        org=obj.get("organism") or {}
        if not isinstance(primary,str) or not isinstance(org.get("taxonId"),int):
            raise RuntimeError("unexpected identity schema")
        isoforms=[]
        for c in obj.get("comments",[]):
            if c.get("commentType")=="ALTERNATIVE PRODUCTS":
                for i in c.get("isoforms",[]):
                    isoforms.extend(i.get("isoformIds",[]))
        features=[f for f in obj.get("features",[]) if f.get("featureId")]
        feature_ids=[f["featureId"] for f in features]
        if token==acc: substatus="parent_token"
        elif "-PRO_" in token:substatus="feature_exact" if token.split("-",1)[1] in feature_ids else "feature_unresolved"
        else:substatus="isoform_exact" if token in isoforms else "isoform_unresolved"
        identity="primary" if acc==primary else "secondary" if acc in secondary else "unresolved_returned_identity"
        emit("PILOT_RECORD",json.dumps({"requested":token,"queried_parent":acc,
             "primary":primary,"identity":identity,"substate":substatus,
             "taxon_id":org["taxonId"],"scientific_name":org.get("scientificName"),
             "entry_type":obj.get("entryType"),"release":release,"release_date":release_date,
             "bytes":len(body),"sha256":hashlib.sha256(body).hexdigest(),
             "top_level_keys":sorted(obj),"isoform_ids":isoforms,
             "feature_count":len(features),"matching_feature":[{"id":f.get("featureId"),"type":f.get("type"),"location":f.get("location")} for f in features if token.endswith("-"+str(f.get("featureId")))],
             "gene_names":[g.get("geneName",{}).get("value") for g in obj.get("genes",[])]},sort_keys=True))
        results[identity]+=1;results[substatus]+=1
        time.sleep(0.25)
    if len(releases)!=1:raise RuntimeError("no unique successful release: "+repr(releases))
    emit("PILOT_RELEASE",sorted(releases))
    emit("PILOT_COUNTS",dict(results))
    emit("P1A04M_UNIPROT_IDENTITY_PILOT_COMPLETE; NO_FULL_ACQUISITION_OR_MERGES")
if __name__=="__main__":main()
