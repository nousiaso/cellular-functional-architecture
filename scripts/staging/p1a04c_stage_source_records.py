#!/usr/bin/env python3
"""P1A.04C source-faithful staging. No semantic merging or graph inference."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import xml.etree.ElementTree as ET
import zipfile

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/"data/raw"
GO=RAW/"go_2026-08-05/human_indexed_models"
BIO=RAW/"reactome_v98/primary/biopax.zip"
OUT=ROOT/"results/P1A04C_SOURCE_STAGING_V1/source_records.sqlite"
RDF="{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
XML="{http://www.w3.org/XML/1998/namespace}"
BP="{http://www.biopax.org/release/biopax-level3.owl#}"
BIO_HASH="5ff50502d2a35e42887dae29a1636491faaeec7bea23d4b7371584f1eb7885f3"
INV_HASH="03418cc71f0195d49f89f3bcbad09ca1e6416867c58afa0d0fe48304a4034d09"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def j(x):
    return json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(",",":"))

def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)

DDL="""
CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
CREATE TABLE bp_record(
 ordinal INTEGER PRIMARY KEY, rdf_id TEXT UNIQUE, class_uri TEXT NOT NULL,
 xml_payload TEXT NOT NULL);
CREATE TABLE bp_property(
 record_ordinal INTEGER NOT NULL REFERENCES bp_record(ordinal),
 property_ordinal INTEGER NOT NULL, predicate_uri TEXT NOT NULL,
 target_raw TEXT, target_local_id TEXT, literal_text TEXT, attributes_json TEXT NOT NULL,
 PRIMARY KEY(record_ordinal,property_ordinal));
CREATE TABLE go_model(
 model_id TEXT PRIMARY KEY, payload_json TEXT NOT NULL, sha256 TEXT NOT NULL,
 primary_taxon TEXT, indexed_taxa_json TEXT NOT NULL);
CREATE TABLE go_record(
 model_id TEXT NOT NULL REFERENCES go_model(model_id), field TEXT NOT NULL,
 ordinal INTEGER NOT NULL, source_id TEXT, payload_json TEXT NOT NULL,
 PRIMARY KEY(model_id,field,ordinal));
CREATE TABLE go_association(
 model_id TEXT NOT NULL REFERENCES go_model(model_id),
 activity_ordinal INTEGER NOT NULL, activity_id TEXT NOT NULL,
 field TEXT NOT NULL, ordinal INTEGER NOT NULL, predicate TEXT,
 target_id TEXT NOT NULL, payload_json TEXT NOT NULL,
 PRIMARY KEY(model_id,activity_ordinal,field,ordinal));
"""
def stage_bp(db,stream):
    stack=[]
    base=None
    classes=Counter()
    n=0
    def norm(s):
        if s is None:
            return None
        if s.startswith("#"):
            return s[1:]
        if base and s.startswith(base.rstrip("#")+"#"):
            return s[len(base.rstrip("#"))+1:]
        return s
    for event,e in ET.iterparse(stream,events=("start","end")):
        if event=="start":
            stack.append(e)
            if len(stack)==1:
                need(e.tag==RDF+"RDF","unexpected BioPAX root")
                base=e.get(XML+"base")
                db.execute("INSERT INTO metadata VALUES (?,?)",("biopax_xml_base",j(base)))
            continue
        if len(stack)==2:
            n+=1
            ident=norm(e.get(RDF+"ID") or e.get(RDF+"about"))
            db.execute("INSERT INTO bp_record VALUES (?,?,?,?)",
                       (n,ident,e.tag,ET.tostring(e,encoding="unicode")))
            rows=[]
            for k,c in enumerate(e):
                target=c.get(RDF+"resource")
                rows.append((n,k,c.tag,target,norm(target),c.text,j(c.attrib)))
            db.executemany("INSERT INTO bp_property VALUES (?,?,?,?,?,?,?)",rows)
            classes[e.tag]+=1
            stack[-2].remove(e)
            e.clear()
            if n%100000==0:
                print("BIOPAX_STAGED",n,flush=True)
        stack.pop()
    return classes

def stage_go(db,mid,obj,digest,taxa):
    need(obj.get("id")=="gomodel:"+mid,"GO model identity "+mid)
    db.execute("INSERT INTO go_model VALUES (?,?,?,?,?)",
               (mid,j(obj),digest,obj.get("taxon"),j(taxa)))
    ids={}
    for field in ("activities","molecules","objects","provenances"):
        records=obj[field]
        need(isinstance(records,list),"GO field "+mid+" "+field)
        seen=set()
        for k,r in enumerate(records):
            need(isinstance(r,dict),"GO record not object "+mid)
            ident=r.get("id")
            if field in ("activities","molecules"):
                need(isinstance(ident,str) and ident not in seen,
                     "missing or duplicate local GO ID "+mid)
                seen.add(ident)
            db.execute("INSERT INTO go_record VALUES (?,?,?,?,?)",(mid,field,k,ident,j(r)))
        ids[field]=seen
    for k,a in enumerate(obj["activities"]):
        for field,target_key,target_field in (
            ("causal_associations","downstream_activity","activities"),
            ("molecular_associations","molecule","molecules")):
            rels=a.get(field,[])
            need(isinstance(rels,list),"non-list GO associations "+mid)
            for i,r in enumerate(rels):
                target=r.get(target_key)
                need(isinstance(target,str) and target in ids[target_field],
                     "unresolved local GO target "+mid+" "+str(target))
                db.execute("INSERT INTO go_association VALUES (?,?,?,?,?,?,?,?)",
                           (mid,k,a["id"],field,i,r.get("predicate"),target,j(r)))

def scalar(db,sql,args=()):
    return db.execute(sql,args).fetchone()[0]

def validate(db,full):
    db.executescript("""
    CREATE INDEX bp_property_pred ON bp_property(predicate_uri);
    CREATE INDEX bp_property_target ON bp_property(target_local_id);
    CREATE INDEX go_record_source ON go_record(model_id,field,source_id);
    CREATE INDEX go_association_pred ON go_association(predicate);
    """)
    unresolved=db.execute("""
      SELECT p.predicate_uri,p.target_raw,COUNT(*)
      FROM bp_property p LEFT JOIN bp_record r ON p.target_local_id=r.rdf_id
      WHERE p.target_raw IS NOT NULL AND r.ordinal IS NULL
      GROUP BY p.predicate_uri,p.target_raw
    """).fetchall()
    print("BIOPAX_UNRESOLVED_REFERENCES",unresolved)
    if full:
        need(unresolved==[("{http://www.w3.org/2002/07/owl#}imports",
                           "http://www.biopax.org/release/biopax-level3.owl",1)],
             "unexpected BioPAX unresolved refs")
        need(scalar(db,"SELECT COUNT(rdf_id) FROM bp_record")==548688,"BioPAX ID count")
        for cls,n in {"BiochemicalReaction":16159,"Catalysis":6871,"Control":3271,
                      "Protein":33782,"Complex":17654,"Degradation":20,
                      "TemplateReaction":8,"TemplateReactionRegulation":3}.items():
            actual=scalar(db,"SELECT COUNT(*) FROM bp_record WHERE class_uri=?",(BP+cls,))
            need(actual==n,"BioPAX class count "+cls)
        need(scalar(db,"SELECT COUNT(*) FROM go_model")==1520,"GO model count")
        for f,n in {"activities":7385,"molecules":7322,"objects":27304,"provenances":1520}.items():
            need(scalar(db,"SELECT COUNT(*) FROM go_record WHERE field=?",(f,))==n,"GO count "+f)
        for f,n in {"causal_associations":5423,"molecular_associations":10068}.items():
            need(scalar(db,"SELECT COUNT(*) FROM go_association WHERE field=?",(f,))==n,"GO association count "+f)
    for table in ("bp_record","bp_property","go_model","go_record","go_association"):
        print("STAGED_TABLE",table,scalar(db,"SELECT COUNT(*) FROM "+table))
    print("GO_ASSOCIATION_COUNTS",db.execute("SELECT field,predicate,COUNT(*) FROM go_association GROUP BY field,predicate ORDER BY field,predicate").fetchall())
    need(db.execute("PRAGMA foreign_key_check").fetchone() is None,"foreign key failure")
    need(scalar(db,"PRAGMA integrity_check")=="ok","SQLite integrity failure")

def run():
    need(sha(BIO)==BIO_HASH,"BioPAX checksum changed")
    inv=GO/"MODEL_INVENTORY.json"
    need(sha(inv)==INV_HASH,"GO inventory checksum changed")
    receipts=json.loads(inv.read_text())
    need(len(receipts)==1520,"inventory cardinality")
    OUT.parent.mkdir(parents=True,exist_ok=True)
    need(not OUT.exists(),"output exists; do not overwrite")
    lock=OUT.parent/".staging.lock"
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    os.close(fd)
    try:
        with tempfile.TemporaryDirectory(prefix=".staging-",dir=OUT.parent) as tmp:
            p=Path(tmp)/"source_records.sqlite"
            db=sqlite3.connect(p)
            try:
                db.execute("PRAGMA foreign_keys=ON")
                db.executescript(DDL)
                db.executemany("INSERT INTO metadata VALUES (?,?)",[
                    ("schema_version","P1A04C_v1"),
                    ("biopax_zip_sha256",BIO_HASH),
                    ("go_inventory_sha256",INV_HASH),
                    ("scope","source records only; not certified canonical graph")])
                with zipfile.ZipFile(BIO) as z,z.open("Homo_sapiens.owl") as f:
                    counts=stage_bp(db,f)
                print("BIOPAX_CLASS_COUNTS",dict(counts),flush=True)
                for r in receipts:
                    mid=r["model"]["model_id"]
                    source=GO/(mid+".json")
                    need(source.stat().st_size==r["bytes"] and sha(source)==r["sha256"],
                         "GO checksum changed "+mid)
                    stage_go(db,mid,json.loads(source.read_text()),r["sha256"],r["indexed_taxa"])
                validate(db,True)
                db.commit()
            finally:
                db.close()
            # Publish complete database only after validation. No overwrite.
            need(not OUT.exists(),"output appeared during staging")
            os.link(p,OUT)
        print("OUTPUT",OUT.relative_to(ROOT),"BYTES",OUT.stat().st_size,"SHA256",sha(OUT))
        print("P1A04C_STAGING_PASS; CANONICALIZATION_AND_SOURCE_CERTIFICATION_PENDING")
    finally:
        lock.unlink()

def self_test():
    import io
    db=sqlite3.connect(":memory:")
    db.execute("PRAGMA foreign_keys=ON")
    db.executescript(DDL)
    xml=b'<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#" xml:base="http://example.org/#"><bp:Protein rdf:ID="P1"/><bp:BiochemicalReaction rdf:ID="R1"><bp:left rdf:resource="#P1"/></bp:BiochemicalReaction></rdf:RDF>'
    stage_bp(db,io.BytesIO(xml))
    obj={"id":"gomodel:x","taxon":"NCBITaxon:9606","activities":[
        {"id":"a","causal_associations":[{"predicate":"RO:example","downstream_activity":"b"}]},
        {"id":"b","molecular_associations":[{"predicate":"RO:example","molecule":"m"}]}],
        "molecules":[{"id":"m","term":"UniProtKB:example"}],"objects":[],"provenances":[]}
    stage_go(db,"x",obj,"fixture",["NCBITaxon:9606"])
    need(scalar(db,"SELECT COUNT(*) FROM bp_record")==2,"fixture BP")
    need(scalar(db,"SELECT COUNT(*) FROM go_association")==2,"fixture GO")
    validate(db,False)
    db.close()
    print("SELF_TEST_PASS (synthetic source-record and reference fixtures)")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--execute",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
    elif args.execute:
        run()
    else:
        print("PLAN_ONLY",OUT.relative_to(ROOT))
        print("Reads checksum-pinned BioPAX and GO-CAM inputs; creates new local SQLite DB only with --execute")
