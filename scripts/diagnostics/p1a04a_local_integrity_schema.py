#!/usr/bin/env python3
"""Read-only P1A.04A local integrity and schema reconnaissance; no extraction/network."""
import sys
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import tarfile
import xml.etree.ElementTree as ET
import zipfile

sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"acquisition"))
import p1a03_acquire_core as core
import p1a03b_acquire_human_gocam as go

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def compact(x, limit=700):
    return json.dumps(x,ensure_ascii=False,sort_keys=True)[:limit]

def top(counter,n=30):
    return counter.most_common(n)

def verify_core():
    with core.CANDIDATES.open() as f:
        rows=[r for r in csv.DictReader(f,delimiter="\t") if r["id"]!="R14"]
    with go.EXPECTED.open() as f:
        expected={r["id"]:r for r in csv.DictReader(f,delimiter="\t")}
    inventory=json.loads((core.RAW/"P1A03_CORE_RECEIPTS.json").read_text())
    require(len(rows)==len(expected)==len(inventory)==19,"core cardinality")
    require(set(expected)=={r["id"] for r in rows},"core ID universe")
    paths={}
    for row, agg in zip(rows,inventory):
        ident=row["id"]
        p=core.destination(row)
        rec=json.loads(p.with_name(p.name+".receipt.json").read_text())
        e=expected[ident]
        require(rec==agg,f"core aggregate mismatch {ident}")
        require(all(rec.get(k)==row[k] for k in row),f"core source identity {ident}")
        require(p.stat().st_size==rec["bytes"]==int(e["bytes"]),f"core size {ident}")
        require(core.sha256(p)==rec["sha256"]==e["sha256"],f"core SHA256 {ident}")
        require(rec["payload_check"]==e["payload_check"],f"core check receipt {ident}")
        require(row["resource"]!="Reactome" or
                rec["reactome_version_before"]==rec["reactome_version_after"]=="98",
                f"core version receipt {ident}")
        paths[ident]=p
    print("CORE_LOCAL_REHASH_PASS",len(paths),sum(p.stat().st_size for p in paths.values()))
    return paths

def shape(value,path,counts,examples):
    counts[path][type(value).__name__]+=1
    if isinstance(value,dict):
        for k,v in value.items():
            shape(v,path+"."+k,counts,examples)
    elif isinstance(value,list):
        for item in value:
            shape(item,path+"[]",counts,examples)
    elif value is not None and path not in examples:
        examples[path]=str(value)[:110]

def verify_models():
    ids,memberships=go.load_indices()
    directory=go.DIRECTORY
    require(not (directory/".acquisition.lock").exists(),"acquisition active or stale lock")
    inventory=directory/"MODEL_INVENTORY.json"
    aggregate=json.loads(inventory.read_text())
    require(len(aggregate)==len(ids)==1520,"model count")
    require([r["model"]["model_id"] for r in aggregate]==ids,"model inventory IDs/order")
    expected={i+".json" for i in ids}|{i+".receipt.json" for i in ids}|{"MODEL_INVENTORY.json"}
    require({p.name for p in directory.iterdir()}==expected,"unexpected or missing model directory entries")
    groups=defaultdict(lambda:defaultdict(Counter))
    examples=defaultdict(dict)
    fieldkeys=defaultdict(Counter)
    additional_shapes=Counter()
    taxon_comparison=Counter()
    empty_molecules=[]
    nonhuman=[]
    samples={}
    for ident,agg in zip(ids,aggregate):
        p=directory/(ident+".json")
        rec=json.loads((directory/(ident+".receipt.json")).read_text())
        require(rec==agg,f"model aggregate mismatch {ident}")
        require(rec["sha256"]==core.sha256(p) and rec["bytes"]==p.stat().st_size,f"model hash/size {ident}")
        require(rec["indexed_taxa"]==sorted(memberships[ident]),f"model taxon receipt {ident}")
        require(rec["url"]==go.BASE+ident+".json" and rec["release_archive"]=="2026-08-05",f"model source {ident}")
        require(rec["model"]==go.inspect(p,ident),f"model summary {ident}")
        obj=json.loads(p.read_text())
        ns=go.namespace(ident)
        for field in ("activities","molecules","objects","provenances"):
            require(isinstance(obj[field],list),f"not list {ident} {field}")
            for item in obj[field]:
                # Count immediate record key signatures; inspect nested structure on samples.
                sig=tuple(sorted(item)) if isinstance(item,dict) else (type(item).__name__,)
                fieldkeys[(ns,field)][sig]+=1
            shape(obj[field],field,groups[ns],examples[ns])
            if obj[field] and (ns,field) not in samples:
                samples[(ns,field)]=(ident,obj[field][0])
        primary=obj.get("taxon")
        additional=obj.get("additional_taxa")
        additional_shapes[type(additional).__name__]+=1
        if additional is None or (isinstance(additional,list) and all(isinstance(x,str) for x in additional)):
            declared=set(additional or [])|({primary} if isinstance(primary,str) else set())
            taxon_comparison["exact_set_match" if declared==memberships[ident] else "set_mismatch"]+=1
            if declared!=memberships[ident] and taxon_comparison["set_mismatch"]<=8:
                print("TAXON_SET_EXAMPLE",ident,sorted(declared),sorted(memberships[ident]))
        else:
            taxon_comparison["uninterpreted_shape"]+=1
            if taxon_comparison["uninterpreted_shape"]<=3:
                print("ADDITIONAL_TAXA_EXAMPLE",ident,compact(additional))
        if primary!=go.HUMAN:
            nonhuman.append(ident)
            print("NONHUMAN_PRIMARY",ident,primary,compact(additional),compact(obj.get("title"),200))
        if not obj["molecules"]:
            empty_molecules.append(ident)
            if len(empty_molecules)<=6:
                print("EMPTY_MOLECULES_EXAMPLE",ident,compact(obj.get("title"),200),compact(obj["activities"][0] if obj["activities"] else None))
    print("MODEL_LOCAL_REHASH_PASS",len(ids),"BYTES",sum(r["bytes"] for r in aggregate))
    print("MODEL_INVENTORY_SHA256",core.sha256(inventory))
    print("MODEL_TAXON_COMPARISON",dict(taxon_comparison),"ADDITIONAL_SHAPES",dict(additional_shapes))
    print("NONHUMAN_PRIMARY_COUNT",len(nonhuman),"EMPTY_MOLECULES_COUNT",len(empty_molecules))
    for ns in sorted(groups):
        for field in ("activities","molecules","objects","provenances"):
            print("GO_RECORD_SIGNATURES",ns,field,compact(top(fieldkeys[(ns,field)],8),1600))
            if (ns,field) in samples:
                ident,item=samples[(ns,field)]
                print("GO_RECORD_EXAMPLE",ns,field,ident,compact(item,1800))
        paths=groups[ns]
        print("GO_SCHEMA_PATH_COUNT",ns,len(paths))
        for path,counts in sorted(paths.items(),key=lambda x:(-sum(x[1].values()),x[0]))[:45]:
            print("GO_SCHEMA_PATH",ns,path,dict(counts),repr(examples[ns].get(path,"")))
    return set(ids)

def local(tag):
    return tag.rsplit("}",1)[-1]

def xml_scan(stream,label):
    # Full streaming syntax/element inventory for selected member only; clear processed elements.
    counts=Counter()
    root=None
    for event,elem in ET.iterparse(stream,events=("start","end")):
        if root is None:
            root=elem.tag
        if event=="end":
            counts[elem.tag]+=1
            elem.clear()
    print("XML_SCAN",label,"ROOT",root,"TOTAL_ELEMENTS",sum(counts.values()),"DISTINCT_TAGS",len(counts))
    print("XML_TAG_COUNTS",label,compact(top(counts,45),4000))

def archive_audit(paths,ids):
    for ident in ("R01","R10","R13"):
        with zipfile.ZipFile(paths[ident]) as z:
            members=[i for i in z.infolist() if not i.is_dir()]
            print("ZIP_INVENTORY",ident,len(members),"UNCOMPRESSED_BYTES",sum(i.file_size for i in members))
            print("ZIP_MEMBER_EXAMPLES",ident,[(i.filename,i.file_size) for i in members[:12]])
            print("ZIP_SUFFIXES",ident,dict(Counter(Path(i.filename).suffix for i in members)))
            if ident=="R01":
                human=[i for i in members if "homo_sapiens" in i.filename.lower() or "homo sapiens" in i.filename.lower()]
                print("BIOPAX_HUMAN_NAME_CANDIDATES",[(i.filename,i.file_size) for i in human[:15]],"COUNT",len(human))
                if len(human)==1 and Path(human[0].filename).suffix.lower() in (".owl",".rdf",".xml"):
                    with z.open(human[0]) as f:
                        xml_scan(f,"R01:"+human[0].filename)
                else:
                    print("BIOPAX_XML_SCAN_DEFERRED_MEMBER_SELECTION")
            if ident=="R13":
                # Filename stem matches are candidates, not proof of equivalent records.
                stems={Path(i.filename).stem for i in members}
                print("REACTOME_GOCAM_FILENAME_STEM_OVERLAP",len(stems & ids),"EXAMPLES",sorted(stems & ids)[:10])
    with tarfile.open(paths["R03"],"r:gz") as t:
        members=[i for i in t.getmembers() if i.isfile()]
        print("SBML_TAR_INVENTORY",len(members),"UNCOMPRESSED_BYTES",sum(i.size for i in members))
        print("SBML_MEMBER_EXAMPLES",[(i.name,i.size) for i in members[:10]])
        candidates=[i for i in members if Path(i.name).suffix.lower() in (".xml",".sbml")]
        print("SBML_XML_CANDIDATES",len(candidates))
        if candidates:
            sample=min(candidates,key=lambda i:(i.size,i.name))
            with t.extractfile(sample) as f:
                xml_scan(f,"R03:"+sample.name)
    with paths["R02"].open("rb") as f:
        print("NATIVE_DUMP_PREFIX_HEX",f.read(48).hex())
    print("NATIVE_DUMP_NOT_OPENED")

def text_audit(paths):
    for ident in ("R04","R05","R06","R07","R08","R09","R11","R12","H01"):
        with paths[ident].open(encoding="utf-8") as f:
            widths=Counter()
            for n,line in enumerate(f,1):
                if n<=3:
                    print("TEXT_SAMPLE",ident,n,repr(line.rstrip("\r\n")[:1500]))
                widths[len(line.rstrip("\r\n").split("\t"))]+=1
        print("TEXT_TAB_WIDTHS",ident,"ROWS",sum(widths.values()),top(widths,12))
    for ident in ("G03","G05","G04"):
        opener=gzip.open if ident=="G04" else open
        with opener(paths[ident],"rt",encoding="utf-8") as f:
            for i,line in enumerate(f):
                if i>=15:
                    break
                print("GO_HEADER",ident,repr(line.rstrip()[:500]))

def main():
    print("CFA_P1A04A_LOCAL_INTEGRITY_SCHEMA_V1",flush=True)
    paths=verify_core()
    ids=verify_models()
    archive_audit(paths,ids)
    text_audit(paths)
    print("P1A04A_RECONNAISSANCE_COMPLETE; SEMANTIC_CERTIFICATION_PENDING")

if __name__=="__main__":
    main()
