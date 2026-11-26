# P1A.03 core acquisition results (user-reported cluster evidence)

- Commit executed: 374323c1eac4dde1fe9a46570749b0ec13c1999a.
- Offline self-test passed; 19 of 19 core endpoints acquired.
- Total bytes: 992633336 (sum independently recomputed from reported sizes).
- Reactome final version: 98; ZIP/GZIP checks passed on relevant files. Native dump only nonempty, text files only prefix checks.
- Core reported SHA-256 values are recorded in manifests/P1A03_CORE_REPORTED_SHA256.tsv. These are user-reported from terminal output, not independently rehashed off-cluster. Original local receipt JSONs remain authoritative for timestamps/URLs.
- GO taxon index: dictionary with 67 keys; NCBITaxon:9606 maps to a list of 1520 IDs. Directory preflight previously found 2078 unique filenames. Uniqueness, subset membership, multitaxon membership and model-level taxa remain to be audited.

Next: acquire all verified human-indexed IDs without excluding model status or host-pathogen contexts. Check model identity, retain actual taxon/status and counts of raw activity/molecule/object/provenance entries. Those counts are not yet deduplicated molecular operations. Run plan, pilot, then complete collection using three concurrent workers. No graph construction. Other taxa, provenance overlap and detailed model schema are reviewed before primary inclusion. The raw core files and per-model receipts are local only; no unpublished result data is committed. Individual quota and site policy remain user-side checks.
