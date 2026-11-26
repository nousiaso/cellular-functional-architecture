# P1A.02 source feasibility — preliminary, 2026-11-26

The cluster checkout of 8ae40eafa902ae002fb4741797955ad0535d3e7d was reported synchronized and clean. This validates Git transport, not individual storage quota or the scientific protocol.

Official source documentation reviewed:
- https://reactome.org/download/ and https://reactome.org/download/current/ (v98 filenames)
- https://reactome.org/about/news/299-version-98-released
- https://reactome.org/dev/content-service/ (version endpoint)
- https://zenodo.org/records/24067822 (v98 full archive; 4.7 GB, MD5 01b3820e03ae20e2a55305c5625652ff, CC BY 4.0)
- https://geneontology.org/docs/download-go-cams/ (JSON recommended; archived since June 2026)
- https://release.geneontology.org/2026-08-05/ (versioned candidate; model and taxon indexes, ontology and GAF)
- https://www.genenames.org/download/ (HGNC)

Source selection is provisional. The Reactome current URL is mutable; verify version before AND after acquisition, record timestamps and SHA256, and prefer immutable archives if the version changes or equivalence is uncertain. A version API result alone cannot prove byte-level correspondence. The GO archive date and ontology internal version must be recorded separately. HGNC archive choice, ChEBI and any extra UniProt mapping releases are pending identity coverage. The live directory includes a Reactome GO-CAM export; it is a cross-format control, not independent corroboration. Keep GO hierarchy, GAF, Reactome pathway names and memberships out of primary graph fitting except where explicitly justified for identity/semantics. No functional benchmark or expression data is accessed in this phase.

Run scripts/diagnostics/p1a02_source_preflight.py from a checkout. It makes HEAD requests and bounded GETs for the version and model index/pilot JSONs, writes no data. HEAD errors are not proof of absence (some servers disallow HEAD); resolve before bulk acquisition. JSON pilot is schema reconnaissance, not a representative sample. Source_manifest remains empty until actual acquisition. No source bytes or checksums have been certified. Review output and local quota before download. Do not commit institutional paths or raw logs to public repository.
