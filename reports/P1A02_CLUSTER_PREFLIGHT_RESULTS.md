# P1A.02 cluster preflight results

Evidence: user-provided terminal output, source script at 6e655a320268396059d975bd55bfc9aac7ed395e.

- Reactome database endpoint returned exactly 98; curl exit 0.
- 20/20 candidate HEAD requests returned HTTP 200, curl exit 0.
- GO-CAM JSON HTML index parser found 2078 distinct model filenames.
- Three deterministic pilots parsed as objects, with common top-level keys: activities, comments, date_modified, id, molecules, objects, provenances, status, taxon, title.
- 5323da1800000002: human NEDD4-model1; production; 2 activities, 1 molecule, 13 objects, 1 provenance.
- 66a3e0bb00001633: human ADM2/CALCRL-RAMP3 signalling; production; 5 activities, 4 molecules, 17 objects, 1 provenance.
- 6a6bcaed00000108: Drosophila NADPH regeneration; production; 3 activities, 6 molecules, 19 objects, 1 provenance.

Accessibility: PASS. Source integrity, license verification, archive equivalence, model coverage and scientific certification: PENDING. The 2078 count is all matched filenames, not a human-model count. Three pilots do not establish full schema or prevalence.

## Next acquisition tranche
Acquire 19 core endpoints (R01-R13, G01-G05, H01); R14 is fallback only. Keep role-separated directories. Mutable snapshots are byte-pinned by local SHA256 and receipts, not claimed to be immutable upstream releases. Reactome version checked before and after each new Reactome transfer and the full run. This is not proof of byte equivalence with the Zenodo archive. No archive extraction, graph building, or individual GO-CAM bulk transfer yet. Taxon index schema diagnostic determines subsequent model list. Source manifest certification remains pending. User must confirm storage allocation/site policy before execute. No raw institutional logs or paths committed.
