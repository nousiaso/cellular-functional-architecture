# P1A.03B corrected index plan results

User-reported cluster plan at 633a71f0a0a581371796d093027b40323c7e6bb6: successful, no transfers in plan mode.

- HTML 2370 distinct IDs, taxon union 2363, 67 taxon keys.
- Human index 1520 unique IDs, zero duplicates, zero absent from HTML.
- 7 HTML IDs absent from every taxon index; zero indexed IDs absent from HTML.
- All HTML namespaces: hexadecimal 2078, YeastPathways 157, Reactome_HSA 133, MGI 2.
- Human namespaces: hexadecimal 1387, Reactome_HSA 133.
- 126 human-indexed models also indexed under another taxon; 1394 human-only indexed. These are index memberships, not molecular organism assignments.
- Plan-only selected 1520.

Status: index plan PASS; pilot/full transfer pending. Audit the seven unindexed models using bounded read-only JSON requests; report their actual internal IDs, taxon and status. Do not silently classify them as nonhuman or add them to the pinned index universe without review. Retain the 126 multitaxon models for subsequent entity-level taxonomy review. Reactome-style IDs indicate potential source dependence, not exact equivalence or exhaustive provenance classification. No model identity assumptions beyond the tested hexadecimal pilots are certified.
