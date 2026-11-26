# P1A.03B GO-CAM identifier correction v1.1

Evidence: user cluster read-only diagnostic after the v1 acquisition plan failed before transfers.

- All 1520 human taxon entries are strings and unique; all occur in the HTML directory.
- Complete HTML inventory: 2370 distinct JSON filenames (not 2078).
- Original hexadecimal regex excluded 292 filenames: 133 human-index entries, 2 mouse-index entries, 157 yeast-index entries. Examples include R-HSA-1237112, MGI_MGI_107676 and YeastPathways_ARGSPECAT-PWY.
- 1387 human IDs match the old regex; 133 do not. The examples show Reactome-style IDs. The exact namespace counts for the complete human set will be printed by v1.1.
- Human-indexed fraction: 1520/2370 = 64.14%, not 73.15%. This is directory membership, not a measure of gene or mechanism coverage.
- Reactome-style IDs flag potential source dependence, not proof of identical content; retain records and audit provenance/event crosswalks before graph integration.

Correction: use a conservative safe filename syntax accepting the observed namespaces; fail explicitly for unexpected JSON filenames rather than silently omit. Require the observed 2370 total and 1520 unique matched human IDs from checksum-verified indexes. The pilot now samples up to five IDs per namespace. Internal model ID equality remains strict and has not been verified for Reactome-style records; pilot will expose discrepancies explicitly. No individual models were acquired by the failed plan. Offline fixtures include synthetic Reactome ID structure, not an assertion about real Reactome JSON contents. Do not run the full transfer if plan or pilot fails.

Core acquisition receipts remain unchanged.
