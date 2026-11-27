# Certified molecular source inventory

**Part 1A: FROZEN_APPROVED, 2026-11-27.** Based on user-reported cluster outputs; no direct cluster access. Exact approval, manifest and historical policy hashes: [P1A_SOURCE_LOCK_V1.json](../manifests/P1A_SOURCE_LOCK_V1.json). The 2,085-row local manifest and candidate JSON are immutable. No additional downloads or reruns needed. No graph inference or functional validation performed.

## Inventory and key checks
| Component | Result |
|---|---:|
| Core payloads and receipts | 19/19 verified |
| Supplementary reference/inventory hashes | 6/6 verified |
| Derived tables | 13/13 verified |
| Reactome v98 operations / states | 16,187 / 61,418 |
| BP typed assertions | 441,156 |
| Complex component / alternative member links | 40,083 / 36,113 |
| GO archive 2026-08-05 models / activities / molecules | 1,520 / 7,385 / 7,322 |
| GO typed assertions | 15,491 (5,423 causal, 10,068 molecular) |
| Targeted UniProt 2026_04 parent responses / tokens | 521 / 565 |
| Active referenced ChEBI 257 union | 2,233 |
| Manifest records | 2,085 |

SBML comparison: all 2,889 models; 75,233 reaction appearances, 160,266 species appearances, 16,187 unique reactions matched; all 16,179 comparable material input/output sets, 16,187 catalytic/positive/negative sets and 4,400 explicit comparable coefficients agree. Eight template material sets not comparable. 10,145 native controls, unique pairs 6,871 catalyst, 2,162 positive, 1,111 negative. Initial regulatory decoding defect resolved by corrected full rerun. Dependent formats, not independent biological validation.
Scope: 13,624 no listed disease context, 2,052 disease-only, 95 mixed, 416 unmapped. Pathway ancestry unchanged; no exclusions. Raw bridging states 8,098; 345,409 ordered nonself pair occurrences, top two ~60.5%. Cross-source chemical records: 149 exact ID/location, 2,550 missing GO location, 144 nonmatching explicit location, 503 no BP identifier. 682 GO activity occurrences correspond to 669 Reactome events and are dependent. Detailed identity, scope, evidence and molecular composition decisions in approved source policy.

## Limitations and reproducibility
Native R02 dump bytes verified, semantics deferred; R14 full archive not acquired; seven unindexed GO models held separately. Sparse coefficients/locations, ambiguous genes, incomplete complexes, inactive accession and correlated evidence remain explicit. See [source policy](../protocol/P1A_SOURCE_REPRESENTATION_POLICY.md) and [limitations](../registries/LIMITATIONS_REGISTRY.tsv).
Full historical acquisition, diagnostic and audit records remain recoverable in Git history at pre-cleanup commit 2c7a2dfa2fdee19d53b382adc026080bd0809365; not duplicated in the active tree. For historical replay, use a separate checkout of that commit, never overwrite the frozen candidate. Original source certification snapshot is the candidate commit in the lock. Current active inference protocol governs downstream work.
Performance: corrected comparison 3428 s wall, 37.263 s CPU, peak 538.1 MiB; dominant SQLite full scans. Final scope preflight 123 s, CPU 1.636 s; manifest build 111 s, CPU 5.557 s, peak 52.3 MiB. Prefer compact verified tables and profile before parallelizing. All completed jobs exit 0:0.

## Part 1B.01 — source-scoped harmonization (PASS)
User-reported run at bd1bb46bba5b4482660223ba0e008940ac63ce9c; stage SHA256 17a13fb13da7c14fdf296219989c85273782245f0c97eac4a2596e356ffe076d. Source lock and 13 derived tables verified, 1,520 GO models reconstructed, exit 0. Runtime 47.35 s, CPU 10.75 s, peak RSS 341796 KiB. Output manifest SHA256 bb532cdca6d16e6357d06d3e76de6054275c44346c638680a1cfdbba6bcec714 at results/P1B01_SOURCE_SCOPED_HARMONIZATION_V1/manifest.json. Nodes 235,140: BP 125,298 auxiliary, 10,145 controls, 16,187 operations, 61,418 states; GO 7,385 activities, 7,385 enabler observations, 7,322 molecules. Assertions: 441,156 BP, 15,491 GO native, 7,385 GO enabler observations. Direct control unique pairs 6,871 catalyst, 2,162 positive, 1,111 negative. No filters/merges. Template product 8, no typed template input recorded; do not invent one. BP publication materialization still pending. Counts are source-scoped records, not distinct biological mechanisms. Next: read-only compatibility feasibility on actual state/composition and identity records before selecting 1B.02 expansion rules. No graph projection yet.

## Part 1B.02 — feasibility and implementation
User-reported diagnostic: 40,083 AND links (34,916 equal and 5,167 different source location IDs); 36,113 OR links (35,470 equal and 643 different). Feature IDs differ on 18,582 AND and 24,837 OR links. 6,022 OR owners, 354 nested; 40,324 owner-terminal pairs, maximum 1,301. Exact material bridge occurrences 345,409; direct OR child-to-owner 10,790 and owner-to-child 5,392; AND 68,307 and 154,802 respectively. These are counts of possible occurrences, not unique edges; nucleotide-bound complexes dominate large AND interfaces. Factored implementation preserves direct and nested OR witnesses, component prerequisites, identifier/gene candidates and recorded location uncertainty; no graph edges. Execution pending.
