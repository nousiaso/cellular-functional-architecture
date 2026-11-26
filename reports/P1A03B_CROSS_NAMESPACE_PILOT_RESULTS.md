# P1A.03B cross-namespace acquisition pilot

Evidence: user-reported cluster execution of p1a03b_acquire_human_gocam.py v1.1 at repository state 7b984ab79129ddbeab6376b1a4e5b253c167b12d (expected from preceding synchronization; exact local HEAD not re-reported).

- Selected and acquired 10/10 models (all new): five hexadecimal IDs and five Reactome_HSA IDs.
- Pilot IDs: 5323da1800000002, 59cc772c00000000, 59dc728000000351, 5e72450500002719, 5e72450500004019, R-HSA-1237112, R-HSA-140180, R-HSA-1475029, R-HSA-1482798, R-HSA-1482801.
- All ten passed strict internal ID validation. Primary model taxonomy NCBITaxon:9606 and status production for all ten.
- One multitaxon-indexed model; one occurrence of additional_taxa field. Aggregate output does not independently establish these are the same record.
- All ten contain list-valued, nonempty activities, molecules, objects and provenances fields. Raw entry totals: activities 81, molecules 104, objects 229, provenances 10. Comments present in nine.
- Total payload bytes 275100. Per-file receipts and checksums retained locally by acquisition script; not independently inspected here.
- Reconciled index unchanged: 2370 HTML IDs, 2363 taxon-index union, 1520 unique human-indexed IDs, 126 human-indexed multitaxon models, 7 unindexed scope-review records.

Pilot acquisition PASS. Authorize full acquisition of the existing 1520 indexed candidates using the unchanged v1.1 script. Reuse ten pilot files after receipt/checksum validation. Do not claim that pilot field completeness generalizes to all models or that Reactome and GO-CAM are independent. Full inventory, source integrity certification, entity-level taxonomy, provenance and seven unindexed records' definitive disposition remain pending.
