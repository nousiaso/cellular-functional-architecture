# CFA Part 1 mathematical protocol

Version: 0.1.0 — PROPOSED, not locked. Date: 2026-11-26.

## 1. Estimands and limits
Part 1 estimates reproducible organization of *curated molecular knowledge*, not physiological flux, cell-specific activity, causal effects, or cellular decision rules. A graph cluster is a candidate molecular architecture; a task requires an operational output and independent functional evidence. Missing annotations are unknown, not evidence of absence. Database completeness and curation bias limit all claims.

The core hypothesis of independent convergence is tested only in later parts. In Part 1 we assess whether mechanistic relationships support nontrivial, reproducible multiscale organization beyond membership overlap. Do not condition the primary graph on the three lineages' transcriptomes, later archetypes or experimental screen results.

## 2. Source and identity policy
Proposed primary sources: Reactome and GO-CAM. HGNC and other reference mappings support identity resolution; curated pathway hierarchies and gene-set annotations are held out from primary fitting where possible. Version, URI, hashes, species, license, format, download timestamp and provenance required. Check each real export schema before parser design. Separate human, inferred orthology, host/pathogen and ambiguous events. No silent loss of unknown or unmapped entities. Retain physical state, modification, compartment, isoform, complex composition, entity sets, reaction direction and evidence groups. Do not identify all forms of a protein merely by gene symbol. Distinguish an alternative entity set (OR) from a required complex (AND); unknown logic remains unknown.

## 3. Native graph and candidate projections
Use a typed directed incidence representation with entities E and reactions/activities R. Retain original event IDs, source records and relation types. Let I[e,r] and O[e,r] be binary known input/output indicators. Preserve stoichiometry separately where documented. Unknown stoichiometry is not assumed to be one; binary incidence does not imply stoichiometry. Catalyst, positive/negative regulation, curated precedence and shared machinery are separate layers. Explicit reversible events require documented direction handling. GO-CAM activities and Reactome events are not assumed identical; cross-source mapping must be auditable and ambiguous links retained separately.

A candidate material-dependency projection is A_flow = O^T D I, with D a diagonal entity-weight matrix. Remove self loops for this projection only; retain native records. An edge denotes a potential output-to-input dependency, not observed flow or reaction feasibility. Match molecular state and compartment; explicit transport events bridge compartments. Keep signed causal regulation in separate positive and negative layers. Do not symmetrize or collapse relation types without explicit justification.

Compare binary unweighted incidence, specificity weighting w_e = log((N+1)/(d_e+1)) with N the number of reactions and d_e the number involving e, and a predeclared currency-entity policy. Native graph retains all entities. Audit high-degree entities, degree distributions, edge multiplicity, isolated events, components, source coverage, and loss from projection. A projection is a derived, non-independent set of edges.

## 4. Model candidates (not selected)
Baseline: curated pathway gene overlap, size- and degree-aware similarity; pathway-only baseline must be explicitly distinguished from reaction-level shared-machinery baseline. Compare directed Infomap and Leiden (only on compatible positive layers, with direction/weight transformations recorded), and degree-corrected hierarchical stochastic block models. A nested hard partition is not an overlapping ontology: compare link communities or mixed membership if supported. Model directed, signed and weighted layers with appropriate likelihoods; a binary SBM is only a simplified example. Avoid arbitrary weighted layer sums. Select by blocked validation, stability and complexity; do not use later functional outcomes for selection.

For a probabilistic candidate, description length Sigma = -log P(data | partition, parameters) + model coding cost. Only compare description lengths on the same observations with compatible coding conventions; use held-out predictive scores for comparisons across representations. Evaluate whether modules represent assortative communities, complementary roles, or generic network hubs. A hierarchy may require directed relationships or overlapping composition rather than a single tree.

## 5. Hypothesis tests and nulls
H1: actual typed mechanistic structure is nontrivial relative to constrained incidence nulls.
H2: mechanistic relationships provide incremental information beyond shared genes/machinery and degrees.
H3: candidate architectures are robust to evidence-group resampling and reasonable parameter changes.
H4: multiple levels of organization are supported beyond a flat model, measured on identical held-out observations.
H5: candidate architectures contain complementary operations beyond near-duplicate annotations, with directed dependency evidence and bounded redundancy.

Define the specific statistics and admissible randomizations only after schema/coverage audit and before fitting final models. Randomize native incidence while preserving role-specific degrees, species, compartments and other feasible structural constraints; separate nulls for causal edges and shared machinery. Preserve or block equivalent events and evidence groups. Check swap mixing, null diversity and constraint satisfaction. Degree matching alone cannot guarantee biological plausibility. For prediction, hold out native evidence/near-duplicate event groups BEFORE projection; remove all descendants of held-out evidence and account for cross-source duplication. Unobserved edges are not confirmed negatives; interpret prediction as recovery of curated evidence, not biological truth. Compare on identical held-out targets, candidate universes and split groups; record uncertainty. No random projected-edge split as primary validation.

Empirical p = (1 + number of null statistics >= observed)/(B+1). Repeat model selection and module matching inside null runs if used in the test. Set test families, effect sizes, permutation budget and FDR method before confirmatory analysis. Use dependence-aware inference or conservative correction as warranted. No fixed sample size or q threshold is certified yet. Stability assessed via aligned memberships (Jaccard, split/merge-aware summaries), evidence-block resampling, and seed sensitivity. Report low coverage and non-identifiability rather than force a partition.

## 6. Functional specification and freeze
For each module retain stable ID, member events, molecular roles, evidence, input/output candidates, complementary suboperations, alternatives, uncertainty, overlap, and stability. Separate mechanistically connected modules from operationally interpretable candidate tasks; do not force functional names. Keep catalogue, graph build, mappings and protocol checksums. Preserve pre-existing broad terms as annotation/benchmark rather than proof. No empirical task validation claim from Part 1.

## 7. Acceptance gates
G0 source integrity and provenance; G1 identity and extraction fidelity (native-vs-export checks, synthetic fixtures, missingness); G2 graph feasibility and null feasibility; G3 predeclared model comparison and robustness; G4 module and hierarchy certification with explicit scope and limits; G5 frozen catalogue and protocol before opening Part 2. Gate failures trigger a documented review, not automatic threshold relaxation. No numerical cutoffs are fixed until feasibility evidence is available. Data collection and exploratory diagnostics are not confirmatory tests.

## 8. Methodological references for review
- Peixoto, 2014, hierarchical block structures: https://doi.org/10.1103/PhysRevX.4.011047
- Peixoto, 2015, layered networks: https://doi.org/10.1103/PhysRevE.92.042807
- Ahn et al., 2010, link communities: https://doi.org/10.1038/nature09182
- PAC, 2022, pathway communities: https://doi.org/10.1093/nargab/lqac044
- Hart et al., 2015, Pareto task inference (later phase): https://doi.org/10.1038/nmeth.3254
