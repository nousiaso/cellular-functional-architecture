# Part 1 mathematical protocol — v0.3.0

Approved inference direction, 2026-11-27. Confirmatory implementation not frozen. The [source lock](../manifests/P1A_SOURCE_LOCK_V1.json) and approved source representation remain unchanged.

## 1. Goal and estimand
Infer organization of curated molecular knowledge beyond known pathways, then identify complementary mechanistic assemblies. Distinguish structural role modules, composite architectures and experimentally supported functions. Label-excluded discovery is not function-blind: topology reflects curation. Missing assertions are unknown, not biological negatives. No flux, cell activity or causal-function claims from Part 1.
Reactome BioPAX is primary discovery. Preserve GO-CAM separately for eligible nonimported cross-source corroboration; GO-only and pooled views are secondary. Audit shared publications/imports and unknown provenance. Exclude curated precedence, pathway memberships/hierarchy and gene sets from primary fitting; use them for benchmarks, pathway-conditioned tests and separate sensitivity analyses. No expression or screen data before catalogue freeze.

## 2. Molecular representation
H=(E,R,L,P): states E, operations/activities R, typed assertions L, provenance P. Preserve source-local IDs, compartment, modifications, molecular reference, required components (AND), alternative members (OR), direction, explicit coefficients and uncertainty. Binary incidence is participation, not stoichiometry. Do not merge by gene, chemical ID, fingerprint or citation. Human feasibility H/F/U follows approved source policy; distinguish foreign machinery from human host-response requirements. Disease annotation alone is not exclusion; compare explicit scope strata. Preserve many-to-many gene links and complex logic.

Let I[e,r], O[e,r], C[e,r], Pplus[e,r], Pminus[e,r] denote recorded material input/output, catalyst and positive/negative regulation. Separate template/product, composition, native causal predicates and context.

Candidate exact-state projections (nonnegative diagonal weights D):
```text
Aflow  = O^T D I
Acat   = O^T Dcat C
Aplus  = O^T Dplus Pplus
Aminus = O^T Dminus Pminus
Ashared = C^T Dshared C   (shared machinery baseline)
```
Remove self pairs only in projected view. These are possible dependencies, not flux or measured regulation. Compare unweighted, inverse frequency w_e=log((N+1)/(d_e+1)), bridge normalization w_e=1/b_e for b_e>0, and prespecified currency-molecule sensitivity. Here d_e counts distinct incident operations under a registered definition; b_e=nout(e)*nin(e)-sum_r O[e,r]I[e,r]. Native molecules retained. Measure molecular bridge frequency, cross-pathway contribution, degree and class; separate unmapped operations. No biologically hand-picked weights.

Expanded connectivity:
- OR: K_or[e,f]=1 if compatible recorded alternative realizations intersect; distinguish possible overlap from guaranteed inclusion R(e) subset R(f) with complete recorded realizations. A_or=O^T K_or I. Require state/location compatibility and realization witnesses; keep exact/added-only separate.
- AND: K_and[e,f]=1 if e is a recorded required component of f under compatible context. A_component=O^T K_and I expresses a prerequisite, not production of the whole complex. Retain other requirements and explicit assembly events; no invented transformations.
- Context: U[q,r] records context; U^T Dq U is compatibility, not causal dependency. Do not infer pH, timing or transport from location. Explicitly assess domain coverage for metabolism, signalling, trafficking, ECM, barrier and host-response mechanisms.
- Cross-source: strict exact identifier and compatible explicit context; missing-location candidates separate. No automatic equivalence. GO mappings not used in primary discovery.
Keep layers separate. Do not sum arbitrary weights, treat inhibition as negative random-walk probability or count projected links as independent observations.

## 3. Two-level inference
**Level 1: structural role modules.** Prioritize native typed heterogeneous/bipartite incidence and degree-corrected directed block models. Compare projected models, directed Infomap, compatible Leiden objectives and overlap alternatives. Illustrative count model:
```text
Aij^k ~ Poisson(theta_out[i,k]*theta_in[j,k]*Omega^k[z_i,z_j])
Sigma = -log P(data|model) + coding_cost
```
z_i is a structural block, not a functional label. Normalize degree parameters for identifiability; define admissible pairs and appropriate binary/count/weighted likelihoods. Independence of derived edges/layers is not established. Compare description length only on identical observations and compatible codes; use common held-out targets across representations. No fixed hierarchy, number of blocks or algorithm.

**Level 2: composite assemblies.** For hard memberships Z[i,a]=1(z_i=a), W^k=Z^T A^k Z summarizes block interfaces. A block path A->B->C does not prove a continuous mechanism: different operations inside B may be disconnected. Define F=(V,L,I_boundary,O_boundary,Gamma), with operation set V, supported typed links L, candidate boundaries and explicit molecular/contextual requirements Gamma. Search bounded directed paths, branches, feedback and alternatives with operation-level interface witnesses; allow overlapping assemblies and reused blocks. Path witnesses do not prove simultaneous realizability; unknown requirements stay unknown.
Evaluate mechanistic support, complementarity, boundary coherence, redundancy and complexity separately. Candidate constrained objective: maximize mechanistic support subject to bounds on size/redundancy, minimum complementarity and interface requirements. Define search space, thresholds and uncertainty treatment before confirmatory testing; no arbitrary weighted score or attractive biological name selection.
Candidate residual W_ab-E_null(W_ab|degree,pathway,context), standardized only for nonzero null SD. Test on holdout or account for selection under the null.

## 4. Primary tests and curation controls
H02: incremental predictive information beyond curated pathways, degree, shared machinery and gene overlap. H05: operation-witnessed complementary cross-pathway assemblies beyond pathway-conditioned expectation and annotation redundancy. H03 robustness required; H01 basic organization diagnostic; H04 hierarchy secondary; H06 temporal extension.
For overlapping pathway membership X[i,p], Cij=sum_p X[i,p]X[j,p] at specified levels. Compare hard partitions using AMI/NMI only on well-defined hard reference partitions: NMI=2I(Z;C)/(H(Z)+H(C)); use overlap-aware co-membership for overlapping pathways and assemblies. Report size, level, mapping coverage and matched universe. Neither high nor low agreement alone establishes novelty.
Candidate held-out comparison:
```text
logit(p0_ij) = alpha + beta_out log(1+kout_i)
              + beta_in log(1+kin_j) + beta_path Cij + gamma^T xij
logit(p1_ij) = logit(p0_ij) + architecture_features
DeltaS = S1 - S0
```
Baseline xij includes justified role/context and machinery/overlap features. Define a native-assertion decoder, risk set, sampling and score; compare on identical targets, train/calibrate on discovery/development. Unrecorded links are not validated negatives; prediction concerns curated knowledge. Use historical features only for temporal targets.

Null ladder: N0 degree and operation size; N1 adds role, molecular class, compartment, source and eligibility; N2 adds broad pathway constraints; N3 finer pathway/evidence constraints where feasible. Overlapping memberships need full profiles or pathway-conditioned incidence constraints; single-pathway swaps may violate other memberships. Audit mixing, accepted swaps, diversity, constraints and attainable resolution. Do not weaken nulls post hoc. Randomize native incidence before projection; block imported and near-duplicate assertions and shared evidence. Hold out native evidence groups before any projection, removing their derived links; never rely on random projected-edge splits.
For valid upper-tail randomization p=(1+#null>=observed)/(B+1). Repeat relevant selection under null when testing selected statistics on the same graph; otherwise freeze on training and test independent holdout. Pilot cost and mixing before choosing B; define test families, effect sizes, error control and stopping rules before confirmatory tests. Assess split/merge-aware stability, evidence resampling and representation sensitivity. If constraints make a test infeasible, limit the claim.

## 5. Corroboration and temporal prediction
Freeze Reactome assemblies and correspondence rules before scoring eligible nonimported GO-CAM causal assertions. Candidate T_GO=fraction of mapped assertions whose operations share a frozen assembly. Compare matched pathway/degree/context baselines; account for mapping uncertainty, unmatched assertions, direction, assembly size/overlap, model/publication dependence. Partial corroboration, not independent experiment.
Plan historical release discovery, separate development where feasible and later evaluation. Use separate versioned manifests; no overwrite of v98 lock. Choose release windows after bounded archive/schema feasibility without outcome tuning. v98 aggregates already inspected: retrospective, not untouched. Future unexamined releases preferred. Separate added assertions on known nodes, new operations, corrections, splits/merges and deletions. New nodes need a frozen assignment rule using information disjoint from evaluation targets; otherwise restrict to known nodes. Historical mappings/pathways/degrees only in training. New curation is not necessarily new biology. No additional acquisition yet.

## 6. Execution gates
| Stage | Deliverable |
|---|---|
| 1B | Certified harmonized catalogue, candidates, scope and evidence/benchmark metadata |
| 1C | Native typed graph and distinct candidate projection layers |
| 1D | Connectivity, hub, domain, pathway and null-feasibility diagnostics; bounded model/assembly/test design frozen before confirmatory fitting |
| 1E | Structural role inference and separate witnessed assembly inference |
| 1F | H02/H05 tests, robustness, temporal and source corroboration where feasible |
| 1G | Molecular boundaries, required/alternative mechanisms, context, evidence and uncertainty; unresolved architectures allowed |
| 1H | Freeze modules and assemblies separately, witnesses, gene links, support, pathway agreement, parameters and hashes |

Only then independently infer expression programmes (Part 2), test convergence with donor-level replication and disjoint-gene coordination (Part 3), test measured perturbation phenotypes beyond existing pathways and technical baselines (Part 4), and assess specialization and intercellular cooperation (Parts 5–6). Expression proxies are not measured activity; colocalization is not causation. Exact downstream models remain open.
Gate failure requires documented review, not threshold relaxation. No new inference, numerical cutoffs or statistical results are certified.

## References
- Peixoto, hierarchical blocks: https://doi.org/10.1103/PhysRevX.4.011047
- Peixoto, layered networks: https://doi.org/10.1103/PhysRevE.92.042807
- Ahn et al., overlapping link communities: https://doi.org/10.1038/nature09182
- PAC, pathway communities: https://doi.org/10.1093/nargab/lqac044
- Hart et al., cellular task inference: https://doi.org/10.1038/nmeth.3254
