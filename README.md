# Cellular Functional Architecture (CFA)

**Goal:** discover complementary molecular assemblies beyond existing pathway organization, then test whether they correspond to independently derived cellular programmes and measured functions in human macrophages, fibroblasts and endothelial cells.

**Status:** Part 1A source lock approved; Part 1B frozen with explicit approval; Part 1C in progress. Inference direction v0.3.0 approved; algorithms and confirmatory thresholds not frozen.

## Design
1. Preserve molecular states, operations, typed roles, composition and provenance.
2. Infer structural role modules from Reactome alone, excluding curated precedence and pathway membership from primary fitting.
3. Assemble complementary modules using operation-level molecular interface witnesses.
4. Test incremental information and complementarity against pathway-conditioned baselines and feasible nulls.
5. Evaluate later-release curation and eligible nonimported GO-CAM corroboration; freeze candidate architectures.
6. Test independent expression convergence and functional perturbation evidence, then cellular specialization and tissue coordination.

A statistical module is not a function. Curated topology is not function-blind. Part 1 does not establish physiological activity or causality.

## Essential files
| File | Purpose |
|---|---|
| [Mathematical protocol](protocol/PART1_MATHEMATICAL_PROTOCOL.md) | Current scientific and statistical design |
| [Harmonization plan](protocol/PART1B_HARMONIZATION_PLAN.md) | Current implementation gates |
| [Source certification](reports/SOURCE_CERTIFICATION.md) | Verified inventory, findings, limitations and provenance |
| [Harmonization lock](manifests/P1B_HARMONIZATION_LOCK_V1.json) | Approved certificate and upstream fingerprints |
| [Source lock](manifests/P1A_SOURCE_LOCK_V1.json) | Immutable approved source and historical policy fingerprints |
| [Source policy](protocol/P1A_SOURCE_REPRESENTATION_POLICY.md) | Approved state/identity representation (source scope only) |
| [Licensing](protocol/SOURCE_LICENSE_POLICY.md) | Data attribution |
| [Execution](protocol/SLURM_EXECUTION_POLICY.md) | Cluster requirements |
| [Scripts](scripts/README.md) | Runnable stages and retained dependencies |
| [Registries](registries/) | Current decisions, hypotheses, milestones and open limitations |

## Execution
Parts 1B.01–1B.04 passed; do not rerun. See the certification report for manifest fingerprints, counts and freeze scope. Next: Part 1C construction according to approved v0.3.0 design (native incidence first; separate candidate projections). Run substantial processing through the SLURM wrapper from repository root. Raw data, results and logs stay outside Git. The approved source lock relies on historical Git objects: use a full clone (not a shallow snapshot).
