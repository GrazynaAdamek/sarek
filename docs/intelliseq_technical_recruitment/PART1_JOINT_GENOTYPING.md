# Part 1 — Extend nf-core/sarek with joint genotyping for a large WGS cohort

Running `--tools deepvariant --joint_genotype` produces one multi-sample, optionally VEP-annotated VCF for the whole cohort.

**Without annotation:**
```bash
--tools deepvariant --joint_genotype
```
Outputs one multi-sample VCF at `variant_calling/deepvariant/joint_variant_calling/`.

**With VEP annotation:**
```bash
--tools deepvariant,vep --joint_genotype
```
Additionally produces a VEP-annotated VCF at `annotation/vep/joint_variant_calling/`.
`vep` must be added explicitly to `--tools`; `--joint_genotype` does not imply it
(doing so would conflict with `--tools ...,snpeff` and fail if no VEP cache is configured).

---

## High-level architecture

```
per-sample CRAM ──► DEEPVARIANT_RUNDEEPVARIANT (scatter per interval)
                          │
                          ├── per-sample VCF  ──► (existing path, unchanged)
                          └── per-sample gVCF + index
                                   │
                                   ▼  (only if --joint_genotype)
                     BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT
                          │
                          ├── group ALL samples' gVCFs into one cohort set
                          ├── fan out across the same `intervals` used for
                          │   variant calling (scatter)
                          ├── GLNEXUS  (--config DeepVariant, per-interval gVCFs, no --bed)
                          ├── BCFTOOLS_VIEW   (BCF → VCF.gz)
                          ├── MERGE_GLNEXUS_VCF (GATK4 MergeVcfs, gather)
                          │   or TABIX_TABIX if only one interval
                          ▼
                  one multi-sample VCF (id: joint_variant_calling,
                                         patient: all_samples,
                                         variantcaller: deepvariant)
                          │
                          ▼ (replaces vcf_deepvariant in vcf_all)
                  POST_VARIANTCALLING ──► VCF_ANNOTATE_ALL (VEP)
                          │
                          ▼
        annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz
```

---

## Deliverable 1 — Nextflow process module for joint genotyping

**File:** [`modules/nf-core/glnexus/main.nf`](../../modules/nf-core/glnexus/main.nf)

### Why GLnexus

Sarek had no step to combine DeepVariant's per-sample gVCFs into a cohort VCF. GATK's `GenomicsDBImport`/`GenotypeGVCFs` — already implemented for HaplotypeCaller in `subworkflows/local/bam_joint_calling_germline_gatk` and triggered via `--joint_germline` - can't be reused here: they assume HaplotypeCaller's gVCF conventions (`<NON_REF>` allele, PL/AD semantics) and are not validated against DeepVariant's CNN-based likelihood model.

**GLnexus** (`glnexus_cli --config DeepVariant`) is the correct tool: it merges DeepVariant gVCFs into a multi-sample BCF in one step with a preset tuned for DeepVariant's format.

Following nf-core best practices, the GLnexus module is sourced from `nf-core/modules` and minimally patched (`modules/nf-core/glnexus/glnexus.diff`) only for pipeline-specific needs (`.tbi` index staging, retry-safety). Tool arguments (`--config DeepVariant`) are kept out of the module and set via `ext.args` in `conf/modules/deepvariant_joint_genotype.config`, keeping the module close to upstream and easy to update.

### Module interface

| Aspect | Detail |
|---|---|
| **Container** | `community.wave.seqera.io/library/bcftools_glnexus:cf380f1a6410f606` (Wave/Seqera; same image serves Docker and Singularity) |
| **Inputs** | `tuple val(meta), path(gvcfs), path(tbis), path(custom_config)` — gVCF files, their `.tbi` indexes (staged so GLnexus can find them), optional config; `tuple val(meta2), path(bed)` — optional interval BED |
| **Outputs** | `tuple val(meta), path("*.bcf"), emit: bcf` — merged cohort BCF; `path "versions.yml", emit: versions` — classic versions.yml; topic `versions_glnexus` (Nextflow's newer versions mechanism, picked up automatically by `workflows/sarek/main.nf`) |
| **versions.yml** | Written in both `script:` and `stub:` blocks; also emitted via the topic mechanism — both paths report the GLnexus version |

### Patches over upstream

Patches are documented in [`modules/nf-core/glnexus/glnexus.diff`](../../modules/nf-core/glnexus/glnexus.diff):

- **`path(tbis)` input added** — the gVCF `.tbi` indexes aren't referenced on the command line, but GLnexus needs them present on disk next to the gVCFs.
- **`ulimit -n 65536`** — GLnexus opens every sample's gVCF simultaneously; at ~1000 samples the default open-file limit (1024) would be exceeded.
- **`rm -rf GLnexus.DB`** — `glnexus_cli` refuses to start if its scratch DB directory already exists; this guard makes `task.attempt` retries safe.
- **`path "versions.yml", emit: versions`** — writes and emits a classic `versions.yml`, redundant with `versions_glnexus` (topic-based) but added to satisfy the literal versions.yml deliverable requirement alongside the newer mechanism.

`--config DeepVariant` is supplied via `ext.args` in
[`conf/modules/deepvariant_joint_genotype.config`](../../conf/modules/deepvariant_joint_genotype.config), keeping the module close to upstream.

---

## Deliverable 2 — Subworkflow integrating joint genotyping across all cohort samples

**File:** [`subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf`](../../subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf)

### Scatter/gather design

**Problem:** the target cohort is ~1000 WGS samples. A naive design would run GLnexus *once*, on the whole genome, across all 1000 gVCFs — a single huge, long-running, memory-hungry, non-resumable job. That's a poor fit for production (one OOM/preemption near the end re-runs everything) and for cost on spot/preemptible instances.

**Solution:** mirror the scatter/gather pattern `joint_germline` already uses for HaplotypeCaller - group the cohort's **per-interval, per-sample** gVCFs by interval, so each GLnexus invocation only ever sees the gVCFs for one region:

- `BAM_VARIANT_CALLING_DEEPVARIANT` exposes a new `gvcf_tbi_intervals` emit — the **per-interval** gVCF/tbi for each sample (before they get merged into the whole-genome `gvcf`/`gvcf_tbi`), built the same way HaplotypeCaller's `gvcf_tbi_intervals` is. This is an additive change: existing `gvcf` and `gvcf_tbi` emits are untouched, so no downstream processes are affected.
- `BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT` groups these by `intervals_name` across all samples (`groupTuple()`), so GLnexus runs **once per interval** — each invocation processing only that region's gVCFs for the whole cohort.
- Each per-interval result (BCF → VCF via `BCFTOOLS_VIEW`) is gathered back into a single cohort VCF with `GATK4_MERGEVCFS` (aliased `MERGE_GLNEXUS_VCF`) — the same merge module DeepVariant itself uses for `MERGE_DEEPVARIANT_VCF`/`MERGE_DEEPVARIANT_GVCF`.
- If there is only one interval (`num_intervals <= 1`): GLnexus runs once on the whole genome and the `BCFTOOLS_VIEW` output is indexed directly with `TABIX_TABIX` - no merge step needed.

With 24 intervals (one per chromosome), 24 GLnexus jobs run in parallel, each handling all ~1000 samples but only for its chromosome. Failures/preemptions only cost one interval's worth of work and retry (`task.attempt`) independently.

### Subworkflow interface

```groovy
workflow BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT {
    take:
    gvcf_tbi_intervals  // [ meta, gvcf, tbi, intervals ] per-interval, per sample
    dict

    main:
    // 1. group all samples' per-interval gVCFs by interval (cohort-wide)
    // 2. GLNEXUS per interval (no --bed needed, inputs already per-interval)
    // 3. BCFTOOLS_VIEW: BCF -> VCF.gz
    // 4. branch on num_intervals:
    //      >1  -> MERGE_GLNEXUS_VCF (GATK4 MergeVcfs) gathers per-interval VCFs
    //      <=1 -> TABIX_TABIX indexes the single VCF directly
    // 5. remap meta -> { id: joint_variant_calling, patient: all_samples,
    //                     variantcaller: deepvariant }
    emit:
    genotype_vcf
    genotype_index
    versions
}
```

### Critical correctness detail: `intervals_name` stripping

`intervals_name` must be stripped from `meta` before `groupTuple()` — otherwise each interval has a distinct key, the tuple never collapses, and `MERGE_GLNEXUS_VCF` receives one-element lists instead of the full per-interval set:

---

## Deliverable 3 — Integration into the VARIANT_CALLING workflow via `--joint_genotype`

**Key design choice: overwrite, don't branch.** `BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT` replaces `vcf_deepvariant`/`tbi_deepvariant` in place — before they are mixed into `vcf_all` — so every downstream consumer (annotation, QC, publishing) receives the joint cohort VCF with zero changes to those consumers.

**Changed files** ([`nextflow.config`](../../nextflow.config), [`nextflow_schema.json`](../../nextflow_schema.json), [`workflows/sarek/main.nf`](../../workflows/sarek/main.nf), [`bam_variant_calling_germline_all/main.nf`](../../subworkflows/local/bam_variant_calling_germline_all/main.nf), [`bam_variant_calling_deepvariant/main.nf`](../../subworkflows/local/bam_variant_calling_deepvariant/main.nf), [`utils_nfcore_sarek_pipeline/main.nf`](../../subworkflows/local/utils_nfcore_sarek_pipeline/main.nf)) plus the new [`conf/modules/deepvariant_joint_genotype.config`](../../conf/modules/deepvariant_joint_genotype.config).

**New emits on `bam_variant_calling_deepvariant`** — `gvcf_tbi` (whole-genome gVCF index, previously missing) and `gvcf_tbi_intervals` (per-interval gVCF + index before merge, mirroring HaplotypeCaller's existing emit). Both are additive; existing emits are unchanged.

**`deepvariant_joint_genotype.config`** scopes all process overrides to `.*:BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT:*`, so `BCFTOOLS_VIEW`, `TABIX_TABIX`, and `GATK4_MERGEVCFS` used elsewhere in the pipeline are unaffected. GLnexus memory uses attempt-based escalation (`64.GB * task.attempt`); intermediate per-interval BCFs and VCFs are not published.


---

## Deliverable 4 — Integration with the ANNOTATE subworkflow for a VEP-annotated multi-sample VCF

**Satisfied with zero changes to the annotation subworkflow.**

Because Deliverable 3 replaces `vcf_deepvariant` upstream of `vcf_all`, the joint cohort VCF flows through the unchanged `VCF_ANNOTATE_ALL` subworkflow (called at [`workflows/sarek/main.nf:562`](../../workflows/sarek/main.nf#L562)) exactly like any other VCF.

### End-to-end behavior

| Flags | Result |
|---|---|
| `--tools deepvariant` | unchanged: N per-sample VCFs |
| `--tools deepvariant --joint_genotype` | one multi-sample VCF at `variant_calling/deepvariant/joint_variant_calling/joint_variant_calling.vcf.gz` |
| `--tools deepvariant,vep --joint_genotype` | additionally, one VEP-annotated multi-sample VCF at `annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz` |
| `--joint_genotype` without `deepvariant` in `--tools` | warning logged, flag has no effect |

---

## Deliverable 5 — nf-test tests for the module and subworkflow

### Module test

**File:** [`tests/modules/glnexus/main.nf.test`](../../tests/modules/glnexus/main.nf.test)

Based on the upstream nf-core module test, relocated to `tests/modules/` so it is picked up by CI. The upstream location (`modules/nf-core/glnexus/tests/`) is excluded by the `modules/nf-core/**/tests/*` ignore glob in `nf-test.config` — necessary because the locally patched module diverges from upstream and its test must actually run.

| Case | Mode | What it exercises |
|---|---|---|
| `vcfs, []` | real | baseline merge without indexes, bed, or custom config |
| `vcfs + tbis, []` | real | the new `tbis` input (staged `.tbi` indexes) |
| `vcfs, [], custom_config` | real | optional `custom_config` input |
| `vcfs, bed` | real | optional `--bed` interval restriction |
| `vcfs, bed - stub` | `-stub` | stub wiring |

Each case asserts `process.success`, snapshots `bcf` + `versions_glnexus`, and checks the legacy `versions` (`versions.yml`) output.

```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/modules/glnexus/ \
    --profile debug,test,docker
```

### Subworkflow test

**File:** [`subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test`](../../subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test)

| Test | Mode | What it validates |
|---|---|---|
| `deepvariant joint genotyping - no intervals` | real | GLNEXUS → BCFTOOLS_VIEW → TABIX_TABIX; VCF content hash (`variantsMD5`) is deterministic |
| `deepvariant joint genotyping - multiple intervals - stub` | `-stub` | scatter → group → merge wiring; `intervals_name` stripped before `groupTuple()` so both per-interval VCFs collapse into **one** merged output |

**Determinism note:** Test 1 snapshots VCF content via `variantsMD5` (the nft-vcf plugin), not a raw file md5. `BCFTOOLS_VIEW` writes a `Date=` line into the VCF header on every run, making raw file md5 non-deterministic. `tests/lib/UTILS.groovy:60` establishes this as the repo-wide convention for all VCF snapshots. The companion `tests/nextflow.config` in the same directory is required: it mirrors `conf/modules/deepvariant_joint_genotype.config` to give per-interval outputs distinct prefixes, preventing stage-in collisions into `MERGE_GLNEXUS_VCF`.

Run (first pass creates the snapshot; second pass confirms determinism):
```bash
NXF_SYNTAX_PARSER=v1 nf-test test \
    subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test \
    --profile debug,test,docker --update-snapshot

NXF_SYNTAX_PARSER=v1 nf-test test \
    subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test \
    --profile debug,test,docker
```

### Pipeline-level test

**File:** [`tests/joint_calling_deepvariant.nf.test`](../../tests/joint_calling_deepvariant.nf.test)

End-to-end pipeline test mirroring the analogous HaplotypeCaller test.

| Test | Mode | What it validates |
|---|---|---|
| `nucleotides_per_second 20` | real | full pipeline with scatter (multiple intervals); snapshots output files and VCF content |
| `nucleotides_per_second 100` | real | same pipeline with different interval sizing; confirms output is deterministic across runs |
| `--joint_genotype without --tools deepvariant` | `-stub` | pipeline succeeds and produces no joint VCF — `--joint_genotype` is a no-op (warning only) when deepvariant is not in `--tools` |

```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/joint_calling_deepvariant.nf.test \
    --profile debug,test,docker --verbose
```

---

## Deliverable 6 — Small test profile that runs locally end-to-end

### Mini-genome profile (CI-speed, no external data)

**File:** [`conf/test_joint_genotyping.config`](../../conf/test_joint_genotyping.config)

Follows the same pattern as the existing `test` and `test_joint_germline` profiles. Reuses the existing mini-genome fixtures and `tests/csv/3.0/mapped_joint_bam.csv` (2 samples).

```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping,docker \
    --outdir results_jg_mini
```

**Profile naming rationale:** The specification requested a profile named `test`, but that name is already taken by the existing mini-genome CI profile. To avoid conflict, the joint-genotyping run is provided as a separate composable profile (`test_joint_genotyping`) layered on top of `test` — an additive, non-breaking change that leaves the existing test suite unaffected.

Expected output: `results_jg_mini/variant_calling/deepvariant/joint_variant_calling/joint_variant_calling.vcf.gz`

### Realistic 1000 Genomes profile (3 samples, chr20, with VEP)

**File:** [`conf/test_joint_genotyping_1000g.config`](../../conf/test_joint_genotyping_1000g.config)
**Data script:** [`scripts/prepare_testdata_1000g_chr20.sh`](../../scripts/prepare_testdata_1000g_chr20.sh)

Demonstrates the full deliverable: VEP-annotated multi-sample VCF from real WGS alignments. Test data is **not committed** to the repository — the script downloads and subsets 1000 Genomes phase-3 chr20 BAMs for three GBR individuals (HG00096, HG00097, HG00099).

**Chr20 rationale:** mid-sized (~63 Mb in GRCh37), gene-rich enough to produce real variants for VEP to annotate, and small enough (~20–50 MB per sample after subsetting) to keep the test tractable without dedicated data infrastructure.

```bash
# Step 1: generate test data (requires samtools, ~2 GB download)
bash scripts/prepare_testdata_1000g_chr20.sh

# Step 2: run the pipeline
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping_1000g,docker \
    --outdir results_jg_1000g
```

VEP cache is read from `s3://annotation-cache/vep_cache/` by default. To use a local cache instead, pass `--vep_cache /path/to/cache`.

Expected output: `results_jg_1000g/annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz`

---

## Future addition: cohort QC report

A dedicated cohort-level QC report would be a valuable addition for production use. The pipeline already runs `bcftools stats` and vcftools on the joint VCF and feeds results into MultiQC — adding `--samples -` to `bcftools stats` would immediately surface per-sample variant counts, Ti/Tv, and het/hom ratios as a sortable table in the existing MultiQC HTML (one config line change). A second phase could produce a standalone interactive HTML report (Rmarkdown/Quarto) with an allele frequency spectrum, per-sample outlier plots, and a table of clinically significant variants cross-referenced against gnomAD and ClinVar.

---

## Known technical debt

**Nextflow v1 → v2 syntax migration.** All commands require `NXF_SYNTAX_PARSER=v1` because Nextflow ≥ 24.x defaults to the v2 parser, which rejects v1-style constructs still present in the pipeline. The env var is a temporary compatibility shim — the pipeline needs a full syntax audit (`nf-core pipelines lint` flags most issues) before `NXF_SYNTAX_PARSER=v1` can be dropped. New code added here (GLnexus module, `BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT`) is written to be v2-compatible.

**Topic-based `versions.yml`.** The pipeline mixes the classic `path "versions.yml", emit: versions` with the newer Nextflow topic channel (`emit: topic: 'versions'`). The GLnexus module emits both in parallel as a transitional measure. The long-term migration is to replace all classic emits and their `mix`/`dump_software_versions` wiring with topic-based emission — removing boilerplate from every module and subworkflow.
