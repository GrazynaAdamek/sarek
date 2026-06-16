# Technical Requirements: Fulfillment Summary

This document walks through every deliverable in the
[TECHNICAL RECRUITMENT TASK](TECHNICAL%20RECRUITMENT%20TASK.md) and describes how each one was satisfied, with pointers to the relevant files.

---

## PART 1 — Pipeline implementation

**Goal:** Running `--tools deepvariant --joint_genotype` (optionally with `vep`) must produce one multi-sample, VEP-annotated VCF for the whole cohort.

---

### Deliverable 1 — Nextflow process module for joint genotyping with container, inputs, outputs, and versions.yml

**Satisfied by:** [`modules/nf-core/glnexus/main.nf`](../../modules/nf-core/glnexus/main.nf)

The tool chosen is **GLnexus** (`glnexus_cli --config DeepVariant`), the standard companion
tool for DeepVariant gVCFs. The rationale for this choice over GATK's `GenomicsDBImport` /
`GenotypeGVCFs` is detailed in [JOINT_GENOTYPING_IMPLEMENTATION.md §2](JOINT_GENOTYPING_IMPLEMENTATION.md#2-why-glnexus).

The module is the official `nf-core/modules` GLnexus module, installed via
`nf-core modules install glnexus` and then locally patched (diff recorded in
[`modules/nf-core/glnexus/glnexus.diff`](../../modules/nf-core/glnexus/glnexus.diff)).

| Aspect | Detail |
|---|---|
| **Container** | `community.wave.seqera.io/library/bcftools_glnexus:cf380f1a6410f606` (Wave/Seqera; same image serves Docker and Singularity) |
| **Inputs** | `tuple val(meta), path(gvcfs), path(tbis), path(custom_config)` — gVCF files, their `.tbi` indexes (staged so GLnexus can find them), optional config; `tuple val(meta2), path(bed)` — optional interval BED |
| **Outputs** | `tuple val(meta), path("*.bcf"), emit: bcf` — merged cohort BCF; `path "versions.yml", emit: versions` — classic versions.yml; topic `versions_glnexus` (Nextflow's newer versions mechanism, picked up automatically by `workflows/sarek/main.nf`) |
| **versions.yml** | Written in both `script:` and `stub:` blocks; also emitted via the topic mechanism — both paths report the GLnexus version |

Patches applied over upstream (documented in `glnexus.diff` and [§5 of JOINT_GENOTYPING_IMPLEMENTATION.md](JOINT_GENOTYPING_IMPLEMENTATION.md#5-new-files)):

- `path(tbis)` input added — GLnexus requires index files present on disk next to gVCFs.
- `ulimit -n 65536` — prevents default open-file-limit exhaustion at ~1000 samples.
- `rm -rf GLnexus.DB` — makes retries safe (GLnexus refuses to start if the scratch DB already exists).
- Classic `versions.yml` emit added to satisfy this deliverable alongside the topic-based output.

`--config DeepVariant` is supplied via `ext.args` in
[`conf/modules/deepvariant_joint_genotype.config`](../../conf/modules/deepvariant_joint_genotype.config),
keeping the module close to upstream.

---

### Deliverable 2 — Subworkflow integrating joint genotyping across all cohort samples

**Satisfied by:** [`subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf`](../../subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf)

Full design description: [JOINT_GENOTYPING_IMPLEMENTATION.md §4 & §2](JOINT_GENOTYPING_IMPLEMENTATION.md#2-new-subworkflow-bam_joint_calling_germline_deepvariant-scattergather).

The subworkflow implements a **per-interval scatter / gather** pattern:

1. All samples' per-interval gVCFs (already split by the calling intervals) are grouped by
   interval with `groupTuple()`, so each GLnexus invocation sees only one interval's gVCFs
   for the whole cohort.
2. `GLNEXUS` runs once per interval (no `--bed` needed — inputs are already interval-restricted).
3. `BCFTOOLS_VIEW` converts each per-interval BCF to a compressed VCF.
4. If there are multiple intervals: `MERGE_GLNEXUS_VCF` (aliased `GATK4_MERGEVCFS`) gathers
   the per-interval VCFs into one cohort VCF.
5. If there is only one interval: `TABIX_TABIX` indexes it directly.
6. Output meta is rewritten to `{ id: joint_variant_calling, patient: all_samples, variantcaller: deepvariant }` for compatibility with downstream channels.

Key correctness detail: `intervals_name` is stripped from the grouping key before
`groupTuple()` in the gather step — without this, every interval would remain a separate
group and the merge would never fire (see the comment in the subworkflow at line 41–44).

The subworkflow takes:
- `gvcf_tbi_intervals` — per-interval, per-sample DeepVariant gVCFs (exposed by the patched
  `BAM_VARIANT_CALLING_DEEPVARIANT` subworkflow, see Deliverable 3).
- `dict` — sequence dictionary for `GATK4_MERGEVCFS`.

---

### Deliverable 3 — Integration into the VARIANT_CALLING workflow via `--joint_genotype`

**Satisfied by changes across four files:**

| File | Change |
|---|---|
| [`nextflow.config:88`](../../nextflow.config#L88) | `joint_genotype = false` default parameter added next to `joint_germline` |
| [`nextflow_schema.json:446`](../../nextflow_schema.json#L446) | Boolean schema entry in the `variant_calling` group, with `help_text` noting the `--tools deepvariant` requirement |
| [`nextflow.config:780`](../../nextflow.config#L780) | `includeConfig 'conf/modules/deepvariant_joint_genotype.config'` added next to the DeepVariant config |
| [`workflows/sarek/main.nf:428`](../../workflows/sarek/main.nf#L428) | `params.joint_genotype` passed into `BAM_VARIANT_CALLING_GERMLINE_ALL` alongside `params.joint_germline` |
| [`subworkflows/local/bam_variant_calling_germline_all/main.nf:49,123`](../../subworkflows/local/bam_variant_calling_germline_all/main.nf#L123) | New `joint_genotype` boolean `take:`; inside the DeepVariant block, if true, calls `BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT` and overwrites `vcf_deepvariant` / `tbi_deepvariant` with the joint output |
| [`subworkflows/local/bam_variant_calling_deepvariant/main.nf:31,77,97`](../../subworkflows/local/bam_variant_calling_deepvariant/main.nf#L31) | Added `gvcf_tbi` and `gvcf_tbi_intervals` emits so the per-interval, per-sample gVCFs are exposed to the joint calling subworkflow |
| [`subworkflows/local/utils_nfcore_sarek_pipeline/main.nf:251,272`](../../subworkflows/local/utils_nfcore_sarek_pipeline/main.nf#L251) | `jointGenotypeWithoutDeepvariant()` validation — logs a warning if `--joint_genotype` is set without `--tools deepvariant`, following the same convention as other flag-combination checks |

Overwriting `vcf_deepvariant` before it is mixed into `vcf_all` is the key design decision:
it means every downstream consumer (`POST_VARIANTCALLING`, the annotation channel) receives
the single multi-sample VCF automatically, without any further changes.  This mirrors exactly
how `--joint_germline` swaps in GATK's joint output.

---

### Deliverable 4 — Integration with the ANNOTATE subworkflow for a VEP-annotated multi-sample VCF

**Satisfied with zero changes to the annotation subworkflow.**

Because Deliverable 3 replaces `vcf_deepvariant` upstream of `vcf_all`, the joint cohort VCF
flows through the unchanged `VCF_ANNOTATE_ALL` subworkflow (called at
[`workflows/sarek/main.nf:562`](../../workflows/sarek/main.nf#L562)) exactly like any other
VCF.  Running with `--tools deepvariant,vep --joint_genotype` automatically produces:

```
annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz
```

The expected output paths for all flag combinations are summarised in
[JOINT_GENOTYPING_IMPLEMENTATION.md §7](JOINT_GENOTYPING_IMPLEMENTATION.md#7-end-to-end-behavior).

---

### Deliverable 5 — nf-test tests for the module and subworkflow

**Three test files:**

#### Module test
**File:** [`tests/modules/glnexus/main.nf.test`](../../tests/modules/glnexus/main.nf.test)

The test is placed under `tests/modules/` (not `modules/nf-core/glnexus/tests/`) so that it
sits **outside** the `nf-test.config` ignore glob (`modules/nf-core/**/tests/*`) and is
therefore executed by the standard runner.  This is necessary because the module has been
locally patched.

Test cases:

| Case | Mode | What it exercises |
|---|---|---|
| `vcfs, []` | real | baseline merge without indexes, bed, or custom config |
| `vcfs + tbis, []` | real | the new `tbis` input (staged `.tbi` indexes) |
| `vcfs, [], custom_config` | real | optional `custom_config` input |
| `vcfs, bed` | real | optional `--bed` interval restriction |
| `vcfs, bed - stub` | `-stub` | stub wiring |

Run / regenerate snapshot:
```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/modules/glnexus/ \
    --profile debug,test,docker --update-snapshot
```

#### Subworkflow test
**File:** [`subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test`](../../subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test)

| Test | Mode | What it validates |
|---|---|---|
| `deepvariant joint genotyping - no intervals` | real | GLNEXUS → BCFTOOLS_VIEW → TABIX_TABIX; VCF content hash (`variantsMD5`) is deterministic |
| `deepvariant joint genotyping - multiple intervals - stub` | `-stub` | scatter → group → merge wiring; `intervals_name` stripped before `groupTuple()` so both per-interval VCFs collapse into **one** merged output |

Run:
```bash
NXF_SYNTAX_PARSER=v1 nf-test test \
    subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test \
    --profile debug,test,docker --update-snapshot

# Run a second time (without --update-snapshot) to confirm determinism
NXF_SYNTAX_PARSER=v1 nf-test test \
    subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test \
    --profile debug,test,docker
```

#### Pipeline-level test
**File:** [`tests/joint_calling_deepvariant.nf.test`](../../tests/joint_calling_deepvariant.nf.test)

End-to-end pipeline test mirroring the analogous HaplotypeCaller test:
```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/joint_calling_deepvariant.nf.test \
    --profile debug,test,docker --verbose
```

---

### Deliverable 6 — Small test profile that runs locally end-to-end

**Two test profiles provided:**

#### Mini-genome profile (CI-speed, no external data)
**File:** [`conf/test_joint_genotyping.config`](../../conf/test_joint_genotyping.config)

Reuses the existing mini-genome fixtures and `tests/csv/3.0/mapped_joint_bam.csv`
(2 samples).  VEP excluded (mini-genome has no annotation cache).

```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping,docker \
    --outdir results_jg_mini
```

> **Profile naming note:** the task phrases this as `-profile test,docker`; the actual
> invocation is `-profile test,test_joint_genotyping,docker`.  The extra profile enables
> `--joint_genotype` additively, leaving the shared `test` profile unchanged so the
> existing test suite is unaffected.

Expected output: `results_jg_mini/variant_calling/deepvariant/joint_variant_calling/joint_variant_calling.vcf.gz`

#### Realistic 1000 Genomes profile (3 samples, chr20, with VEP)
**File:** [`conf/test_joint_genotyping_1000g.config`](../../conf/test_joint_genotyping_1000g.config)
**Data script:** [`scripts/prepare_testdata_1000g_chr20.sh`](../../scripts/prepare_testdata_1000g_chr20.sh)

Demonstrates the full deliverable: VEP-annotated multi-sample VCF from real WGS alignments.
Test data is **not committed** to the repository (per the task requirement) — the script
downloads and subsets 1000 Genomes phase-3 chr20 BAMs for three GBR individuals
(HG00096, HG00097, HG00099).

```bash
# Step 1: generate test data (requires samtools, ~2 GB download)
bash scripts/prepare_testdata_1000g_chr20.sh

# Step 2: run the pipeline
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping_1000g,docker \
    --outdir results_jg_1000g
```

Expected output: `results_jg_1000g/annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz`

VEP cache: uncomment `download_cache = true` in `conf/test_joint_genotyping_1000g.config`, or
set `vep_cache` to a pre-downloaded cache directory.

---

## PART 2 — Infrastructure challenge

**Goal:** A prototype Nextflow config for GCP or AWS addressing data-size challenges at
~1000-sample WGS scale (5 TB gVCFs, spot/preemptible instances).

**Satisfied by:**
- [`conf/cloud_aws.config`](../../conf/cloud_aws.config) — AWS Batch prototype
- [`conf/cloud_gcp.config`](../../conf/cloud_gcp.config) — GCP Batch prototype
- [`docs/intelliseq_technical_recruitment/CLOUD_SCALE_DESIGN.md`](CLOUD_SCALE_DESIGN.md) — design document

The design document addresses four challenges:

| Challenge | Solution |
|---|---|
| **Data volume** (5 TB of gVCFs to stage) | Seqera Fusion filesystem — mounts S3/GCS as a virtual POSIX filesystem, serving reads via HTTP range requests; no prior copy needed |
| **GLnexus memory at 1000-sample scale** | Per-interval scatter (already the architecture of Part 1) limits each job to ~200 GB of gVCFs and ~64 GB RAM; attempt-based memory escalation (`64.GB * task.attempt`) handles variance across chromosomes |
| **Spot/preemptible interruption** | Retry on preemption exit codes; `rm -rf GLnexus.DB` makes retries safe; AWS escalates to on-demand queue on final retry; GCP relies on Spot retries (no per-task on-demand escalation possible — documented as a platform limitation) |
| **Work directory hygiene** | `cleanup = true` for successful runs; S3/GCS lifecycle rules (30-day expiry) for failed/interrupted runs |

Usage (AWS):
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test_joint_genotyping_1000g,cloud_aws,docker \
    -work-dir s3://your-bucket/sarek-work \
    --outdir   s3://your-bucket/sarek-results \
    --aws_spot_queue     sarek-spot-queue \
    --aws_ondemand_queue sarek-ondemand-queue \
    --aws_region         eu-west-1
```

Usage (GCP):
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test_joint_genotyping_1000g,cloud_gcp,docker \
    -work-dir gs://your-bucket/sarek-work \
    --outdir   gs://your-bucket/sarek-results \
    --gcp_project your-gcp-project \
    --gcp_region  europe-west4
```

---

## Additional: Cohort QC report design

Beyond the task requirements, a cohort-level QC report design is provided in
[`COHORT_REPORT_DESIGN.md`](COHORT_REPORT_DESIGN.md).  It describes a phased approach:

- **Phase 1** (one config line): add `-s -` to `bcftools stats` on the joint VCF to generate
  per-sample statistics visible in the existing MultiQC report.
- **Phase 2**: a dedicated interactive HTML report (Rmarkdown/Quarto) with SFS, Ti/Tv per
  sample, call-rate distribution, and clinical variant tables (ClinVar / gnomAD).
