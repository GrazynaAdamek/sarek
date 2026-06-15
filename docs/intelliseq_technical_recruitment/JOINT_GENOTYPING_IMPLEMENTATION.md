# Joint Genotyping for the DeepVariant Path (GLnexus)

## 1. Goal

Running sarek with:

```bash
--tools deepvariant --joint_genotype
```

(optionally with `vep` added to `--tools`) takes the per-sample DeepVariant
gVCFs produced for a cohort and merges/genotypes them into **one multi-sample
VCF**, which then flows through the existing annotation subworkflow to produce
a single VEP-annotated, multi-sample VCF for the whole cohort.

## 2. Why GLnexus

DeepVariant calls each sample independently and (optionally) emits a gVCF per
sample, but sarek has no step that combines these into a cohort-level VCF.
GATK's `--joint_germline` already solves this for HaplotypeCaller via
GenomicsDBImport + GenotypeGVCFs + VQSR — but that pipeline is GATK-specific
and overkill for DeepVariant gVCFs.

**GLnexus** (`glnexus_cli --config DeepVariant`) is the standard companion
tool for DeepVariant: it merges per-sample gVCFs directly into one
multi-sample BCF in a single step, with a config preset tuned for
DeepVariant's gVCF conventions. No GenomicsDB step is needed. There is no
nf-core module for it, so a small local module was added
(`modules/local/glnexus/`), following the conventions of other "no-conda"
binary tool modules in this repo (e.g. `deepvariant/rundeepvariant`).

## 3. High-level architecture

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
                          ├── GLNEXUS  (--config DeepVariant, --bed <interval>)
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

## 4. Scatter/gather design (the scalability decision)

**Problem:** the target cohort is ~1000 WGS samples. A naive design would run
GLnexus *once*, on the whole genome, across all 1000 gVCFs — a single huge,
long-running, memory-hungry, non-resumable job. That's a poor fit for
production (one OOM/preemption near the end re-runs everything) and for cost
on spot/preemptible instances.

**Solution:** mirror the scatter/gather pattern DeepVariant itself already
uses for variant calling, reusing the *same* `intervals` channel:

- `gvcf_tbi` from all samples is grouped into one cohort-wide list
  (`meta = [id: 'joint_variant_calling']`) and then `.combine(intervals)`'d,
  so GLnexus runs **once per interval** (e.g. once per chromosome), each
  invocation processing the *whole cohort* but restricted to that region via
  `--bed`.
- Each per-interval result (BCF → VCF via `BCFTOOLS_VIEW`) is then gathered
  back into a single cohort VCF with `GATK4_MERGEVCFS` (aliased
  `MERGE_GLNEXUS_VCF`) — the exact same merge module DeepVariant itself uses
  for `MERGE_DEEPVARIANT_VCF`/`MERGE_DEEPVARIANT_GVCF`.
- If there are no intervals (`num_intervals <= 1`), GLnexus runs once on the
  whole genome and the BCFTOOLS_VIEW output is indexed directly with
  `TABIX_TABIX` — no merge step needed.

This means with, say, 24 intervals (one per chromosome), 24 GLnexus jobs run
in parallel, each handling all ~1000 samples but only for its chromosome,
instead of one job processing the entire genome for 1000 samples at once.
Failures/preemptions only cost one interval's worth of work and retry
(`task.attempt`) independently.

## 5. New files

### `modules/local/glnexus/main.nf`

```groovy
process GLNEXUS {
    ...
    container "quay.io/mlin/glnexus:v1.3.1"
    input:
    tuple val(meta), path(gvcfs), path(tbis), path(intervals)
    output:
    tuple val(meta), path("${prefix}.bcf"), emit: bcf
    path "versions.yml", emit: versions
    ...
}
```

Key implementation details and why:

- **`path(tbis)`** — the gVCF `.tbi` indexes aren't referenced on the command
  line, but GLnexus needs them present on disk next to the gVCFs, so they're
  declared as a staged input.
- **`--bed ${intervals}`** (conditional) — enables the per-interval scatter
  described above; omitted entirely when no intervals are used.
- **`ulimit -n 65536`** — GLnexus opens every sample's gVCF simultaneously;
  at ~1000 samples the default open-file limit (1024) would be exceeded.
- **`--mem-gbytes` = 90% of `task.memory`** — GLnexus's internal cache size is
  tied to the task's memory allocation, with 10% headroom left for the
  process itself and htslib buffers (avoids OOM at the cache boundary).
- **gVCF manifest file (`gvcf.list`)** — paths are written to a file and
  expanded via `$(cat gvcf.list)` rather than inlined directly, for
  reproducibility/debuggability. At ~1000 samples (~50KB of paths) this is
  still well under typical `ARG_MAX` (~2MB), so no `xargs`/batching is needed.
- **`rm -rf GLnexus.DB`** — `glnexus_cli` refuses to start if its scratch DB
  directory already exists; this guard makes `task.attempt` retries safe.
- No `conda`/`environment.yml` — like `DEEPVARIANT_RUNDEEPVARIANT`, GLnexus is
  a complex statically-linked binary not packaged for conda; container-only,
  matching the existing precedent in this repo.

### `subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf`

```groovy
workflow BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT {
    take:
    gvcf_tbi   // [ meta, gvcf, tbi ] per sample
    dict
    intervals  // [ intervals, num_intervals ] or [ [], 0 ]

    main:
    // 1. group all samples into one cohort set, fan out across intervals
    // 2. GLNEXUS per interval (--bed restricted)
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

A subtle but important detail: the per-interval `meta.intervals_name` (used
for `ext.prefix` on `GLNEXUS`/`BCFTOOLS_VIEW`) must be **stripped before
grouping** for the merge step:

```groovy
vcf_to_merge = vcf_out.intervals
    .map{ meta, vcf -> [ groupKey(meta - meta.subMap('intervals_name'), meta.num_intervals), vcf ] }
    .groupTuple()
```

Without this, every interval would produce a distinct grouping key (because
`intervals_name` differs per interval), `groupTuple()` would never actually
group anything, and `MERGE_GLNEXUS_VCF` would run once per interval on
single-element lists instead of merging the cohort VCF back together.

## 6. Wiring changes

### `subworkflows/local/bam_variant_calling_deepvariant/main.nf`

Added a `gvcf_tbi` emit (the gVCF index), built the same way the existing
`tbi` emit is built for the VCF — via a new `gvcf_tbi_out` branch on
`DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index`, mixed with
`MERGE_DEEPVARIANT_GVCF.out.tbi` for the multi-interval case. Previously the
subworkflow emitted the gVCF itself but not its index, which the joint
genotyping subworkflow needs.

### `subworkflows/local/bam_variant_calling_germline_all/main.nf`

- New `take:` boolean `joint_genotype` (default `false`), alongside the
  existing `joint_germline`.
- Inside the existing `if (tools.contains('deepvariant'))` block, after the
  normal per-sample `vcf_deepvariant`/`tbi_deepvariant` are captured:

```groovy
if (joint_genotype) {
    gvcf_tbi_deepvariant = BAM_VARIANT_CALLING_DEEPVARIANT.out.gvcf
        .join(BAM_VARIANT_CALLING_DEEPVARIANT.out.gvcf_tbi, failOnMismatch: true)

    BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT(gvcf_tbi_deepvariant, dict, intervals)

    vcf_deepvariant = BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT.out.genotype_vcf
    tbi_deepvariant = BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT.out.genotype_index
    versions = versions.mix(BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT.out.versions)
}
```

This is exactly the same shape as how `joint_germline` swaps
`vcf_haplotypecaller`/`tbi_haplotypecaller` for the joint-genotyped GATK
output — by overwriting `vcf_deepvariant`/`tbi_deepvariant` *before* they're
mixed into `vcf_all`/`tbi_all`, every downstream consumer
(`POST_VARIANTCALLING`, `vcf_to_annotate`, `VCF_ANNOTATE_ALL`) automatically
receives the single multi-sample VCF instead of N per-sample VCFs — **no
changes needed in the annotation subworkflow at all**.

`dict` and `intervals` were already `take:` parameters of this subworkflow
(used to scatter DeepVariant itself), so they're simply forwarded — the joint
genotyping scatter reuses the exact same interval split as variant calling.

### `workflows/sarek/main.nf`

- Passes `params.joint_genotype` into `BAM_VARIANT_CALLING_GERMLINE_ALL`,
  alongside `params.joint_germline`.

### `subworkflows/local/utils_nfcore_sarek_pipeline/main.nf`

- Adds a `jointGenotypeWithoutDeepvariant()` validation function, called from
  `validateInputParameters()` alongside `genomeExistsError()` and
  `sparkAndBam()` — the existing convention for pipeline-level
  parameter-combination checks run during initialization:

```groovy
// Warn (not error) so shared/default param files setting joint_genotype=true don't break runs with other --tools
def jointGenotypeWithoutDeepvariant() {
    if (params.joint_genotype && !(params.tools && params.tools.split(',').contains('deepvariant'))) {
        log.warn("--joint_genotype is set but '--tools deepvariant' is not. Joint genotyping with GLnexus will not be run.")
    }
}
```

This catches the likely user error of setting `--joint_genotype` without
`--tools deepvariant` (the flag would otherwise silently do nothing, since
the new code path only runs inside the `deepvariant` tools branch).

### `nextflow.config` / `nextflow_schema.json`

- New param `joint_genotype = false`, placed next to `joint_germline`.
- New schema entry (boolean, `variant_calling` group), with `help_text`
  noting the `--tools deepvariant` requirement.
- New `includeConfig 'conf/modules/deepvariant_joint_genotype.config'`,
  placed next to `includeConfig 'conf/modules/deepvariant.config'`.

### `conf/modules/deepvariant_joint_genotype.config` (new)

Per-process overrides scoped to
`.*:BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT:*` (so they don't affect any other
use of `BCFTOOLS_VIEW`/`GATK4_MERGEVCFS`/`TABIX_TABIX` elsewhere in the
pipeline):

- **`GLNEXUS`**: `cpus = params.max_cpus`, `memory = 64.GB * task.attempt`,
  `ext.prefix` differentiates per-interval outputs
  (`joint_variant_calling.<interval>`) vs. the single-interval case
  (`joint_variant_calling`), `ext.when` gated on
  `tools.contains('deepvariant') && joint_genotype`, intermediate BCF not
  published.
- **`BCFTOOLS_VIEW`**: `--output-type z` (BCF → `vcf.gz`), same prefix
  scheme, not published (intermediate per-interval VCF).
- **`MERGE_GLNEXUS_VCF`**: prefix `joint_variant_calling`, published to
  `variant_calling/deepvariant/joint_variant_calling/`.
- **`TABIX_TABIX`**: published to the same directory (no-intervals case).

The `GLNEXUS` memory/cpu settings are the primary scalability knob for ~1000
WGS samples — bump per-interval memory as needed; `task.attempt` retry
escalation absorbs occasional OOMs on larger chromosomes (e.g. chr1/chr2).

## 7. End-to-end behavior

| Flags | Result |
|---|---|
| `--tools deepvariant` | unchanged: N per-sample VCFs |
| `--tools deepvariant --joint_genotype` | one multi-sample VCF at `variant_calling/deepvariant/joint_variant_calling/joint_variant_calling.vcf.gz` |
| `--tools deepvariant,vep --joint_genotype` | additionally, one VEP-annotated multi-sample VCF at `annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz` |
| `--joint_genotype` without `deepvariant` in `--tools` | warning logged, flag has no effect |

## 8. Known limitations / explicitly out of scope for this change

- **Tests**: nf-test coverage for the new module/subworkflow and a small
  multi-sample (`-profile test`) joint-genotyping test profile were
  deliberately deferred.
- **Cloud-scale infra (Part 2)**: GCP/AWS Batch-specific `nextflow.config`
  profile tuning (spot/preemptible handling, GCS/S3 staging of ~5GB gVCFs ×
  1000 samples, work-dir lifecycle) is a separate follow-up; the scatter/gather
  design above is the foundation that makes that follow-up tractable (failures
  are scoped to one interval, not the whole cohort).
