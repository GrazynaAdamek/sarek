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

This isn't just a stylistic preference: `GenomicsDBImport`/`GenotypeGVCFs` assume
GATK HaplotypeCaller's gVCF conventions (its genotype-likelihood model,
`<NON_REF>` symbolic allele handling, `PL`/`AD` field semantics). DeepVariant's
gVCFs encode genotype likelihoods from a different (CNN-based) model, so GATK's
joint-genotyping tools aren't designed or validated to merge them correctly.

**GLnexus** (`glnexus_cli --config DeepVariant`) is the standard companion
tool for DeepVariant: it merges per-sample gVCFs directly into one
multi-sample BCF in a single step, with a config preset tuned for
DeepVariant's gVCF conventions. No GenomicsDB step is needed. GLnexus's
`--config DeepVariant` preset is purpose-built for DeepVariant's gVCF
representation — using it (instead of GATK's tools) is a correctness
requirement (matching tool to gVCF format), not a stylistic choice.

The official `nf-core/modules` GLnexus module
(`modules/nf-core/glnexus/`) is used, installed via `nf-core modules install`
and locally patched (`modules/nf-core/glnexus/glnexus.diff`, applied via
`nf-core modules patch`) to add the `.tbi` index staging and retry-safety
handling this pipeline needs (see section 5). `--config DeepVariant` is
supplied via `ext.args` in `conf/modules/deepvariant_joint_genotype.config`
rather than hardcoded in the module, keeping the module close to upstream.

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

## 4. Scatter/gather design (the scalability decision)

**Problem:** the target cohort is ~1000 WGS samples. A naive design would run
GLnexus *once*, on the whole genome, across all 1000 gVCFs — a single huge,
long-running, memory-hungry, non-resumable job. That's a poor fit for
production (one OOM/preemption near the end re-runs everything) and for cost
on spot/preemptible instances.

**Solution:** mirror the scatter/gather pattern `joint_germline` already uses
for HaplotypeCaller — group the cohort's **per-interval, per-sample** gVCFs by
interval, so each GLnexus invocation only ever sees the gVCFs for one region:

- `BAM_VARIANT_CALLING_DEEPVARIANT` exposes a new `gvcf_tbi_intervals` emit —
  the **per-interval** gVCF/tbi for each sample (before they get merged into
  the whole-genome `gvcf`/`gvcf_tbi` used elsewhere), built the same way
  HaplotypeCaller's `gvcf_tbi_intervals` is (joining the per-interval gVCF/tbi
  with `cram_intervals`).
- `BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT` groups these by `intervals_name`
  across all samples (`groupTuple()`), so GLnexus runs **once per interval**
  (e.g. once per chromosome), each invocation processing only that region's
  gVCFs for the whole cohort — **no `--bed` needed**, since the input gVCFs
  are already restricted to that interval.
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
Because each job only stages the small per-interval gVCFs (not the
whole-genome gVCF restricted via `--bed`), there's no redundant staging of
each sample's full gVCF once per interval. Failures/preemptions only cost one
interval's worth of work and retry (`task.attempt`) independently.

## 5. New files

### `modules/nf-core/glnexus/main.nf` (+ `glnexus.diff`)

Installed from `nf-core/modules` (`nf-core modules install glnexus`) and
locally patched (`nf-core modules patch glnexus`, recorded in
`modules/nf-core/glnexus/glnexus.diff`):

```groovy
process GLNEXUS {
    ...
    container "${ ... 'community.wave.seqera.io/library/bcftools_glnexus:...' }"
    input:
    tuple val(meta), path(gvcfs), path(tbis), path(custom_config)
    tuple val(meta2), path(bed)
    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    tuple val("${task.process}"), val('glnexus'), eval(...), topic: versions, emit: versions_glnexus
    path "versions.yml", emit: versions
    ...
}
```

Patch details and why:

- **`path(tbis)`** (added) — the gVCF `.tbi` indexes aren't referenced on the
  command line, but GLnexus needs them present on disk next to the gVCFs, so
  they're declared as a staged input alongside `gvcfs`.
- **`ulimit -n 65536`** (added) — GLnexus opens every sample's gVCF
  simultaneously; at ~1000 samples the default open-file limit (1024) would
  be exceeded.
- **`rm -rf GLnexus.DB`** (added) — `glnexus_cli` refuses to start if its
  scratch DB directory already exists; this guard makes `task.attempt`
  retries safe.
- **`path "versions.yml", emit: versions`** (added) — writes and emits a
  classic `versions.yml`, redundant with `versions_glnexus` below but kept
  to satisfy the literal "...versions.yml" deliverable requirement (see
  note below).
- Everything else (container/Wave management, conda `environment.yml`,
  `--mem-gbytes` defaulting, optional `--bed` via the second input tuple, the
  topic-based `versions` output) is kept as shipped upstream, so
  `nf-core modules update` remains viable.
- `--config DeepVariant` is **not** in the module — it's supplied via
  `ext.args` in `conf/modules/deepvariant_joint_genotype.config`.
- The `versions_glnexus` output uses Nextflow's `topic: versions` mechanism;
  sarek's `workflows/sarek/main.nf` already does
  `versions.mix(channel.topic("versions"))`, so this is picked up
  automatically without any extra wiring in the subworkflow.
- Note: this `topic: versions` mechanism is newer than the classic
  `versions.yml` emit used by the other modules in this subworkflow
  (`BCFTOOLS_VIEW`, `MERGE_GLNEXUS_VCF`, `TABIX_TABIX`), which are still
  collected via explicit `versions.mix(MODULE.out.versions)` calls. Both are
  collected correctly by the pipeline; the topic-based approach is newer and
  may be adopted by other nf-core modules over time via
  `nf-core modules update`.
- As an additional patch on top of upstream, the module also writes a
  classic `versions.yml` file and emits it via `emit: versions` (mixed into
  `versions` in the subworkflow like the other modules). This is redundant
  with `versions_glnexus`/`topic: versions` above — both end up reporting
  the same `glnexus` version into the pipeline's aggregated versions report
  — but it directly satisfies the task's literal "...with container, inputs,
  outputs, and versions.yml" requirement without removing the newer
  topic-based mechanism the upstream module ships with.

### `subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf`

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

Added two new emits:

- `gvcf_tbi` — the gVCF index, built the same way the existing `tbi` emit is
  built for the VCF (via a new `gvcf_tbi_out` branch on
  `DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index`, mixed with
  `MERGE_DEEPVARIANT_GVCF.out.tbi` for the multi-interval case). Previously
  the subworkflow emitted the gVCF itself but not its index.
- `gvcf_tbi_intervals` — the **per-interval, per-sample** gVCF/tbi (before
  merging into the whole-genome `gvcf`/`gvcf_tbi`), built the same way
  HaplotypeCaller's `gvcf_tbi_intervals` is — joining
  `DEEPVARIANT_RUNDEEPVARIANT.out.gvcf`/`.out.gvcf_index` with
  `cram_intervals` (already computed for the scatter):

```groovy
gvcf_tbi_intervals = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
    .join(DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index, failOnMismatch: true)
    .join(cram_intervals, failOnMismatch: true)
    .map{ meta, gvcf, tbi, cram, crai, intervals -> [ meta, gvcf, tbi, intervals ] }
```

This is what the joint genotyping subworkflow consumes, so each GLnexus
invocation only ever sees gVCFs already restricted to one interval (see
section 4).

### `subworkflows/local/bam_variant_calling_germline_all/main.nf`

- New `take:` boolean `joint_genotype` (default `false`), alongside the
  existing `joint_germline`.
- Inside the existing `if (tools.contains('deepvariant'))` block, after the
  normal per-sample `vcf_deepvariant`/`tbi_deepvariant` are captured:

```groovy
if (joint_genotype) {
    BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT(
        BAM_VARIANT_CALLING_DEEPVARIANT.out.gvcf_tbi_intervals,
        dict
    )

    vcf_deepvariant = BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT.out.genotype_vcf
    tbi_deepvariant = BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT.out.genotype_index
    versions = versions.mix(BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT.out.versions)
}
```

This is exactly the same shape as how `joint_germline` swaps
`vcf_haplotypecaller`/`tbi_haplotypecaller` for the joint-genotyped GATK
output (which is fed directly from
`BAM_VARIANT_CALLING_HAPLOTYPECALLER.out.gvcf_tbi_intervals`, with no join at
the call site either) — by overwriting `vcf_deepvariant`/`tbi_deepvariant`
*before* they're mixed into `vcf_all`/`tbi_all`, every downstream consumer
(`POST_VARIANTCALLING`, `vcf_to_annotate`, `VCF_ANNOTATE_ALL`) automatically
receives the single multi-sample VCF instead of N per-sample VCFs — **no
changes needed in the annotation subworkflow at all**.

`dict` was already a `take:` parameter of this subworkflow, so it's simply
forwarded.

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

### Annotation at cohort scale (potential future optimization)

In the current design, joint genotyping is scattered/gathered per interval,
but **annotation is not**: `MERGE_GLNEXUS_VCF` gathers the per-interval VCFs
into one whole-genome multi-sample VCF first, and VEP then runs **once** on
that single merged file (via the unchanged shared `VCF_ANNOTATE_ALL`
subworkflow). This is intentional for this task — it directly produces the
required *single* "multi-sample VEP-annotated VCF" with zero changes to the
shared annotation subworkflow, and runs fine on the small test data.

At ~1000-sample WGS scale this single VEP job is a potential bottleneck.
Worth noting precisely *why*:

- **Parsing the large VCF is not itself the problem.** VEP streams the VCF
  record-by-record; it does not load the whole file into memory, so a very
  large multi-sample VCF won't exhaust RAM just from being read.
- The real costs of one giant annotation job are:
  1. **Wall-clock as a single serial job** — a cohort VCF can contain tens of
     millions of variant *sites*, and VEP throughput (even forked) is on the
     order of thousands of records/sec, so this can run for hours.
  2. **No resumability on preemption** — on spot/preemptible instances, a
     long non-checkpointed job that is killed near the end re-does all of its
     work. This is the same reliability/cost concern that motivated scattering
     GLnexus in the first place (section 4).
  3. **Object-store I/O** — Nextflow must localize the whole merged VCF from
     GCS/S3 to the worker before VEP starts and delocalize the result
     afterward; one huge file means one large serial transfer plus
     bgzip/tabix of a large output.

Possible approaches (not implemented here):

- **Scatter VEP per interval, gather after.** VEP annotates each record
  independently by genomic position against the cache/reference — there is no
  cross-record or cross-chunk context in the standard consequence annotation
  sarek runs. So annotating the per-interval VCFs (the `vcf_out.intervals`
  channel that already exists *before* `MERGE_GLNEXUS_VCF`) and merging the
  annotated results afterward is **result-equivalent** to annotating the merged
  file (intervals partition the genome and never split a variant). This reuses
  the existing GLnexus scatter and gives parallelism, resumability, and
  parallel localization. The cost: VEP has real per-invocation overhead (cache
  region load, fork setup, container start), so the chunks should be the
  calling intervals (dozens), not thousands of tiny pieces, to keep per-job
  runtime dominant over startup. Structurally, this would mean either
  restructuring the shared `VCF_ANNOTATE_ALL` to scatter internally (affects
  *every* sarek path, not just the joint path) or adding a joint-path-specific
  annotate-then-merge — both larger than this task's "small test data" scope,
  which is why it's deferred.
- **More CPUs for the joint-path VEP job (low-effort lever).** VEP's built-in
  multithreading (`--fork`) parallelizes annotation across CPUs within a single
  job. Note that sarek's module **already** runs `--fork ${task.cpus}`
  (`modules/nf-core/ensemblvep/vep/main.nf`), so this is not an unused switch to
  flip — `--fork` automatically tracks the task's allocated CPUs. The actual
  lever is therefore raising `cpus` for the cohort VEP job: today
  `withName: 'ENSEMBLVEP_VEP'` in `conf/modules/annotate.config` sets `ext.args`
  and `publishDir` but does not override `cpus`, so the huge cohort VCF is
  annotated with the same default resources as a tiny single-sample VCF. A
  scoped override (e.g. `withName: '.*:VCF_ANNOTATE_ALL:ENSEMBLVEP_VEP'` in
  `deepvariant_joint_genotype.config`, so it doesn't bloat per-sample
  annotation) bumping `cpus` would raise `--fork` with it. This is the cheapest
  win short of restructuring and needs **no pipeline topology change**, but it
  improves wall-clock only — not resumability (a single forked job that is
  preempted still loses all its work).
- **Trim VEP options/plugins** to only what's needed — each plugin adds
  per-record cost — and right-size cpus/memory so `--fork` has cores to use.
- **Disable VEP summary stats (`--no_stats`) for the cohort job.** sarek
  passes `--stats_file` (`conf/modules/annotate.config:37`), so VEP aggregates
  a per-run summary across *all* records in memory and writes an HTML report.
  On a cohort-scale VCF this aggregation is a real memory/runtime cost;
  `--no_stats` (scoped to the joint-path VEP job) cuts both. Tradeoff: you lose
  the per-run VEP summary HTML and its MultiQC panel for that file — so this is
  a deliberate trade, not a free win.
- **`--buffer_size` tuning (minor lever).** VEP annotates variants in buffers
  (default 5000) before each batch lookup; a larger buffer can improve
  throughput (better cache locality, fewer round-trips) at some memory cost.
  Set via `ext.args`. Low impact, listed for completeness.
- **The same single-job cost applies to the other annotation tools, not just
  VEP.** `VCF_ANNOTATE_ALL` runs the merged cohort VCF through the `snpeff`,
  `bcfann` (`BCFTOOLS_ANNOTATE`), and `merge` branches the same way — each is a
  single job on the whole merged file. The per-interval scatter remedy above
  would benefit all of them; the resource/args levers (`--fork`/cpus,
  `--no_stats`, buffer size) are per-tool.

Recommendation: keep the current single-VEP-on-merged-VCF as the Part 1
deliverable (correct, minimal, reuses the shared subworkflow), and treat the
scatter-annotation design above as a Part 2 scaling item — Part 2 explicitly
accepts a written design over an implementation.
