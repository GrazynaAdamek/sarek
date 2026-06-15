# Joint Genotyping for DeepVariant (GLnexus) + VEP annotation

## Context

Sarek's DeepVariant path (`BAM_VARIANT_CALLING_DEEPVARIANT`) currently calls each
sample independently and emits per-sample VCFs (it also produces gVCFs internally
but they are discarded/not exposed). There is no way to combine DeepVariant calls
across a cohort into one multi-sample VCF. GATK's `--joint_germline` already
implements this pattern for HaplotypeCaller (GenomicsDBImport + GenotypeGVCFs),
giving us a template to follow.

The standard companion tool for joint-genotyping DeepVariant gVCFs is **GLnexus**
(`glnexus_cli --config DeepVariant`), which merges per-sample gVCFs directly into
one multi-sample BCF — no GenomicsDB step needed. GLnexus has no nf-core module,
so we add a small local module.

Goal: `--tools deepvariant --joint_genotype` (with `--tools deepvariant,vep`)
produces one multi-sample, VEP-annotated VCF for the whole cohort, following the
same wiring pattern as `joint_germline`.

## 1. New local module: GLNEXUS

`modules/local/glnexus/main.nf` (mirrors structure/conventions of
`modules/nf-core/deepvariant/rundeepvariant/main.nf`, but lives under `local/`
since no nf-core module exists). It accepts an optional BED file so it can be
scattered across genomic intervals, exactly like `DEEPVARIANT_RUNDEEPVARIANT`:

```groovy
process GLNEXUS {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/mlin/glnexus:v1.4.1"

    input:
    tuple val(meta), path(gvcfs), path(tbis), path(intervals)

    output:
    tuple val(meta), path("${prefix}.bcf"), emit: bcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def bed    = intervals ? "--bed ${intervals}" : ''
    """
    ulimit -n 65536

    glnexus_cli \\
        --config DeepVariant \\
        --threads ${task.cpus} \\
        --mem-gbytes ${task.memory.toGiga()} \\
        ${bed} \\
        ${args} \\
        ${gvcfs} > ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(glnexus_cli --version 2>&1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(glnexus_cli --version 2>&1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """
}
```

Add `modules/local/glnexus/meta.yml` describing inputs/outputs (gvcfs list,
tbi list, optional intervals bed, output bcf) following the format of other
local module meta.yml files.

GLnexus needs the per-sample gVCF indexes present on disk (read via
`path(tbis)` so Nextflow stages them) even though they're not referenced in the
command line. The `ulimit -n` bump accounts for ~1000 simultaneously-open gVCF
files; `--mem-gbytes` ties GLnexus's internal cache size to the task's memory
allocation so it doesn't overrun the node.

## 2. New subworkflow: BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT (scatter/gather)

New file `subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf`,
modeled on both `bam_joint_calling_germline_gatk/main.nf` (joint-calling pattern)
and `bam_variant_calling_deepvariant/main.nf` (scatter/gather over `intervals`).
GLnexus runs once per interval across the whole cohort (restricted to that
region via `--bed`), and per-interval joint BCFs are converted to VCF and
merged back together — same shape as `MERGE_DEEPVARIANT_VCF`:

```groovy
include { GLNEXUS                                } from '../../../modules/local/glnexus/main'
include { BCFTOOLS_VIEW                          } from '../../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX                            } from '../../../modules/nf-core/tabix/tabix/main'
include { GATK4_MERGEVCFS as MERGE_GLNEXUS_VCF   } from '../../../modules/nf-core/gatk4/mergevcfs/main'

workflow BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT {
    take:
    gvcf_tbi   // channel: [ meta, gvcf, tbi ]  (per-sample DeepVariant gVCFs)
    dict       // channel: [ meta, dict ]
    intervals  // channel: [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Group all samples into one cohort-wide list, then fan out across intervals
    glnexus_input = gvcf_tbi
        .map{ meta, gvcf, tbi -> [ [ id:'joint_variant_calling' ], gvcf, tbi ] }
        .groupTuple()
        .combine(intervals)
        .map{ meta, gvcf, tbi, intervals_, num_intervals ->
            [ meta + [ num_intervals:num_intervals, intervals_name: intervals_ ? intervals_.baseName : null ], gvcf, tbi, intervals_ ?: [] ]
        }

    GLNEXUS(glnexus_input)

    // BCF -> compressed VCF (per interval, or whole genome if no intervals)
    BCFTOOLS_VIEW(GLNEXUS.out.bcf.map{ meta, bcf -> [ meta, bcf, [] ] }, [], [], [])

    vcf_out = BCFTOOLS_VIEW.out.vcf.branch{
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when scattered across intervals: merge per-interval VCFs back into one.
    // `intervals_name` differs per interval, so it must be stripped from the
    // grouping key first - otherwise every interval gets its own group of size 1
    // and MERGE_GLNEXUS_VCF never actually merges anything across intervals.
    vcf_to_merge = vcf_out.intervals
        .map{ meta, vcf -> [ groupKey(meta - meta.subMap('intervals_name'), meta.num_intervals), vcf ] }
        .groupTuple()

    MERGE_GLNEXUS_VCF(vcf_to_merge, dict)

    // Single-interval / no-intervals case needs its own index
    TABIX_TABIX(vcf_out.no_intervals)

    vcf_merged = Channel.empty().mix(MERGE_GLNEXUS_VCF.out.vcf, vcf_out.no_intervals)
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'deepvariant' ], vcf ] }

    tbi_merged = Channel.empty().mix(MERGE_GLNEXUS_VCF.out.tbi, TABIX_TABIX.out.tbi)
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'deepvariant' ], tbi ] }

    versions = versions.mix(GLNEXUS.out.versions)
    versions = versions.mix(BCFTOOLS_VIEW.out.versions)
    versions = versions.mix(MERGE_GLNEXUS_VCF.out.versions)
    versions = versions.mix(TABIX_TABIX.out.versions)

    emit:
    genotype_vcf   = vcf_merged
    genotype_index = tbi_merged
    versions
}
```

This means with, say, 24 intervals (one per chromosome), 24 GLnexus jobs run in
parallel — each handling all ~1000 samples but only for its chromosome — instead
of one job processing the entire genome for 1000 samples at once.

## 3. Expose DeepVariant gVCF + wire joint genotyping in BAM_VARIANT_CALLING_GERMLINE_ALL

**a) `subworkflows/local/bam_variant_calling_deepvariant/main.nf`**
Already emits `gvcf` but no matching index — add a `gvcf_tbi` emit. This needs
a *new* branch on `DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index` (mirroring the
existing `tbi_out` branch on `.out.vcf_index` at lines 52-56), since the
no-intervals gVCF index comes from `gvcf_index`, not from `gvcf_out` (which
holds the gVCF file itself, not its index):

```groovy
// Figuring out if there is one or more gvcf index(es) from the same sample
gvcf_tbi_out = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index.branch{
    intervals:    it[0].num_intervals > 1
    no_intervals: it[0].num_intervals <= 1
}

gvcf_tbi = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.tbi, gvcf_tbi_out.no_intervals)
    .map{ meta, tbi -> [ meta - meta.subMap('num_intervals') + [ variantcaller:'deepvariant' ], tbi ] }
```

Add `gvcf_tbi` to the `emit:` block.

**b) `subworkflows/local/bam_variant_calling_germline_all/main.nf`**
- Add `include { BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT } from '../bam_joint_calling_germline_deepvariant/main'`
- Add new `take:` param `joint_genotype` (boolean, default false)
- In the DEEPVARIANT block (around line 108-120), after capturing `vcf_deepvariant`/`tbi_deepvariant`, add:

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

`dict` and `intervals` are already `take:` parameters of
`BAM_VARIANT_CALLING_GERMLINE_ALL`, so no new inputs are needed at this level —
they're just forwarded to the new subworkflow (same `intervals` channel used to
scatter DeepVariant itself, so the joint genotyping scatter uses the same
chromosome/interval split).

This mirrors exactly how `joint_germline` swaps `vcf_haplotypecaller`/`tbi_haplotypecaller`
for the joint-genotyped versions, so `vcf_all`/`tbi_all` (and therefore
`POST_VARIANTCALLING` and the annotate channel) automatically carry the
single multi-sample VCF instead of N per-sample VCFs.

## 4. Wire the new flag through workflows/sarek/main.nf

- In `workflows/sarek/main.nf`, find the call to `BAM_VARIANT_CALLING_GERMLINE_ALL`
  and add `params.joint_genotype` as the new argument matching the new `take:`
  parameter (placed alongside the existing `params.joint_germline` argument).
- Add a validation warning (alongside the other `--tools`-related checks in
  `workflows/sarek/main.nf` / `subworkflows/local/utils_nfcore_sarek_pipeline`)
  that logs a warning (and is a no-op) if `--joint_genotype` is set but
  `--tools` does not contain `deepvariant` — mirrors how other tool-specific
  flags are validated, and clarifies the "requires `--tools deepvariant`" note
  already in the schema `help_text`.

## 5. New CLI flag + defaults

- `nextflow.config`: add `joint_genotype = false` near `joint_germline`/`joint_mutect2`
  (line ~88-89), with a short inline comment.
- `nextflow_schema.json`: add a `joint_genotype` boolean entry next to `joint_germline`
  (same `variant_calling` section, ~line 440), e.g.:

```json
"joint_genotype": {
    "type": "boolean",
    "fa_icon": "fas fa-toolbox",
    "description": "Turn on joint genotyping for DeepVariant using GLnexus",
    "help_text": "Merges per-sample DeepVariant gVCFs across the whole cohort into a single multi-sample VCF using GLnexus. Requires `--tools deepvariant`."
}
```

## 6. Module config

New `conf/modules/deepvariant_joint_genotype.config`, included from `nextflow.config`
next to `includeConfig 'conf/modules/joint_germline.config'`:

```groovy
process {
    withName: '.*:BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT:GLNEXUS' {
        cpus       = { params.max_cpus }
        memory     = { 64.GB * task.attempt }
        ext.prefix = { meta.num_intervals <= 1 ? 'joint_variant_calling' : "joint_variant_calling.${meta.intervals_name}" }
        publishDir = [ enabled: false ]
    }

    withName: '.*:BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT:BCFTOOLS_VIEW' {
        ext.args   = '--output-type z'
        ext.prefix = { meta.num_intervals <= 1 ? 'joint_variant_calling' : "joint_variant_calling.${meta.intervals_name}" }
        publishDir = [ enabled: false ]
    }

    withName: '.*:BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT:MERGE_GLNEXUS_VCF' {
        ext.prefix = { 'joint_variant_calling' }
        publishDir = [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/variant_calling/deepvariant/joint_variant_calling/" },
            pattern: "*{vcf.gz,vcf.gz.tbi}"
        ]
    }

    withName: '.*:BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT:TABIX_TABIX' {
        publishDir = [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/variant_calling/deepvariant/joint_variant_calling/" },
            pattern: "*.tbi"
        ]
    }
}
```

The `GLNEXUS` memory/cpu overrides are the main scalability knob for ~1000 WGS
samples — bump per-interval as needed (`task.attempt` retry escalation handles
occasional OOMs on larger chromosomes). Per-interval intermediate files
(`BCFTOOLS_VIEW`, raw `GLNEXUS` bcf) are not published; only the final merged
cohort VCF/index are.

Also extend `conf/modules/deepvariant.config`'s `MERGE_DEEPVARIANT_GVCF` block if a
prefix tweak is needed for the gvcf tbi (likely not — `GATK4_MERGEVCFS` already
emits `.tbi` alongside `.vcf` with the same prefix).

## 7. Annotation (no new code needed)

Because the joint VCF replaces `vcf_deepvariant` before it's mixed into `vcf_all`,
it flows through `POST_VARIANTCALLING` → `vcf_to_annotate` → `VCF_ANNOTATE_ALL`
exactly like any other VCF (same as the `joint_germline` GATK case). Running with
`--tools deepvariant,vep --joint_genotype` will therefore automatically produce
`annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz` as the
final multi-sample annotated output. No changes needed to `vcf_annotate_all` or
`workflows/sarek/main.nf`'s annotation block.

## Files to touch (summary)

- `modules/local/glnexus/main.nf` (new)
- `modules/local/glnexus/meta.yml` (new)
- `subworkflows/local/bam_joint_calling_germline_deepvariant/main.nf` (new)
- `subworkflows/local/bam_variant_calling_deepvariant/main.nf` (add `gvcf_tbi` emit)
- `subworkflows/local/bam_variant_calling_germline_all/main.nf` (new `joint_genotype` take + DEEPVARIANT block)
- `workflows/sarek/main.nf` (pass `params.joint_genotype`)
- `nextflow.config` (new param default + includeConfig)
- `nextflow_schema.json` (new schema entry)
- `conf/modules/deepvariant_joint_genotype.config` (new)

## Verification

- `nextflow run main.nf -profile test,docker --tools deepvariant --joint_genotype --outdir results_jg` (multi-sample test profile, e.g. `test_cache`/`test`) and confirm:
  - Per-sample DeepVariant gVCFs are produced.
  - `GLNEXUS` runs (once per interval, or once overall if no intervals), producing per-interval `joint_variant_calling[.<interval>].bcf`.
  - A single `joint_variant_calling/joint_variant_calling.vcf.gz` (+ `.tbi`) appears containing all sample genotypes (check `bcftools query -l` lists all samples).
  - With `--tools deepvariant,vep --joint_genotype`, a VEP-annotated multi-sample VCF is produced under `annotation/vep/joint_variant_calling/`.
- Run `nf-core modules lint` / `nf-core subworkflows lint` (or `nf-core lint`) on the new module/subworkflow for nf-core compliance.
- `nextflow config -profile test` to confirm schema/param wiring has no errors.
