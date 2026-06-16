# Intelliseq Technical Recruitment — nf-core/sarek extension

## Documents

| File | Contents |
|---|---|
| [`PART1_JOINT_GENOTYPING.md`](PART1_JOINT_GENOTYPING.md) | Implementation walkthrough: GLnexus module, scatter/gather subworkflow, pipeline integration, tests, and test profiles |
| [`PART2_CLOUD_INFRASTRUCTURE.md`](PART2_CLOUD_INFRASTRUCTURE.md) | Cloud-scale design for ~1000-sample WGS cohorts: Fusion, memory, spot interruption, work directory hygiene |
| [`COHORT_REPORT_DESIGN.md`](COHORT_REPORT_DESIGN.md) | Design proposal for a per-sample cohort QC report (phased; not implemented) |
| [`TECHNICAL RECRUITMENT TASK.md`](TECHNICAL%20RECRUITMENT%20TASK.md) | Original task specification |

---

## How to run

> All commands require `NXF_SYNTAX_PARSER=v1` — see [Known technical debt](PART1_JOINT_GENOTYPING.md#known-technical-debt) in Part 1.

### Using the new joint genotyping feature

Add `--joint_genotype` to any existing DeepVariant run:

```bash
# Joint VCF only
--tools deepvariant --joint_genotype

# Joint VCF + VEP annotation
--tools deepvariant,vep --joint_genotype
```

Output: `variant_calling/deepvariant/joint_variant_calling/joint_variant_calling.vcf.gz`
Annotated: `annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz`

---

### 1. Unit tests (GLnexus module)

```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/modules/glnexus/ \
    --profile debug,test,docker
```

### 2. Subworkflow test

First pass creates the snapshot; second pass confirms determinism:

```bash
NXF_SYNTAX_PARSER=v1 nf-test test \
    subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test \
    --profile debug,test,docker --update-snapshot

NXF_SYNTAX_PARSER=v1 nf-test test \
    subworkflows/local/bam_joint_calling_germline_deepvariant/tests/main.nf.test \
    --profile debug,test,docker
```

### 3. Pipeline-level test (stub + mini-genome fixtures, no external data)

```bash
NXF_SYNTAX_PARSER=v1 nf-test test tests/joint_calling_deepvariant.nf.test \
    --profile debug,test,docker --verbose
```

### 4. End-to-end run — mini-genome profile

```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping,docker \
    --outdir results_jg_mini
```

Expected output: `results_jg_mini/variant_calling/deepvariant/joint_variant_calling/joint_variant_calling.vcf.gz`

### 5. End-to-end run — realistic 1000 Genomes profile (3 samples, chr20, VEP)

```bash
# Step 1: download and subset test data (~2 GB, requires samtools)
bash scripts/prepare_testdata_1000g_chr20.sh

# Step 2: run the pipeline
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping_1000g,docker \
    --outdir results_jg_1000g
```

Expected output: `results_jg_1000g/annotation/vep/joint_variant_calling/joint_variant_calling_VEP.ann.vcf.gz`

VEP cache is read from `s3://annotation-cache/vep_cache/` by default. To use a local cache, pass `--vep_cache /path/to/cache`.
