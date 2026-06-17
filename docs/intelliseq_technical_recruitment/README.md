# Intelliseq Technical Recruitment — nf-core/sarek extension

## Documents

| File | Contents |
|---|---|
| [`PART1_JOINT_GENOTYPING.md`](PART1_JOINT_GENOTYPING.md) | Implementation walkthrough: GLnexus module, scatter/gather subworkflow, pipeline integration, tests, and test profiles |
| [`PART2_CLOUD_INFRASTRUCTURE.md`](PART2_CLOUD_INFRASTRUCTURE.md) | Cloud-scale design for ~1000-sample WGS cohorts: Fusion, memory, spot interruption, work directory hygiene |
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

### 5. End-to-end run — 1000 Genomes profiles (3 samples, real WGS data)

Real-life complement to the mini-genome profile above: same pipeline path, but on real WGS data, so runtime is longer (minutes vs. seconds) and disk usage is higher (BAMs + reference + VEP cache vs. a few KB).

Prepare test data once (downloads ~1.4 GB chr20 BAMs and reference):

```bash
bash scripts/prepare_testdata_1000g_chr20.sh
```

Two profiles with distinct purposes:

| Profile | Tests | VEP |
|---|---|---|
| `test_joint_genotyping_1000g` | DeepVariant + GLnexus joint genotyping on full chr20 | no |
| `test_joint_genotyping_1000g_intervals` | Scatter/gather (3 intervals) + joint genotyping + annotation | yes — requires `--vep_cache` |

```bash
# Joint genotyping module (no VEP):
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping_1000g,docker \
    --outdir results_jg_1000g

# Full end-to-end including VEP — download the cache once, then pass --vep_cache:
bash scripts/prepare_testdata_1000g_chr20.sh --vep /path/to/vep_cache

NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test,test_joint_genotyping_1000g_intervals,docker \
    --vep_cache /path/to/vep_cache \
    --outdir results_jg_1000g_intervals
```
