# Part 2 — Cloud-scale infrastructure for joint genotyping (~1000-sample WGS cohort)

**Goal:** A prototype Nextflow config for GCP or AWS addressing data-size challenges at ~1000-sample WGS scale (5 TB gVCFs, spot/preemptible instances).

**Satisfied by:**
- [`conf/cloud_aws.config`](../../conf/cloud_aws.config) — AWS Batch prototype
- [`conf/cloud_gcp.config`](../../conf/cloud_gcp.config) — GCP Batch prototype

The per-interval scatter/gather design from Part 1 is the foundation: instead of one job touching all 5 TB of gVCFs, each parallel GLnexus job touches only ~200 GB (one chromosome's gVCFs for the whole cohort).

| Challenge | Solution |
|---|---|
| **Data volume** (5 TB gVCFs) | Seqera Fusion — mounts S3/GCS as virtual POSIX filesystem via HTTP range requests; no prior copy |
| **GLnexus memory** | Per-interval scatter limits each job to ~64 GB RAM; attempt-based escalation (`64.GB * task.attempt`) handles large chromosomes |
| **Spot interruption** | Retry on preemption exit codes; `rm -rf GLnexus.DB` makes retries safe; AWS escalates to on-demand on final retry |
| **Work directory hygiene** | `cleanup = true` on success; S3/GCS lifecycle rules (30-day expiry) for failed runs |

---

## 1. Scale facts

| Resource | Per sample | 1 000 samples |
|---|---|---|
| CRAM | 20–30 GB | 20–30 TB |
| gVCF | ~5 GB | ~5 TB |
| gVCF per interval (≈ per chromosome) | ~200 MB | ~200 GB per GLnexus job |

---

## 2. Challenge 1 — Data volume: staging 5 TB of gVCFs

Without Fusion, Nextflow copies every input file into the task work directory first: 24 concurrent GLnexus jobs × 200 GB = ~4.8 TB simultaneous S3 egress, plus large instance disks and slow job start times.

**Seqera Fusion** mounts S3/GCS as a virtual POSIX filesystem. GLnexus sees local paths; reads are served via HTTP range requests directly from object storage — no copy, and only the positions GLnexus actually needs are fetched. Disk is used only for the `GLnexus.DB` scratch (~100 GB).

```groovy
fusion { enabled = true; exportStorageCredentials = true }
wave   { enabled = true; endpoint = 'https://wave.seqera.io' }
```

Requires a Seqera Platform account (free tier). Prefer instance types with local SSD — Fusion caches reads on local disk and suffers on slow storage. Without Fusion, set `disk = 250.GB` per GLnexus job for standard staging.

---

## 3. Challenge 2 — GLnexus memory at 1 000-sample scale

| Scope | Approximate RAM |
|---|---|
| Full genome, single job | 300–600 GB |
| Per-chromosome scatter | 30–60 GB |
| Large chromosomes (chr1, chr2) | up to 120 GB |

Per-chromosome scatter keeps each job within `r5.4xlarge` (AWS, 128 GB, ~$0.25/hr Spot) or `n2-highmem-16` (GCP, 128 GB, ~$0.20/hr Spot). Attempt-based escalation handles variance without per-chromosome tuning:

```groovy
memory = { 64.GB * task.attempt }   // 64 → 128 → 192 → 256 GB across attempts 1–4
```

**Open file descriptors.** GLnexus holds all input gVCF files open simultaneously. With 1 000 samples this exceeds the OS default limit of 1 024; `ulimit -n 65536` (added in Part 1) raises it to 65 536. This is a known architectural limitation — at >10 000 samples, [**Hail**](https://hail.is/) (`gvcf_combiner`) is worth evaluating: it uses a distributed merge that avoids the file-handle constraint, at the cost of requiring a Spark cluster.

---

## 4. Challenge 3 — Spot/preemptible interruption

A GLnexus job on a large chromosome may run 2–4 hours, making preemption a real risk. The configs retry on preemption exit codes and leave real failures (`exit 1`) to fail fast:

```groovy
// AWS  (104 = transient connection reset)
errorStrategy = { task.exitStatus in [104, 130, 137, 143] ? 'retry' : 'finish' }

// GCP
errorStrategy = { task.exitStatus in [14, 50001] ? 'retry' : 'finish' }
```

AWS escalates to an on-demand queue on the final retry; GCP has no per-task equivalent — a persistent capacity drought requires rerunning with `google.batch.spot = false` (use `-resume` to reuse completed work).

---

## 5. Challenge 4 — Work directory hygiene

At 1 000 samples, intermediate files (per-interval BCFs, VCFs) can reach tens of TB if left uncleaned.

- **`cleanup = true`** — deletes work subdirectories after a successful run (disables `-resume` for that run).
- **Object storage lifecycle rules** — delete work objects after 30 days, covering failed or abandoned runs:

```
AWS: S3 lifecycle rule  — Prefix: sarek-work/  — Expire after 30 days
GCP: GCS lifecycle rule — Prefix: sarek-work/  — Delete after 30 days
```

---

## 6. How to run

The cloud configs are additive (executor, Fusion, resource overrides). Supply input and reference paths pointing to S3/GCS — the local `test_joint_genotyping_1000g` profile uses `${projectDir}/tests/data/` paths and cannot be used directly on cloud.

```bash
# AWS Batch
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile cloud_aws,docker \
    -work-dir s3://your-bucket/sarek-work \
    --outdir   s3://your-bucket/sarek-results \
    --aws_spot_queue sarek-spot-queue --aws_ondemand_queue sarek-ondemand-queue \
    --aws_region eu-west-1 \
    --input s3://your-bucket/samples.csv --fasta s3://your-bucket/ref/genome.fasta \
    --tools deepvariant --joint_genotype

# GCP Batch
NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile cloud_gcp,docker \
    -work-dir gs://your-bucket/sarek-work \
    --outdir   gs://your-bucket/sarek-results \
    --gcp_project your-gcp-project --gcp_region europe-west4 \
    --input gs://your-bucket/samples.csv --fasta gs://your-bucket/ref/genome.fasta \
    --tools deepvariant --joint_genotype
```

---

## 7. Out of scope

Items identified but not implemented — see [`FUTURE_WORK.md`](FUTURE_WORK.md) for details:

- **Scatter VEP annotation per interval** — VEP currently runs once on the merged cohort VCF (`VCF_ANNOTATE_ALL`, unchanged). At scale this is a bottleneck: millions of variant sites, hours of serial work, no resumability on preemption. The correct fix is to annotate the per-interval VCFs (the `vcf_out.intervals` channel that exists before `MERGE_GLNEXUS_VCF`) and merge annotated results afterward — VEP annotates each record independently by position, so annotate-then-merge is result-equivalent. This requires either restructuring the shared `VCF_ANNOTATE_ALL` subworkflow or adding a joint-path-specific annotate-then-merge, both out of scope here.
- **GPU acceleration for DeepVariant** — dominant wall-clock cost at 1 000-sample scale.
- **Terraform / CDK for compute environment provisioning** — AWS Batch queues, IAM roles, S3 buckets, GCP Batch setup.
