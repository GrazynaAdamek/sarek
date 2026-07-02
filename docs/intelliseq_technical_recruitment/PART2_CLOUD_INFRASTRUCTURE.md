# Part 2 — Cloud-scale infrastructure for joint genotyping (~1000-sample WGS cohort)

**Goal:** A prototype Nextflow config for GCP or AWS addressing data-size challenges at ~1000-sample WGS scale (5 TB gVCFs, spot/preemptible instances).

**Satisfied by:**
- [`conf/cloud_aws.config`](../../conf/cloud_aws.config) — AWS Batch prototype
- [`conf/cloud_gcp.config`](../../conf/cloud_gcp.config) — GCP Batch prototype

The configs are additive — add `-profile cloud_aws,docker` (or `cloud_gcp`) with S3/GCS paths for `-work-dir`, `--outdir`, `--input`, and `--fasta`:

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

The per-interval scatter/gather design from Part 1 is the foundation: instead of one job touching all 5 TB of gVCFs, each parallel GLnexus job touches only ~200 GB (one chromosome's gVCFs for the whole cohort).

| Challenge | Solution |
|---|---|
| **Data volume** (5 TB gVCFs) | Seqera Fusion — mounts S3/GCS as virtual POSIX filesystem via HTTP range requests; no prior copy |
| **GLnexus memory** | Per-interval scatter limits each job to ~64 GB RAM; attempt-based escalation (`64.GB * task.attempt`) handles large chromosomes |
| **Spot interruption** | Retry on preemption exit codes; `rm -rf GLnexus.DB` makes retries safe; AWS escalates to on-demand on final retry |
| **Work directory hygiene** | `cleanup = true` on success; S3/GCS lifecycle rules (30-day expiry) for failed runs |

**VEP annotation is not yet scattered.** `VCF_ANNOTATE_ALL` runs once on the merged cohort VCF — at 1 000-sample scale this is a serial bottleneck (millions of variant sites, hours of wall time, no resumability on preemption). The correct fix is to annotate per-interval VCFs before `MERGE_GLNEXUS_VCF` and merge annotated results afterward; VEP annotates each record independently by position, so annotate-then-merge is result-equivalent. This requires restructuring the shared `VCF_ANNOTATE_ALL` subworkflow.

**GPU acceleration for DeepVariant** is not configured. At 1 000-sample scale, DeepVariant per-sample calling dominates total wall-clock time; GLnexus joint genotyping is fast by comparison. GPU instances (e.g., `p3.2xlarge` on AWS) accelerate DeepVariant 10–20× over CPU and would be the highest-impact optimization, but require a separate GPU-enabled job queue and container image.

---

## Challenge 1 — Data volume: staging 5 TB of gVCFs

Without Fusion, Nextflow copies every input file into the task work directory first: 24 concurrent GLnexus jobs × 200 GB = ~4.8 TB simultaneous S3 egress, plus large instance disks and slow job start times.

**Seqera Fusion** is a virtual filesystem driver built for Nextflow: it intercepts filesystem calls inside the job container and serves them as HTTP range requests to S3/GCS — no full copy, and only the byte ranges GLnexus actually needs are fetched. Local SSD caches reads; disk is used only for `GLnexus.DB` scratch (~100 GB). Requires a Seqera Platform account and Wave. Without Fusion, set `disk = 250.GB` per GLnexus job for standard staging.

---

## Challenge 2 — GLnexus memory at 1 000-sample scale

| Scope | Approximate RAM |
|---|---|
| Full genome, single job | 300–600 GB |
| Per-chromosome scatter | 30–60 GB |
| Large chromosomes (chr1, chr2) | up to 120 GB |

Per-chromosome scatter keeps each job within `r5.4xlarge` (AWS, 128 GB, ~$0.25/hr Spot) or `n2-highmem-16` (GCP, 128 GB, ~$0.20/hr Spot). Attempt-based escalation handles variance without per-chromosome tuning:

```groovy
memory = { 64.GB * task.attempt }   // 64 → 128 → 192 → 256 GB across attempts 1–4
```

The RAM figures are estimates from published GLnexus benchmarks — exact values depend on cohort size and chromosome length and are not known upfront. In practice, we should run a pilot on a small subset (e.g., 10 samples, all chromosomes), record peak RSS per job, then extrapolate linearly to the target cohort size. The attempt-based escalation is a pragmatic substitute for this profiling step.

**Open file descriptors.** GLnexus holds all input gVCF files open simultaneously. With 1 000 samples this exceeds the OS default of 1 024 file descriptors; the configs raise it to 65 536 via `ulimit -n 65536`. At >10 000 samples two options exist: (1) [Hail](https://hail.is/) (`gvcf_combiner`) uses a distributed merge that avoids the constraint entirely, but requires a Spark cluster (EMR on AWS, Dataproc on GCP) — outside the Batch model; (2) a two-stage hierarchical merge — run GLnexus in batches of ~200 samples to produce intermediate joint VCFs, then merge those — stays within Batch but is not officially supported by GLnexus (which is designed to merge gVCFs, not joint VCFs) and changes the statistical model.

---

## Challenge 3 — Spot/preemptible interruption

A GLnexus job on a large chromosome may run 2–4 hours, making preemption a real risk. The configs retry on preemption exit codes and leave real failures (`exit 1`) to fail fast:

```groovy
// AWS  (104 = transient connection reset)
errorStrategy = { task.exitStatus in [104, 130, 137, 143] ? 'retry' : 'finish' }

// GCP
errorStrategy = { task.exitStatus in [14, 50001] ? 'retry' : 'finish' }
```

AWS escalates to an on-demand queue on the final retry; GCP has no per-task equivalent — a persistent capacity drought requires rerunning with `google.batch.spot = false` (use `-resume` to reuse completed work).

---

## Challenge 4 — Work directory hygiene

At 1 000 samples, intermediate files (per-interval BCFs, VCFs) can reach tens of TB if left uncleaned.

- **`cleanup = true`** — deletes work subdirectories after a successful run (disables `-resume` for that run).
- **Object storage lifecycle rules** — delete work objects after 30 days, covering failed or abandoned runs:

```
AWS: S3 lifecycle rule  — Prefix: sarek-work/  — Expire after 30 days
GCP: GCS lifecycle rule — Prefix: sarek-work/  — Delete after 30 days
```

