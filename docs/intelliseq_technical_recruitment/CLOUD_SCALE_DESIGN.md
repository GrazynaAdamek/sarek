# Cloud-scale joint genotyping: design document

## 1. Scale facts

| Resource | Per sample | 1 000 samples |
|---|---|---|
| CRAM | 20–30 GB | 20–30 TB |
| gVCF | ~5 GB | ~5 TB |
| gVCF per interval (≈ per chromosome) | ~200 MB | ~200 GB per GLnexus job |

Per-sample CRAM (20–30 GB) and gVCF (~5 GB) figures assume higher-coverage WGS; lighter
coverage shifts the totals down proportionally but not the scaling argument below.

The joint genotyping pipeline from Part 1 already scatters GLnexus across the
variant-calling intervals (one job per interval — roughly per chromosome with the default
WGS interval set, though the exact count depends on the intervals BED). That scatter is the
foundation that makes cloud-scale tractable: instead of one job touching all 5 TB of gVCFs,
each parallel job touches only ~200 GB.

---

## 2. Challenge 1 — Data volume: staging 5 TB of gVCFs

**Problem.** By default, Nextflow copies every input file from S3/GCS into the task work
directory before a process starts. For a per-chromosome GLnexus job, that means:

```
1 000 samples × ~200 MB gVCF/chromosome = ~200 GB copied per job
24 concurrent jobs                       = ~4.8 TB simultaneous S3 egress
```

That is slow (minutes of copy before computation starts), expensive (S3 egress + large
instance disks), and brittle (a copy failure means restarting the job from scratch).

**Solution: Seqera Fusion filesystem.**

Fusion mounts S3 or GCS as a virtual POSIX filesystem backed by HTTP range requests.
From GLnexus's perspective the gVCF files are local paths — but reads are served directly
from object storage, with no prior copy. GLnexus reads only the positions it needs (it
does not scan whole gVCF files sequentially), so the effective I/O is much smaller than
the nominal 200 GB per job.

```
Nextflow sees:  /fusion/s3/bucket/sample1.chr20.g.vcf.gz  (virtual local path)
Actual I/O:     HTTP GET with Range header → S3 → GLnexus process
Disk used:      ~100 GB for GLnexus.DB scratch only (not for input files)
```

Fusion is configured identically for both AWS and GCP:

```groovy
fusion { enabled = true; exportStorageCredentials = true }
wave   { enabled = true; endpoint = 'https://wave.seqera.io' }
```

Note that Fusion also caches reads on local disk, so it benefits from fast NVMe/SSD scratch.
The `disk = 100.GB` setting in the configs sizes the scratch volume but does not pin its
type — on compute environments where the instance disk is slow, Fusion throughput suffers,
so prefer instance types with local SSD for the GLnexus jobs.

It requires a Seqera Platform account (free tier, subject to Seqera's limits) for the Wave
service that builds Fusion-enabled containers on the fly. The same container images the
pipeline already uses (Docker/Singularity) are supported — no Dockerfile changes needed.

**Fallback (without Fusion).** Standard Nextflow staging works but requires large instance
disks (`disk = 250.GB` per GLnexus job) and adds significant wall-clock time before
computation starts. Use this if running without a Seqera Platform account.

---

## 3. Challenge 2 — GLnexus memory at 1 000-sample scale

GLnexus memory usage scales with the number of variant sites in the region, not file count.
Per-chromosome figures from the GLnexus documentation and benchmarks:

| Scope | Samples | Approximate RAM |
|---|---|---|
| Full genome (naive, single job) | 1 000 | 300–600 GB |
| Per-chromosome job (scatter) | 1 000 | 30–60 GB |
| Large chromosomes (chr1, chr2) | 1 000 | up to 120 GB |

The per-chromosome scatter keeps each job within the range of `r5.4xlarge` (AWS, 128 GB RAM,
~$0.25/hr Spot) or `n2-highmem-16` (GCP, 128 GB RAM, ~$0.20/hr Spot).

Memory is configured with attempt-based escalation so OOMs on larger chromosomes
automatically retry with more capacity:

```groovy
memory = { 64.GB * task.attempt }   // 64 → 128 → 192 → 256 GB across attempts 1–4 (maxRetries = 3)
```

The largest chromosomes (~120 GB) are covered by the second attempt (128 GB); the remaining
escalation is headroom. No manual tuning per chromosome is required; the retry loop covers
the variance.

The `ulimit -n 65536` guard already added to the module (Part 1) prevents the default
open-file limit (1 024) from being exhausted at 1 000 samples.

---

## 4. Challenge 3 — Spot/preemptible interruption

Spot instances can be reclaimed at any time. A GLnexus job on a large chromosome may run
for 2–4 hours, making preemption a real operational risk.

**Exit codes that signal preemption (not bugs):**

| Cloud | Exit code | Meaning |
|---|---|---|
| AWS | 137 | SIGKILL (spot termination) |
| AWS | 130 | SIGINT |
| AWS | 143 | SIGTERM |
| GCP | 14 | PREEMPTED |
| GCP | 50001 | Batch Spot reclaim (observed; not in the documented exit-code table) |

The configs retry on these preemption/transient codes rather than failing the pipeline.
Note that generic exit code `1` is deliberately **not** retried — it signals a real failure
(malformed gVCF, a GLnexus bug, an OOM surfacing as 1) that retrying would only delay, at the
cost of hours per attempt:

```groovy
// AWS  (104 = transient connection reset)
errorStrategy = { task.exitStatus in [104, 130, 137, 143] ? 'retry' : 'finish' }

// GCP
errorStrategy = { task.exitStatus in [14, 50001] ? 'retry' : 'finish' }
```

**Safe restarts.** The `rm -rf GLnexus.DB` guard added to the module in Part 1 removes the
stale scratch database that `glnexus_cli` would otherwise refuse to overwrite on retry.
Without it, every preempted GLnexus job would fail permanently on re-run.

**On-demand escalation (AWS).** If a spot capacity zone is consistently preempting, all
retries against Spot will fail. AWS escapes this automatically by switching to an on-demand
job queue on the final retry:

```groovy
// AWS: two named queues, escalate to on-demand on the final attempt
queue = { task.attempt > 2 ? params.aws_ondemand_queue : params.aws_spot_queue }
```

**GCP has no equivalent per-task escape.** On GCP Batch, Spot is a *global* setting
(`google.batch.spot`), and the executor exposes no per-task spot directive — so a single
retry cannot be promoted to on-demand the way the AWS queues allow. A persistent capacity
drought is therefore handled only by the Spot retries (`maxRetries = 3`); escaping it
entirely means rerunning with `google.batch.spot = false` (use `-resume` so completed work
is reused and only the stuck jobs rerun on-demand). This asymmetry between the two clouds is
inherent to GCP Batch, not an oversight in the config.

---

## 5. Challenge 4 — Work directory hygiene

The S3/GCS work directory accumulates intermediate files: per-interval gVCF copies (if not
using Fusion), per-interval BCFs, per-interval VCFs, and the merged cohort VCF. At 1 000
samples this can reach tens of TB if not cleaned up.

Two complementary mechanisms:

**`cleanup = true`** (in both configs) — Nextflow deletes the work subdirectories of
successfully completed processes at the end of a successful pipeline run. This is the primary
cleanup path. Note this only fires on *success*, and once a successful run is cleaned up it
can no longer be `-resume`d — which is fine, because there is nothing left to resume. The
lifecycle rules below are what cover *failed* runs, where `-resume` still matters.

**Object storage lifecycle rules** — handle runs that fail, are killed, or are abandoned:

```
# AWS: S3 lifecycle rule (set via console or Terraform)
Prefix:     sarek-work/
Action:     Expire (delete) objects after 30 days

# GCP: GCS lifecycle rule
Prefix:     sarek-work/
Action:     Delete objects older than 30 days
```

The 30-day window is long enough to resume a failed run (Nextflow's `-resume` reads the
`.nextflow/cache` database, which references the existing work directories) but short enough
to avoid unbounded storage growth.

---

## 6. How to run

### AWS Batch

```bash
# Prerequisites: two AWS Batch job queues (spot + on-demand), Seqera Platform account

# For the 1000G 3-sample chr20 test:
bash scripts/prepare_testdata_1000g_chr20.sh

NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test_joint_genotyping_1000g,cloud_aws,docker \
    -work-dir s3://your-bucket/sarek-work \
    --outdir   s3://your-bucket/sarek-results \
    --aws_spot_queue     sarek-spot-queue \
    --aws_ondemand_queue sarek-ondemand-queue \
    --aws_region         eu-west-1
```

### GCP Batch

```bash
# Prerequisites: GCP project with Batch API enabled, Seqera Platform account

bash scripts/prepare_testdata_1000g_chr20.sh

NXF_SYNTAX_PARSER=v1 nextflow run main.nf \
    -profile test_joint_genotyping_1000g,cloud_gcp,docker \
    -work-dir gs://your-bucket/sarek-work \
    --outdir   gs://your-bucket/sarek-results \
    --gcp_project your-gcp-project \
    --gcp_region  europe-west4
```

The cloud configs are additive: they set the executor, Fusion, and per-process resource
overrides. The test profile (`test_joint_genotyping_1000g`) supplies input, reference, and
VEP parameters. Both are required together.

---

## 7. Out of scope

The following were identified but not implemented (Part 2 accepts a design over an
implementation for these):

- **Scatter VEP annotation per interval** — annotating per-interval VCFs before the merge
  step would give parallelism and resumability for the VEP job, which is currently a single
  serial job on the merged cohort VCF. Design already in
  `JOINT_GENOTYPING_IMPLEMENTATION.md §9`.
- **GPU acceleration for DeepVariant** — `DEEPVARIANT_RUNDEEPVARIANT` supports GPU containers
  (the Wave/Seqera ecosystem can build GPU-enabled images). With 1 000 samples this is the
  dominant wall-clock cost upstream of joint genotyping. Requires a GPU-enabled compute
  environment and `accelerator` config in `conf/base.config` — the conditional pattern is
  already there for GPU-labelled processes (`base.config:31-33`, `withLabel: process_gpu`;
  `PARABRICKS_FQ2BAM` at `base.config:78` uses the same idea).
- **Terraform / CDK for compute environment provisioning** — the AWS Batch compute
  environments (spot + on-demand queues, IAM roles, S3 bucket, lifecycle rules) and the
  equivalent GCP Batch setup are infrastructure-as-code work outside the pipeline scope.
