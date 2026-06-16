# Cohort QC report for joint genotyping output — design proposal

## 1. Goal

A joint multi-sample VCF for ~1000 WGS samples needs a dedicated cohort-level QC report
that goes beyond the aggregate statistics already captured in the standard MultiQC report.
The target audience is a bioinformatician or clinical analyst who needs to:

1. Confirm the cohort run succeeded and looks healthy (< 30 seconds)
2. Identify any per-sample outliers that need re-sequencing or flagging (< 5 minutes)
3. Understand what variants are common in this cohort and whether any are clinically or
   epidemiologically noteworthy

---

## 2. What the pipeline already generates

The pipeline already runs `VCF_QC_BCFTOOLS_VCFTOOLS` on the joint VCF
(`workflows/sarek/main.nf:504–512`) and those outputs feed into MultiQC:

- `bcftools stats` — aggregate statistics (one row for the whole cohort: total SNVs,
  indels, Ti/Tv, quality distribution)
- vcftools filter summary and Ti/Tv summaries

**Gap:** `bcftools stats` is called without the `-s -` flag, so the PSC (per-sample counts)
section is not generated. MultiQC's bcftools module would render a 1000-row per-sample
table if that section were present, but today it produces only the single aggregate row.
The plots in Sections C, D, and E below are not generated at all by the current pipeline.

---

## 3. Proposed report format

**Primary:** standalone interactive HTML — a self-contained single file that can be opened
in any browser, attached to a ticket, or stored in object storage alongside the VCF output.

**Phased delivery:**
- **Phase 1 (low effort):** add `--samples -` to `bcftools stats` for the joint VCF. The
  per-sample table (Section B) and some C plots appear in the existing MultiQC report at
  no additional infrastructure cost.
- **Phase 2 (medium effort):** add an Rmarkdown/Quarto step producing a dedicated cohort
  report HTML with interactive plots (SFS, indel length, population comparison) that
  MultiQC does not generate natively.

| Candidate | Interactivity | Portability | Additional effort beyond Phase 1 |
|---|---|---|---|
| **MultiQC** (Phase 1) | Good — sortable tables, hover | Single HTML | Add `-s -` flag only |
| **Rmarkdown / Quarto** (Phase 2) | High — plotly, DT | Single HTML | R container + new process + script |
| **datavzrd** | High — configurable per column | Single HTML | Config file + new process |
| PDF | None | Portable | High — impractical for 1000-row tables |

---

## 4. Report sections

### Section A — Cohort summary (at-a-glance)

Single summary table, one row for the entire cohort. A reviewer should be able to sign off
the run in under 30 seconds from this table.

| Metric | Expected range (30× WGS) |
|---|---|
| Total samples | — |
| Total variant sites | ~5–8M for European ancestry |
| SNVs | ~4.5–6M |
| Indels | ~0.8–1.2M |
| Multi-allelic sites | < 5% of total |
| Ti/Tv (cohort SNVs) | 2.0–2.1 |
| Singleton rate (% of all sites seen in exactly one sample) | 20–40% typical |
| Mean call rate across samples | ≥ 99% healthy; < 95% flagged |
| Samples below 95% call rate | Should be 0 |

---

### Section B — Per-sample QC table (core of the report)

Sortable, filterable table, one row per sample (1 000 rows). Rows exceeding configured
thresholds are highlighted in red.

| Column | Expected range | Why it matters |
|---|---|---|
| Sample ID | — | identity |
| SNV count | 4–5M (30× WGS) | deviation flags coverage issues or sample swap |
| Indel count | 800K–1M | same |
| Ti/Tv ratio | 2.0–2.1 | < 1.8 or > 2.3 flags contamination or alignment errors |
| Het/Hom ratio | 1.5–2.0 | >> 2 → contamination; << 1 → inbreeding or wrong sample |
| Singleton count (% of sample's variants) | varies | elevated rate relative to cohort peers → contamination or private variants |
| Call rate | ≥ 99% | < 95% flags coverage or processing failure |
| Mean depth (from CRAM QC, if available) | 20–30× | cross-reference to explain low call rate |
| PASS filter rate | ≥ 95% | low rate flags systematic quality issue |

**Data source:** `bcftools stats --samples -` PSC section (Phase 1 addition).

---

### Section C — Distribution plots (outlier detection)

**C1. Variant count per sample — bar chart (sorted descending)**
1 000 bars. Outliers visible as bars far above or below the bulk.
Hover shows exact count and sample ID.

**C2. Ti/Tv ratio per sample — dot plot with reference band**
One dot per sample, sorted by Ti/Tv. Reference band at 2.0 ± 0.15.
Dots outside the band colored red.

**C3. Het/Hom ratio — histogram**
Expected bell curve centered ~1.8. Long right tail (contamination) or left tail
(inbreeding / mapping artifact) are immediately visible.

**C4. Call rate — histogram**
Expected: tight distribution near 99%. A second cluster at < 95% indicates a batch of
samples with a systematic processing issue.

**C5. Allele frequency spectrum (site frequency spectrum, SFS)**
X-axis: minor allele frequency bins (0–5%, 5–10%, …, 45–50%).
Y-axis: number of variant sites per bin.
Expected: high singleton/rare-variant peak decaying to a plateau.
A flat spectrum or missing rare-variant peak flags filtering problems or population
stratification artifacts.

**C6. Indel length distribution — histogram**
Counts of insertions vs. deletions by length (−10 to +10 bp).
Expected: slight deletion bias (more −1 than +1). Extreme imbalance flags aligner bias.

**C7. Variant count per chromosome — grouped bar chart**
One bar per chromosome, colored by type (SNV / indel). Sanity check: count should scale
roughly with chromosome length. A deviant chromosome flags alignment or interval errors.

---

### Section D — Cohort-level variant QC

**D1. FILTER status breakdown — stacked bar**
PASS / low-quality / other filters, per variantcaller.
A healthy cohort should have > 95% PASS.

**D2. Variant type breakdown — pie or donut chart**
SNV / insertion / deletion / MNV / complex.
Reference: ~80–85% SNV for WGS. A large MNV fraction can indicate aligner artifacts.

---

### Section E — Most frequent variants: cohort frequency vs. population databases

This section answers: *which variants are common in this cohort, and do any of them have
known clinical or disease relevance?*

**E1. Top-N variants by cohort allele frequency — annotated table**

Table of the N most frequent variants (e.g. top 500 by cohort AF), with columns:

| Column | Source |
|---|---|
| Variant (chr:pos:ref:alt) | Joint VCF |
| Cohort AF | `bcftools query -f '%AF'` or computed from AC/AN in VCF INFO |
| gnomAD AF (non-Finnish European or all-ancestry) | gnomAD via VEP or `bcftools annotate` |
| Fold enrichment (cohort AF ÷ gnomAD AF) | computed |
| ClinVar significance | ClinVar via VEP |
| Gene | VEP consequence annotation |
| Consequence | VEP (missense / synonymous / stop-gained / …) |
| Disease (OMIM) | VEP OMIM plugin (requires free OMIM API key) or HGNC gene summaries |

Rows with ClinVar Pathogenic or Likely Pathogenic are highlighted.
Rows with fold enrichment > 5× are flagged — they indicate variants substantially more
common in this cohort than in the general population, which may reflect shared ancestry,
a founder effect, or a batch-level technical artifact.

**E2. Cohort AF vs. gnomAD AF — scatter plot (log scale)**

X-axis: gnomAD AF (log₁₀). Y-axis: cohort AF (log₁₀). One dot per variant.
Most dots expected near the diagonal (cohort AF ≈ gnomAD AF).
Dots colored by ClinVar significance: grey = benign/VUS, orange = conflicting, red = pathogenic.

- Points above the diagonal: enriched in this cohort vs. general population
- Points below: depleted (selection signal or ancestry effect)
- Pathogenic variants above the diagonal warrant immediate attention

**E3. Clinically significant variants in cohort — summary table**

Filtered view: ClinVar Pathogenic or Likely Pathogenic variants with cohort AF > 0 (present
in at least one sample). Columns: variant, gene, disease, cohort AF, carrier count, ClinVar
star rating, last ClinVar review date.

This is the "clinical signal" table — which disease-associated variants are present and at
what frequency? **Important caveat** (to be included verbatim in the report): this is an
observational summary from a research pipeline. Clinical interpretation, carrier counseling,
and medical decisions require review by a qualified clinical geneticist or certified
diagnostic laboratory.

**Data source availability:**

| Database | GRCh38 | GRCh37/b37 | Access |
|---|---|---|---|
| gnomAD AF | gnomAD v4 | gnomAD v2.1 | VEP plugin or public VCF download |
| ClinVar | ClinVar VCF (monthly) | same | VEP plugin or NCBI FTP download |
| OMIM | via VEP OMIM plugin | same | Free API key required |
| HGNC gene summaries | — | — | Free, no key |

Note: gnomAD v4 is GRCh38-only. For b37/GRCh37 cohorts (e.g. the 1000G test profile),
use gnomAD v2.1 or liftover the cohort VCF before comparison.

---

### Section F — Software versions and run metadata

Standard nf-core versions table:

- GLnexus version
- DeepVariant version
- bcftools version
- VEP version and cache version
- Reference genome build and contig set
- Pipeline version + git commit hash
- Run date

This section is the reproducibility record and should be included even if the rest of the
report is viewed as preliminary.

---

## 5. Implementation sketch

Only the Phase 1 change is trivial enough to describe precisely; Phase 2 is a design
proposal for a future sprint.

**Phase 1 — enable per-sample stats in MultiQC (one config line):**
Add a scoped `ext.args` override in `conf/modules/deepvariant_joint_genotype.config` for a
`BCFTOOLS_STATS_JOINT` alias (so the `-s -` flag does not affect per-sample `BCFTOOLS_STATS`
calls elsewhere in the pipeline), and call it specifically on the joint VCF after
`BAM_JOINT_CALLING_GERMLINE_DEEPVARIANT`. The output feeds into the existing `reports`
channel and is picked up by MultiQC automatically.
Result: Sections B and C1–C4 appear in the existing MultiQC HTML with no new tools.

**Phase 2 — dedicated cohort HTML report:**
New Nextflow process running an Rmarkdown or Quarto script.
Inputs: joint VCF (or per-sample stats TSV from Phase 1), VEP-annotated VCF.
Output: single-file HTML.
Sections C5–C7 (SFS, indel length, per-chromosome), all of Section D, and all of Section E
require this phase.
Required new tools: R + ggplot2/plotly + DT (or Quarto with observable.js); all available
as Docker containers.

**Sections E2 and E3 additionally require:**
- gnomAD / ClinVar annotation in the VCF (already available if VEP is run with those
  plugins — the 1000G test profile runs VEP)
- A script to extract top-N variants by AF and join with the annotation fields
