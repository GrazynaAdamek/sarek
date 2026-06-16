#!/usr/bin/env bash
# prepare_testdata_1000g_chr20.sh — download 1000 Genomes phase-3 chr20 data
# for the joint-genotyping 1000G test profile (conf/test_joint_genotyping_1000g.config).
#
# What it produces (all under tests/data/, which is gitignored):
#   HG00096.chr20.bam / .bai
#   HG00097.chr20.bam / .bai
#   HG00099.chr20.bam / .bai
#   ref.chr20.fasta / .fai / .dict
#   chr20.bed
#
# The 1000G FTP provides pre-split per-chromosome BAMs with companion .md5 files;
# each chr20 BAM is ~20-50 MB.  MD5 checksums are verified after download.
#
# Requirements: wget, samtools (≥1.15), internet access.
#
# Usage:
#   bash scripts/prepare_testdata_1000g_chr20.sh

set -euo pipefail

OUT="tests/data"
mkdir -p "${OUT}"

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"
SAMPLES=(HG00096 HG00097 HG00099)

# ── chr20 BAMs (pre-split by 1000G; verified against companion .md5) ─────────
echo "==> Downloading chr20 BAMs for ${#SAMPLES[@]} samples..."
for S in "${SAMPLES[@]}"; do
    BAM_NAME="${S}.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    BAM_URL="${BASE_URL}/${S}/alignment/${BAM_NAME}"
    OUT_BAM="${OUT}/${S}.chr20.bam"

    if [[ -f "${OUT_BAM}" ]]; then
        echo "  ${S}: ${OUT_BAM} already exists, skipping."
        continue
    fi

    echo "  ${S}: downloading ${BAM_NAME} ..."
    wget -q "${BAM_URL}"          -O "${OUT_BAM}"
    wget -q "${BAM_URL}.bai"      -O "${OUT_BAM}.bai"
    wget -q "${BAM_URL}.md5"      -O "${OUT_BAM}.md5"

    # Verify md5 — the .md5 file contains "<hash>  <original filename>"
    EXPECTED=$(awk '{print $1}' "${OUT_BAM}.md5")
    ACTUAL=$(md5sum "${OUT_BAM}" | awk '{print $1}')
    if [[ "${EXPECTED}" != "${ACTUAL}" ]]; then
        echo "ERROR: MD5 mismatch for ${OUT_BAM} (expected ${EXPECTED}, got ${ACTUAL})" >&2
        exit 1
    fi
    echo "  ${S}: MD5 OK (${EXPECTED})"
    rm "${OUT_BAM}.md5"
done

# ── GRCh37 single-chromosome reference ───────────────────────────────────────
# We stream chr20 directly from the canonical 1000G GRCh37 reference (human_g1k_v37)
# via samtools HTTP range requests — no full genome is written to disk.
#
# Using the same reference build as the 1000G BAMs (human_g1k_v37, b37 coordinate
# system, no "chr" prefix) is important: a mismatched reference (e.g. hg19 with
# "chr" prefixes, or GRCh38) would cause contig name mismatches and silent variant
# calling errors, making the test invalid as a validation.
#
# Trade-off: because only the chr20 region is fetched, there is no pre-computed
# checksum to verify the output against.  If MD5 verification of the reference is
# required, replace this block with a direct download of the Ensembl GRCh37
# per-chromosome FASTA (which ships with a CHECKSUMS file):
#   https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/
#   Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz  (~60 MB, same coordinate system)

REF_OUT="${OUT}/ref.chr20.fasta"
if [[ -f "${REF_OUT}" ]]; then
    echo "==> ${REF_OUT} already exists, skipping reference extraction."
else
    REF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
    echo "==> Extracting chr20 from GRCh37 reference..."
    echo "    (streaming from ${REF_URL} — may take several minutes)"
    samtools faidx "${REF_URL}" "20" > "${REF_OUT}"
    samtools faidx "${REF_OUT}"
    samtools dict "${REF_OUT}" -o "${OUT}/ref.chr20.dict"
fi

# ── Interval BED ──────────────────────────────────────────────────────────────
BED_OUT="${OUT}/chr20.bed"
if [[ ! -f "${BED_OUT}" ]]; then
    CHR_LEN=$(awk '$1=="20"{print $2}' "${REF_OUT}.fai")
    printf "20\t0\t%s\n" "${CHR_LEN}" > "${BED_OUT}"
    echo "==> Created ${BED_OUT} (${CHR_LEN} bp)"
fi

echo ""
echo "==> Done. Test data written to ${OUT}/:"
ls -lh "${OUT}/"
echo ""
echo "Run the pipeline with:"
echo "  NXF_SYNTAX_PARSER=v1 nextflow run main.nf \\"
echo "      -profile test,test_joint_genotyping_1000g,docker \\"
echo "      --outdir results_jg_1000g"
