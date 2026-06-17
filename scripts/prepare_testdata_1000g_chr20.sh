#!/usr/bin/env bash
# prepare_testdata_1000g_chr20.sh — download 1000 Genomes phase-3 chr20 data
# for the joint-genotyping 1000G test profiles (conf/test_joint_genotyping_1000g.config,
# conf/test_joint_genotyping_1000g_intervals.config).
#
# What it produces (all under tests/data/, which is gitignored):
#   HG00096.chr20.bam / .bai
#   HG00097.chr20.bam / .bai
#   HG00099.chr20.bam / .bai
#   ref.chr20.fasta / .fai / .dict
#   chr20.bed                          (full chr20, single interval)
#   chr20.multi_intervals.bed          (3 equal full-chr20 chunks, unused by current profiles)
#   chr20_subset.multi_intervals.bed   (3 x 500 kbp windows, used by the intervals profile)
#
# The 1000G FTP provides pre-split per-chromosome BAMs; each chr20 BAM is ~20-50 MB.
# No checksum verification is performed on download.
#
# Requirements: docker, internet access.
#
# Usage:
#   bash scripts/prepare_testdata_1000g_chr20.sh [--vep DIR]
#
# Options:
#   --vep DIR   Also download the VEP cache (homo_sapiens GRCh37 v115, ~15 GB) to DIR.
#               Pass the same DIR as --vep_cache when running the pipeline.

set -euo pipefail

# ── Argument parsing ───────────────────────────────────────────────────────────
VEP_CACHE_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vep)
            [[ $# -ge 2 ]] || { echo "ERROR: --vep requires a directory argument" >&2; exit 1; }
            VEP_CACHE_DIR="$2"; shift 2 ;;
        --vep=*)
            VEP_CACHE_DIR="${1#--vep=}"; shift ;;
        *)
            echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -f /.dockerenv ]]; then
    docker run --rm \
        -v "$(pwd):/work" \
        -w /work \
        staphb/samtools:1.21 \
        bash -c "apt-get update -qq && apt-get install -y -qq curl && bash scripts/prepare_testdata_1000g_chr20.sh"

    if [[ -n "$VEP_CACHE_DIR" ]]; then
        mkdir -p "$VEP_CACHE_DIR"
        echo "==> Downloading VEP cache (homo_sapiens GRCh37 v110) to ${VEP_CACHE_DIR} (~15 GB, may take a while)..."
        docker run --rm \
            -v "${VEP_CACHE_DIR}:/cache" \
            community.wave.seqera.io/library/ensembl-vep_perl-math-cdf:1e13f65f931a6954 \
            vep_install \
                --AUTO cf \
                --SPECIES homo_sapiens \
                --ASSEMBLY GRCh37 \
                --CACHE_VERSION 115 \
                --CACHEDIR /cache \
                --NO_UPDATE \
                --NO_HTSLIB
        echo "==> VEP cache written to ${VEP_CACHE_DIR}"
    fi
    exit 0
fi

OUT="tests/data"
mkdir -p "${OUT}"

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"

declare -A BAM_URLS=(
    [HG00096]="${BASE_URL}/HG00096/alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    [HG00097]="${BASE_URL}/HG00097/alignment/HG00097.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"
    [HG00099]="${BASE_URL}/HG00099/alignment/HG00099.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"
)

# ── chr20 BAMs ────────────────────────────────────────────────────────────────
echo "==> Downloading chr20 BAMs for ${#BAM_URLS[@]} samples..."
for S in "${!BAM_URLS[@]}"; do
    OUT_BAM="${OUT}/${S}.chr20.bam"

    if [[ -f "${OUT_BAM}" && -f "${OUT_BAM}.bai" ]]; then
        echo "  ${S}: ${OUT_BAM} already exists, skipping."
        continue
    fi

    BAM_URL="${BAM_URLS[$S]}"
    echo "  ${S}: downloading $(basename "${BAM_URL}") ..."
    curl -fsSL "${BAM_URL}"      -o "${OUT_BAM}"
    curl -fsSL "${BAM_URL}.bai"  -o "${OUT_BAM}.bai"
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
if [[ -f "${REF_OUT}" && -f "${REF_OUT}.fai" ]]; then
    echo "==> ${REF_OUT} already exists, skipping reference download."
else
    # Ensembl GRCh37 per-chromosome FASTA (~60 MB); same b37 coordinate system
    # (no "chr" prefix) as the 1000G BAMs.
    REF_URL="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz"
    echo "==> Downloading chr20 reference from Ensembl GRCh37 (~60 MB)..."
    curl -fsSL "${REF_URL}" | gunzip > "${REF_OUT}"
    samtools faidx "${REF_OUT}"
    samtools dict "${REF_OUT}" -o "${OUT}/ref.chr20.dict"
fi

# ── Interval BEDs ─────────────────────────────────────────────────────────────
CHR_LEN=$(awk '$1=="20"{print $2}' "${REF_OUT}.fai")

BED_OUT="${OUT}/chr20.bed"
if [[ ! -f "${BED_OUT}" ]]; then
    printf "20\t0\t%s\n" "${CHR_LEN}" > "${BED_OUT}"
    echo "==> Created ${BED_OUT} (${CHR_LEN} bp)"
fi

BED_MULTI_OUT="${OUT}/chr20.multi_intervals.bed"
if [[ ! -f "${BED_MULTI_OUT}" ]]; then
    CHUNK=$(( CHR_LEN / 3 ))
    printf "20\t1\t%s\n"                         "${CHUNK}"        >  "${BED_MULTI_OUT}"
    printf "20\t%s\t%s\n"  "$((CHUNK + 1))"      "$((CHUNK * 2))" >> "${BED_MULTI_OUT}"
    printf "20\t%s\t%s\n"  "$((CHUNK * 2 + 1))"  "${CHR_LEN}"     >> "${BED_MULTI_OUT}"
    echo "==> Created ${BED_MULTI_OUT} (3 intervals of ~$((CHR_LEN / 3 / 1000000)) Mbp each)"
fi

BED_SUBSET_MULTI_OUT="${OUT}/chr20_subset.multi_intervals.bed"
if [[ ! -f "${BED_SUBSET_MULTI_OUT}" ]]; then
    printf "20\t10000000\t10500000\n" >  "${BED_SUBSET_MULTI_OUT}"
    printf "20\t30000000\t30500000\n" >> "${BED_SUBSET_MULTI_OUT}"
    printf "20\t55000000\t55500000\n" >> "${BED_SUBSET_MULTI_OUT}"
    echo "==> Created ${BED_SUBSET_MULTI_OUT} (3 × 500 kbp regions spread across chr20)"
fi

echo ""
echo "==> Done. Test data written to ${OUT}/:"
ls -lh "${OUT}/"
echo ""
echo "Joint genotyping module only (no VEP):"
echo "  NXF_SYNTAX_PARSER=v1 nextflow run main.nf \\"
echo "      -profile test,test_joint_genotyping_1000g,docker \\"
echo "      --outdir results_jg_1000g"
echo ""
echo "Full end-to-end incl. VEP (scatter/gather over 1.5 Mbp of chr20; requires a VEP cache,"
echo "download one with: bash scripts/prepare_testdata_1000g_chr20.sh --vep /path/to/vep_cache):"
echo "  NXF_SYNTAX_PARSER=v1 nextflow run main.nf \\"
echo "      -profile test,test_joint_genotyping_1000g_intervals,docker \\"
echo "      --vep_cache /path/to/vep_cache \\"
echo "      --outdir results_jg_1000g_intervals"
