//
// JOINT GERMLINE CALLING FOR DEEPVARIANT
//
// Merge per-sample gVCFs across the whole cohort and jointly genotype with GLnexus,
// scattered/gathered over the same intervals used for variant calling.
//

include { GLNEXUS                              } from '../../../modules/local/glnexus/main'
include { BCFTOOLS_VIEW                        } from '../../../modules/nf-core/bcftools/view/main'
include { TABIX_TABIX                          } from '../../../modules/nf-core/tabix/tabix/main'
include { GATK4_MERGEVCFS as MERGE_GLNEXUS_VCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'

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
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when scattered across intervals: merge per-interval VCFs back into one.
    // intervals_name differs per interval, so it must be stripped from the
    // grouping key first - otherwise every interval gets its own group of size 1
    // and MERGE_GLNEXUS_VCF never actually merges anything across intervals.
    vcf_to_merge = vcf_out.intervals
        .map{ meta, vcf -> [ groupKey(meta - meta.subMap('intervals_name'), meta.num_intervals), vcf ] }
        .groupTuple()

    MERGE_GLNEXUS_VCF(vcf_to_merge, dict)

    // Single-interval / no-intervals case needs its own index
    TABIX_TABIX(vcf_out.no_intervals)

    // Mix intervals and no_intervals channels together
    // Rework meta for variantscalled.csv and annotation tools
    genotype_vcf = Channel.empty().mix(MERGE_GLNEXUS_VCF.out.vcf, vcf_out.no_intervals)
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'deepvariant' ], vcf ] }

    genotype_index = Channel.empty().mix(MERGE_GLNEXUS_VCF.out.tbi, TABIX_TABIX.out.tbi)
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ id:'joint_variant_calling', patient:'all_samples', variantcaller:'deepvariant' ], tbi ] }

    versions = versions.mix(GLNEXUS.out.versions)
    versions = versions.mix(BCFTOOLS_VIEW.out.versions)
    versions = versions.mix(MERGE_GLNEXUS_VCF.out.versions)
    versions = versions.mix(TABIX_TABIX.out.versions)

    emit:
    genotype_vcf    // channel: [ val(meta), [ vcf ] ]
    genotype_index  // channel: [ val(meta), [ tbi ] ]

    versions        // channel: [ versions.yml ]
}
