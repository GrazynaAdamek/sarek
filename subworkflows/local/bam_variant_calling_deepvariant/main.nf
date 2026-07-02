//
// DEEPVARIANT germline calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { DEEPVARIANT_RUNDEEPVARIANT                } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_GVCF } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_DEEPVARIANT_VCF  } from '../../../modules/nf-core/gatk4/mergevcfs/main'

// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow BAM_VARIANT_CALLING_DEEPVARIANT {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    dict          // channel: [optional]  [ meta, dict ]
    fasta         // channel: [mandatory] [ fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    // intervals_name uniquely identifies each (sample, interval) shard so that
    // downstream .join()s on meta can't mis-pair a gVCF with the wrong interval
    // when parallel DeepVariant tasks complete out of order.
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, cram, crai, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals, intervals_name: intervals ? intervals.baseName : null ], cram, crai, intervals ]}

    DEEPVARIANT_RUNDEEPVARIANT(cram_intervals, fasta, fasta_fai, [ [ id:'null' ], [] ], [ [ id:'null' ], [] ])

    // For joint genotyping: per-interval, per-sample gVCFs (before merging)
    gvcf_tbi_intervals = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        .join(DEEPVARIANT_RUNDEEPVARIANT.out.gvcf_index, failOnMismatch: true)
        .join(cram_intervals, failOnMismatch: true)
        .map{ meta, gvcf, tbi, cram, crai, intervals -> [ meta, gvcf, tbi, intervals ] }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_out = DEEPVARIANT_RUNDEEPVARIANT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more gvcf(s) from the same sample
    gvcf_out = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    // intervals_name differs per interval, so it must be stripped from the
    // grouping key first - otherwise every interval gets its own group of size 1
    // and MERGE_DEEPVARIANT_GVCF/VCF never actually merges anything across intervals.
    gvcf_to_merge = gvcf_out.intervals.map{ meta, vcf -> [ groupKey(meta - meta.subMap('intervals_name'), meta.num_intervals), vcf ]}.groupTuple()
    vcf_to_merge = vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta - meta.subMap('intervals_name'), meta.num_intervals), vcf ]}.groupTuple()

    MERGE_DEEPVARIANT_GVCF(gvcf_to_merge, dict)
    MERGE_DEEPVARIANT_VCF(vcf_to_merge, dict)

    // Figuring out if there is one or more tbi(s) from the same sample
    tbi_out = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_index.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Mix intervals and no_intervals channels together
    gvcf = Channel.empty().mix(MERGE_DEEPVARIANT_GVCF.out.vcf, gvcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary fields: num_intervals, intervals_name
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ variantcaller:'deepvariant' ], vcf ] }

    // Mix intervals and no_intervals channels together
    vcf = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.vcf, vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary fields: num_intervals, intervals_name
        .map{ meta, vcf -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ variantcaller:'deepvariant' ], vcf ] }

    tbi = Channel.empty().mix(MERGE_DEEPVARIANT_VCF.out.tbi, tbi_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary fields: num_intervals, intervals_name
        .map{ meta, tbi -> [ meta - meta.subMap('num_intervals', 'intervals_name') + [ variantcaller:'deepvariant' ], tbi ] }

    versions = versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_GVCF.out.versions)
    versions = versions.mix(MERGE_DEEPVARIANT_VCF.out.versions)

    emit:
    gvcf
    gvcf_tbi_intervals
    vcf
    tbi

    versions
}
