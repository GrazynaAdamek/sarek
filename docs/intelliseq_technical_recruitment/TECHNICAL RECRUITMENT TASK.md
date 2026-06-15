TECHNICAL RECRUITMENT TASK

Extend nf-core/sarek with joint genotyping for a large WGS cohort

Nextflow / nf-core  ·  Bioinformatics  ·  Cloud infrastructure 


Background

A Company uses nf-core/sarek as its primary data processing pipeline for germline whole-genome sequencing (WGS) data. Variant calling is performed using DeepVariant, and the pipeline runs on GCP with input, output and temporary data (scratch / workDir) stored in GCS, and the jobs submitted to GCP Batch.

A research partner has commissioned a population-scale study requiring joint genotyping across a cohort of ~1000 whole-genome sequenced patients. Sarek supports only per-sample variant calling with DeepVariant, producing single-sample VCF files. The joint genotyping step — merging population-level allele evidence across all samples — is not yet implemented for the DeepVariant path in the pipeline.

Your task is to design and implement this capability as a production-ready extension of the sarek workflow. 


 PART 1   Pipeline implementation

Implement a joint genotyping subworkflow for the DeepVariant path in sarek. After the modifications the workflow started with --tools deepvariant —joint_genotype should output one multi-sample, VEP-annotated VCF file including genotypes from the entire input cohort. Follow nf-core standards and re-use components where possible. The accepted solution should run on small test data locally, but ideally, it should also address the challenge described in Part 2 of the task.

Deliverables

   - A Nextflow process module for the joint genotyping tool of your choice, with container, inputs, outputs, and versions.yml

   - A subworkflow that integrates the joint genotyping step across all samples in the cohort

   - Integration into sarek’s existing VARIANT_CALLING workflow, activated by a new flag (e.g. --joint_genotype) when --tools deepvariant is set

   - Integration with sarek's existing ANNOTATE subworkflow to produce a VEP-annotated multi-sample VCF as the final output

   - nf-test tests for the module and subworkflow

   - A small test profile (e.g. 3 samples, single chromosome) that runs locally end-to-end with -profile test,docker

Deliver as a fork of nf-core/sarek or a clearly structured PR-ready branch. Don't commit input data to the repo (unless it is tiny) - add a script that will allow creating it for the review. For test data use 1000Genomes alignment files (https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/), subsetted to a single chromosome for manageable runtimes.


 PART 2   Infrastructure challenge

The pipeline will run at scale in the cloud environment, processing GBs of input data (see below). Considering the following input sizes and infrastructure constraints, implement a solution for running joint genotyping reliably and cost-efficiently at this scale. A cloud-tested solution is not expected, due to costs associated with testing. Instead, prepare a prototype of a nextflow.config for GCP or AWS environment, and an implementation addressing the data-size challenges. Alternatively, a written design document highlighting the challenges and how they could be addressed is also acceptable.

COHORT

~1000 samples
	

CRAM

20–30 GB / sample
	

gVCF

~5 GB / sample

COMPUTE

AWS or GCP Batch
	

STORAGE

S3 / GCS
	

INSTANCES

Spot / preemptible