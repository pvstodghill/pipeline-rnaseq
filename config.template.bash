#! /bin/bash

if [ -e /programs/docker/bin/docker1 ] ; then
    export HOWTO_DOCKER_CMD=/programs/docker/bin/docker1
fi

REFERENCE_NAME=DC3000
REFERENCE_GENOME=$(ls inputs/ncbi_dataset/data/GCF_000007805.1/*.fna)
REFERENCE_ANNOTATION_GFF=inputs/regions-from-iaa-rnaseq.gff
REFERENCE_ANNOTATION_GTF=inputs/regions-from-iaa-rnaseq.gtf
REFERENCE_FEATURE=

ILLUMINA=/home/ps27/midden/2020/10-30-5255-rnaseq-brc-downloads/

SAMPLES_NAME[0]=WT1
SAMPLES_TREATMENT[0]=WT
SAMPLES_R1[0]=$(ls ${ILLUMINA}/*_WT_1_*_R1.fastq.gz)
SAMPLES_R2[0]=$(ls ${ILLUMINA}/*_WT_1_*_R2.fastq.gz)

SAMPLES_NAME[1]=WT2
SAMPLES_TREATMENT[1]=WT
SAMPLES_R1[1]=$(ls ${ILLUMINA}/*_WT_2_*_R1.fastq.gz)
SAMPLES_R2[1]=$(ls ${ILLUMINA}/*_WT_2_*_R2.fastq.gz)

SAMPLES_NAME[2]=WT3
SAMPLES_TREATMENT[2]=WT
SAMPLES_R1[2]=$(ls ${ILLUMINA}/*_WT_3_*_R1.fastq.gz)
SAMPLES_R2[2]=$(ls ${ILLUMINA}/*_WT_3_*_R2.fastq.gz)

SAMPLES_NAME[3]=52551
SAMPLES_TREATMENT[3]=5255
SAMPLES_R1[3]=$(ls ${ILLUMINA}/*_5255_1_*_R1.fastq.gz)
SAMPLES_R2[3]=$(ls ${ILLUMINA}/*_5255_1_*_R2.fastq.gz)

SAMPLES_NAME[4]=52552
SAMPLES_TREATMENT[4]=5255
SAMPLES_R1[4]=$(ls ${ILLUMINA}/*_5255_2_*_R1.fastq.gz)
SAMPLES_R2[4]=$(ls ${ILLUMINA}/*_5255_2_*_R2.fastq.gz)

SAMPLES_NAME[5]=52553
SAMPLES_TREATMENT[5]=5255
SAMPLES_R1[5]=$(ls ${ILLUMINA}/*_5255_3_*_R1.fastq.gz)
SAMPLES_R2[5]=$(ls ${ILLUMINA}/*_5255_3_*_R2.fastq.gz)
