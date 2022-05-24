#! /bin/bash

REFERENCE_NAME=DC3000
REFERENCE_GENOME=${PIPELINE}/example/NC_004632.fna
REFERENCE_ANNOTATION_GTF=${PIPELINE}/example/NC_004632.gtf
REFERENCE_FEATURE=gene

SAMPLES_NAME[0]=WT1
SAMPLES_TREATMENT[0]=WT
SAMPLES_R1[0]=${PIPELINE}/example/WT1_R1.fastq.gz
SAMPLES_R2[0]=${PIPELINE}/example/WT1_R2.fastq.gz

SAMPLES_NAME[1]=WT2
SAMPLES_TREATMENT[1]=WT
SAMPLES_R1[1]=${PIPELINE}/example/WT2_R1.fastq.gz
SAMPLES_R2[1]=${PIPELINE}/example/WT2_R2.fastq.gz

SAMPLES_NAME[2]=WT3
SAMPLES_TREATMENT[2]=WT
SAMPLES_R1[2]=${PIPELINE}/example/WT3_R1.fastq.gz
SAMPLES_R2[2]=${PIPELINE}/example/WT3_R2.fastq.gz

SAMPLES_NAME[3]=52551
SAMPLES_TREATMENT[3]=5255
SAMPLES_R1[3]=${PIPELINE}/example/52551_R1.fastq.gz
SAMPLES_R2[3]=${PIPELINE}/example/52551_R2.fastq.gz

SAMPLES_NAME[4]=52552
SAMPLES_TREATMENT[4]=5255
SAMPLES_R1[4]=${PIPELINE}/example/52552_R1.fastq.gz
SAMPLES_R2[4]=${PIPELINE}/example/52552_R2.fastq.gz

SAMPLES_NAME[5]=52553
SAMPLES_TREATMENT[5]=5255
SAMPLES_R1[5]=${PIPELINE}/example/52553_R1.fastq.gz
SAMPLES_R2[5]=${PIPELINE}/example/52553_R2.fastq.gz

ADDITIONAL_GENES=
QC_GENES=${PIPELINE}/example/NC_004632_qc.gtf

# ------------------------------------------------------------------------

# Arguments to FASTP
FASTP_ARGS="--trim_front1 1 --trim_front2 1 --adapter_sequence CTGTCTCTTATACACATCT"

# ------------------------------------------------------------------------

if [ -e /programs/docker/bin/docker1 ] ; then
    export HOWTO_DOCKER_CMD=/programs/docker/bin/docker1
fi

# Uncomment to get packages from HOWTO
PACKAGES_FROM=howto

# uncomment to use conda
# PACKAGES_FROM=conda
# PACKAGES_ENV=rnaseq

# Override the default number of threads (nproc --all)
#THREADS=32
