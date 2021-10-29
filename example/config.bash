#! /bin/bash

if [ -e /programs/docker/bin/docker1 ] ; then
    export HOWTO_DOCKER_CMD=/programs/docker/bin/docker1
fi

REFERENCES_NAME=DC3000
REFERENCES_GENOME=example/NC_004632.fna
REFERENCES_ANNOTATION=example/NC_004632.gff
REFERENCES_FEATURE=gene

SAMPLES_NAME[0]=52551
SAMPLES_R1[0]=example/52551_R1.fastq.gz
SAMPLES_R2[0]=example/52551_R2.fastq.gz

SAMPLES_NAME[1]=52552
SAMPLES_R1[1]=example/52552_R1.fastq.gz
SAMPLES_R2[1]=example/52552_R2.fastq.gz

SAMPLES_NAME[2]=52553
SAMPLES_R1[2]=example/52553_R1.fastq.gz
SAMPLES_R2[2]=example/52553_R2.fastq.gz

SAMPLES_NAME[3]=WT1
SAMPLES_R1[3]=example/WT1_R1.fastq.gz
SAMPLES_R2[3]=example/WT1_R2.fastq.gz

SAMPLES_NAME[4]=WT2
SAMPLES_R1[4]=example/WT2_R1.fastq.gz
SAMPLES_R2[4]=example/WT2_R2.fastq.gz

SAMPLES_NAME[5]=WT3
SAMPLES_R1[5]=example/WT3_R1.fastq.gz
SAMPLES_R2[5]=example/WT3_R2.fastq.gz

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

# fixme - document
#RUN_FASTQC_AFTER_FASTP=yes
