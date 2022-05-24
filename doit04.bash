#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 4. Run Bowtie2
# ------------------------------------------------------------------------

echo 1>&2 '# Indexing genome...'

rm -rf ${BOWTIE2}
mkdir -p ${BOWTIE2}

echo 1>&2 '##' ${BOWTIE2}/${REFERENCE_NAME}+phix.fna
cat ${INPUTS}/${REFERENCE_NAME}.fna ${PIPELINE}/inputs/phix.fna > ${BOWTIE2}/${REFERENCE_NAME}+phix.fna
bowtie2-build -q --threads ${THREADS} \
	 ${BOWTIE2}/${REFERENCE_NAME}+phix.fna \
	 ${BOWTIE2}/${REFERENCE_NAME}+phix.fna

echo 1>&2 '# Aligning the reads to the genomes...'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' ${SAMPLES_NAME[$i]}':' ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.bam
    echo 1>&2 '###' bowtie2...
    if [ "$PE" ] ; then
	bowtie2 --seed 1 --threads ${THREADS} \
		 --end-to-end \
		 -x ${BOWTIE2}/${REFERENCE_NAME}+phix.fna \
		 -1 ${FASTP}/trimmed_${SAMPLES_NAME[$i]}_R1.fastq.gz \
		 -2 ${FASTP}/trimmed_${SAMPLES_NAME[$i]}_R2.fastq.gz \
		 -S ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.sam
    else
	bowtie2 --seed 1 --threads ${THREADS} \
		 --end-to-end \
		 -x ${BOWTIE2}/${REFERENCE_NAME}+phix.fna \
		 -U ${FASTP}/trimmed_${SAMPLES_NAME[$i]}_R1.fastq.gz \
		 -S ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.sam
    fi
    echo 1>&2 '###' samtools sort
    samtools sort -@ ${THREADS} ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.sam -o ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.bam
    rm -f ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.sam
    echo 1>&2 '###' samtools index
    samtools index -@ ${THREADS} ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.bam
done

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
