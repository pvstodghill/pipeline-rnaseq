#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 4. Run Bowtie2
# ------------------------------------------------------------------------

echo 1>&2 '# Indexing genome...'

rm -rf ${BOWTIE2}
mkdir -p ${BOWTIE2}

echo 1>&2 '##' ${BOWTIE2}/genome+phix.fna
cat ${INPUTS}/genome.fna ${PIPELINE}/inputs/phix.fna > ${BOWTIE2}/genome+phix.fna
bowtie2-build -q --threads ${THREADS} \
	 ${BOWTIE2}/genome+phix.fna \
	 ${BOWTIE2}/genome+phix.fna

echo 1>&2 '# Aligning the reads to the genomes...'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${BOWTIE2}/aligned_$i.bam
    echo 1>&2 '###' bowtie2...
    if [ "$PE" ] ; then
	bowtie2 --seed 1 --threads ${THREADS} \
		 --end-to-end \
		 -x ${BOWTIE2}/genome+phix.fna \
		 -1 ${FASTP}/trimmed_${i}_R1.fastq.gz \
		 -2 ${FASTP}/trimmed_${i}_R2.fastq.gz \
		 -S ${BOWTIE2}/aligned_$i.sam
    else
	bowtie2 --seed 1 --threads ${THREADS} \
		 --end-to-end \
		 -x ${BOWTIE2}/genome+phix.fna \
		 -U ${FASTP}/trimmed_${i}_R1.fastq.gz \
		 -S ${BOWTIE2}/aligned_$i.sam
    fi
    echo 1>&2 '###' samtools sort
    samtools sort -@ ${THREADS} ${BOWTIE2}/aligned_$i.sam -o ${BOWTIE2}/aligned_$i.bam
    rm -f ${BOWTIE2}/aligned_$i.sam
    echo 1>&2 '###' samtools index
    samtools index -@ ${THREADS} ${BOWTIE2}/aligned_$i.bam
done

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
