#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs

# ------------------------------------------------------------------------
# Step 2. Run FASTP on Illumina reads
# ------------------------------------------------------------------------

echo 1>&2 '# Clean-up Illumina reads...'

FASTP=data/02_fastp

rm -rf ${FASTP}
mkdir -p ${FASTP}

for i in $SAMPLES_INDICES ; do
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FASTP}/trimmed_${i}_R2.fastq.gz
	fastp \
		 --thread ${THREADS} \
		 --trim_front1 1 --trim_front2 1 --adapter_sequence CTGTCTCTTATACACATCT \
		 --json ${FASTP}/${SAMPLES_NAME[$i]}.json \
		 --html ${FASTP}/${SAMPLES_NAME[$i]}.html \
		 --in1 ${INPUTS}/raw_${i}_R1.fastq.gz \
		 --in2 ${INPUTS}/raw_${i}_R2.fastq.gz \
		 --out1 ${FASTP}/trimmed_${i}_R1.fastq.gz \
		 --out2 ${FASTP}/trimmed_${i}_R2.fastq.gz \
		 --unpaired1 ${FASTP}/unpaired_$i.fastq.gz \
		 --unpaired2 ${FASTP}/unpaired_$i.fastq.gz
    else
	echo 1>&2 '##' $i':' ${FASTP}/trimmed_${i}_R1.fastq.gz
	fastp \
		 --thread ${THREADS} \
		 --trim_front1 1 --adapter_sequence CTGTCTCTTATACACATCT \
		 --json ${FASTP}/${SAMPLES_NAME[$i]}.json \
		 --html ${FASTP}/${SAMPLES_NAME[$i]}.html \
		 --in1 ${INPUTS}/raw_${i}_R1.fastq.gz \
		 --out1 ${FASTP}/trimmed_${i}_R1.fastq.gz
    fi
done

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
