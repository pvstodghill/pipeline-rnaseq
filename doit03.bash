#! /bin/bash

. doit-preamble.bash

FASTP=data/02_fastp

# ------------------------------------------------------------------------
# Step 3. FASTQC, round 2
# ------------------------------------------------------------------------

if [ "$RUN_FASTQC_AFTER_FASTP" ] ; then

    echo 1>&2 '# Running FASTQC on trimmed reads'

    FASTQC2=data/03_fastqc

    rm -rf ${FASTQC2}

    for i in $SAMPLES_INDICES ; do
	echo 1>&2 '##' $i':' ${FASTQC2}/${SAMPLES_NAME[$i]}_R1
	mkdir -p ${FASTQC2}/${SAMPLES_NAME[$i]}_R1
	fastqc -t ${THREADS} \
	       -o ${FASTQC2}/${SAMPLES_NAME[$i]}_R1 \
	       ${FASTP}/trimmed_${i}_R1.fastq.gz
	if [ "$PE" ] ; then
	    echo 1>&2 '##' $i':' ${FASTQC2}/${SAMPLES_NAME[$i]}_R2
	    mkdir -p ${FASTQC2}/${SAMPLES_NAME[$i]}_R2
	    fastqc -t ${THREADS} \
		   -o ${FASTQC2}/${SAMPLES_NAME[$i]}_R2 \
		   ${FASTP}/trimmed_${i}_R2.fastq.gz
	fi
    done

fi

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
