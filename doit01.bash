#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs

# ------------------------------------------------------------------------
# Step 1. FALCO, round 1
# ------------------------------------------------------------------------

echo 1>&2 '# Running FALCO on raw reads'

FALCO1=data/01_falco

rm -rf ${FALCO1}

(
    for i in $SAMPLES_INDICES ; do
	mkdir -p ${FALCO1}/${SAMPLES_NAME[$i]}_R1
	echo falco -q -o ${FALCO1}/${SAMPLES_NAME[$i]}_R1 ${INPUTS}/raw_${i}_R1.fastq.gz
	if [ "$PE" ] ; then
	    mkdir -p ${FALCO1}/${SAMPLES_NAME[$i]}_R2
	    echo falco -q -o ${FALCO1}/${SAMPLES_NAME[$i]}_R2 ${INPUTS}/raw_${i}_R2.fastq.gz
	fi
    done
) | run_parallel

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
