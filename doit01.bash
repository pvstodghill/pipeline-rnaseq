#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs

# ------------------------------------------------------------------------
# Step 1. FALCO, round 1
# ------------------------------------------------------------------------

echo 1>&2 '# Running FALCO on raw reads'

FALCO1=data/01_falco

rm -rf ${FALCO1}

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${FALCO1}/${SAMPLES_NAME[$i]}_R1
    mkdir -p ${FALCO1}/${SAMPLES_NAME[$i]}_R1
    falco -t ${THREADS} \
	   -o ${FALCO1}/${SAMPLES_NAME[$i]}_R1 \
	   ${INPUTS}/raw_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${FALCO1}/${SAMPLES_NAME[$i]}_R2
	mkdir -p ${FALCO1}/${SAMPLES_NAME[$i]}_R2
	falco -t ${THREADS} \
	       -o ${FALCO1}/${SAMPLES_NAME[$i]}_R2 \
	       ${INPUTS}/raw_${i}_R2.fastq.gz
    fi
done

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
