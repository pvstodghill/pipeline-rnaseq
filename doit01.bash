#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 1. FALCO, round 1
# ------------------------------------------------------------------------

echo 1>&2 '# Running FALCO on raw reads'

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
) | run_commands_from_stdin

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
