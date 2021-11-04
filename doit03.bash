#! /bin/bash

. doit-preamble.bash

FASTP=data/02_fastp

# ------------------------------------------------------------------------
# Step 3. FALCO, round 2
# ------------------------------------------------------------------------

if [ "$RUN_FALCO_AFTER_FASTP" ] ; then

    echo 1>&2 '# Running FALCO on trimmed reads'

    FALCO2=data/03_falco

    rm -rf ${FALCO2}

    (
	for i in $SAMPLES_INDICES ; do
	    mkdir -p ${FALCO2}/${SAMPLES_NAME[$i]}_R1
	    echo falco -q -o ${FALCO2}/${SAMPLES_NAME[$i]}_R1 ${FASTP}/trimmed_${i}_R1.fastq.gz
	    if [ "$PE" ] ; then
		mkdir -p ${FALCO2}/${SAMPLES_NAME[$i]}_R2
		echo falco -q -o ${FALCO2}/${SAMPLES_NAME[$i]}_R2 ${FASTP}/trimmed_${i}_R2.fastq.gz
	    fi
	done
    ) | run_parallel

fi

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
