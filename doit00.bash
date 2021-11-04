#! /bin/bash

. doit-preamble.bash

# ------------------------------------------------------------------------
# Step 0. Set up
# ------------------------------------------------------------------------

rm -rf data

echo 1>&2 '# Initializing data/...'
mkdir -p data/tmp

INPUTS=data/00_inputs
mkdir -p ${INPUTS}

echo 1>&2 '# Making copies of reference genome...'

echo 1>&2 '##' ${INPUTS}/genome.fna "<-" "${REFERENCE_GENOME}"
cp "${REFERENCE_GENOME}" ${INPUTS}/genome.fna
echo 1>&2 '##' ${INPUTS}/annotation.gtf "<-" "${REFERENCE_ANNOTATION_GTF}"
cp "${REFERENCE_ANNOTATION_GTF}" ${INPUTS}/annotation.gtf

if [ "${ADDITIONAL_GENES}" ] ; then
    echo 1>&2 '##' ${INPUTS}/annotation.gtf "+<-" "${ADDITIONAL_GENES}"
    cat "${ADDITIONAL_GENES}" >> ${INPUTS}/annotation.gtf
fi

echo 1>&2 '# Making copies of raw reads...'

for i in $SAMPLES_INDICES ; do
    echo 1>&2 '##' $i':' ${INPUTS}/raw_${i}_R1.fastq.gz "<-" "${SAMPLES_R1[$i]}"
    cp "${SAMPLES_R1[$i]}" ${INPUTS}/raw_${i}_R1.fastq.gz
    if [ "$PE" ] ; then
	echo 1>&2 '##' $i':' ${INPUTS}/raw_${i}_R2.fastq.gz "<-" "${SAMPLES_R2[$i]}"
	cp "${SAMPLES_R2[$i]}" ${INPUTS}/raw_${i}_R2.fastq.gz
    fi
done

if [ "$PACKAGES_FROM" = howto ] ; then
    echo 1>&2 '# Ensuring entries in packages.yaml are downloaded...'
    (
	set -x
	./howto/howto -f packages.yaml -p '*'
    )
fi

