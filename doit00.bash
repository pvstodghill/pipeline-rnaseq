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

echo 1>&2 '##' ${INPUTS}/genome.fna "<-" "${REFERENCES_GENOME}"
cp "${REFERENCES_GENOME}" ${INPUTS}/genome.fna
echo 1>&2 '##' ${INPUTS}/annotation.gff "<-" "${REFERENCES_ANNOTATION_GFF}"
cp "${REFERENCES_ANNOTATION_GFF}" ${INPUTS}/annotation.gff
echo 1>&2 '##' ${INPUTS}/annotation.gtf "<-" "${REFERENCES_ANNOTATION_GTF}"
cp "${REFERENCES_ANNOTATION_GTF}" ${INPUTS}/annotation.gtf

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

