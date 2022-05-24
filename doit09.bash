#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 9. Compute some high-level stats
# ------------------------------------------------------------------------

echo 1>&2 '# Generating statistics...'

rm -rf ${STATS}
mkdir -p ${STATS}/temp


cp ${COUNTS}/annotation.gtf ${STATS}/
if [ "${QC_GENES}" ] ; then
    cat "${QC_GENES}" >> ${STATS}/annotation.gtf
fi

# ------------------------------------------------------------------------

init_FEATURECOUNTS_ARGS

FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_biotype"
#FEATURECOUNTS_ARGS+=" -f" # count at feature (exon) level, not the meta-feature (gene) level
FEATURECOUNTS_ARGS+=" -T ${THREADS}"

# FIXME: parameterize
ARGS=
for i in $SAMPLES_INDICES ; do
    ARGS+=" ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.bam"
done

featureCounts $FEATURECOUNTS_ARGS \
    	      -a ${STATS}/annotation.gtf \
    	      -o ${STATS}/raw_counts.txt \
	      $ARGS

(
    echo -n gene_biotype
    for i in $SAMPLES_INDICES ; do
	echo -n $'\t'${SAMPLES_NAME[$i]}
    done
    echo ''
    tail -n+3 ${STATS}/raw_counts.txt \
	| cut -f 1,7-
) > ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
