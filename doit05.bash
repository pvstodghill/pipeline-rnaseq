#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 5. Determine the orientation of the raw reads
# ------------------------------------------------------------------------

if [ "${ORIENTATION}" ] ; then
    echo 1>&2 "# Orientation set mannual:" ${ORIENTATION}
    exit 0
fi

rm -rf ${STRAND}
mkdir -p ${STRAND}

# ------------------------------------------------------------------------

cat ${INPUTS}/annotation.gtf \
    | fgrep $'\t'gene$'\t' \
    | ${PIPELINE}/scripts/sanitize-gtf-for-featureCounts \
		 > ${STRAND}/targets.gtf

# ------------------------------------------------------------------------

init_FEATURECOUNTS_ARGS

FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level

# ------------------------------------------------------------------------

echo 1>&2 '# Running featureCounts: forward'

featureCounts $FEATURECOUNTS_ARGS -s 1 \
    	      -a ${STRAND}/targets.gtf \
    	      -o ${STRAND}/forward.txt \
	      ${BOWTIE2}/aligned_${SAMPLES_NAME[0]}.bam

echo 1>&2 '# Running featureCounts: reverse'

featureCounts $FEATURECOUNTS_ARGS -s 2 \
    	      -a ${STRAND}/targets.gtf \
    	      -o ${STRAND}/reverse.txt \
	      ${BOWTIE2}/aligned_${SAMPLES_NAME[0]}.bam

# ------------------------------------------------------------------------

echo 1>&2 '# Analyzing results'

${PIPELINE}/scripts/rnaseq-strand-analysis.pl \
	   -o ${STRAND}/results.sh \
	   ${STRAND}/forward.txt \
	   ${STRAND}/reverse.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'

