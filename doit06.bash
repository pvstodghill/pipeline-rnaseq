#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs
BOWTIE2=data/04_bowtie2

# ------------------------------------------------------------------------
# Step 6. Make count tables
# ------------------------------------------------------------------------

echo 1>&2 '# Making count tables'

COUNTS=data/06_counts
rm -rf ${COUNTS}
mkdir -p ${COUNTS}

init_FEATURECOUNTS_ARGS

FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level


# FIXME: cp local/annotations.gtf ${COUNTS}/annotation.gtf
fgrep $'\t'gene$'\t' ${INPUTS}/annotation.gtf \
    | ./scripts/sanitize-gtf-for-featureCounts \
	  > ${COUNTS}/annotation.gtf

(
    for i in $SAMPLES_INDICES ; do

	echo featureCounts $FEATURECOUNTS_ARGS \
    	     -a ${COUNTS}/annotation.gtf \
    	     -o ${COUNTS}/counts_$i.txt \
    	     ${BOWTIE2}/aligned_$i.bam

    done
) | run_commands 

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
