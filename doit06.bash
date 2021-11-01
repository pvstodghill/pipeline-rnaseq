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

FEATURECOUNTS_ARGS=
FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level

FEATURECOUNTS_ARGS+=" -O" # Assign reads to all their overlapping features

#FEATURECOUNTS_ARGS+=" -M" # all multi-mapping reads reported alignments will be counted
##FEATURECOUNTS_ARGS+=" --fraction" # Assign fractional counts to features
#FEATURECOUNTS_ARGS+=" --primary" # Count primary alignments only !(0x100)

FEATURECOUNTS_ARGS+=" -s 1" # stranded
#FEATURECOUNTS_ARGS+=" -s 2" # reverse-stranded


FEATURECOUNTS_ARGS+=" -p" # fragments (or pairs) will be counted instead of reads (<v2.0.2)
#FEATURECOUNTS_ARGS+=" -p" # input data contains paired-end reads. (>=v2.0.2)
#FEATURECOUNTS_ARGS+=" --countReadPairs" # Count read pairs (fragments) instead of reads

FEATURECOUNTS_ARGS+=" -B" # Only count read pairs that have both ends aligned.
FEATURECOUNTS_ARGS+=" -P" # Check validity of paired-end distance
FEATURECOUNTS_ARGS+=" -C" # Only count concordant reads

FEATURECOUNTS_ARGS+=" -T ${THREADS}"

FEATURECOUNTS_ARGS+=" "
FEATURECOUNTS_ARGS+=" "


# FIXME: cp local/annotations.gtf ${COUNTS}/annotation.gtf
fgrep $'\t'gene$'\t' ${INPUTS}/annotation.gtf \
      | ./scripts/sanitize-gtf-for-featureCounts \
	    > ${COUNTS}/annotation.gtf

for i in $SAMPLES_INDICES ; do

    featureCounts $FEATURECOUNTS_ARGS \
    		  -a ${COUNTS}/annotation.gtf \
    		  -o ${COUNTS}/counts_$i.txt \
    		  ${BOWTIE2}/aligned_$i.bam

done


# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
