#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs
PROFILES=data/05_profiles
BOWTIE2=data/04_bowtie2
COUNTS=data/06_counts

# ------------------------------------------------------------------------
# Step 8. Compute some high-level stats
# ------------------------------------------------------------------------

echo 1>&2 '# Generating statistics...'

STATS=data/08_stats

rm -rf ${STATS}
mkdir -p ${STATS}/temp

cat ${INPUTS}/annotation.gff inputs/phix.gff \
    | egrep -v '^#' | fgrep -v $'\t'region$'\t' \
			    > ${STATS}/temp/regions.gff

# ------------------------------------------------------------------------

FEATURECOUNTS_ARGS=
FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_biotype"
#FEATURECOUNTS_ARGS+=" -f" # count at feature (exon) level, not the meta-feature (gene) level

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

featureCounts $FEATURECOUNTS_ARGS \
    		  -a ${COUNTS}/annotation.gtf \
    		  -o ${STATS}/raw_counts.txt \
    		  ${BOWTIE2}/aligned_0.bam \
    		  ${BOWTIE2}/aligned_1.bam \
    		  ${BOWTIE2}/aligned_2.bam \
    		  ${BOWTIE2}/aligned_3.bam \
    		  ${BOWTIE2}/aligned_4.bam \
    		  ${BOWTIE2}/aligned_5.bam \


(
    echo gene_biotype$'\t'52551$'\t'52552$'\t'52553$'\t'WT1$'\t'WT2$'\t'WT3
    tail -n+3 ${STATS}/raw_counts.txt \
	| cut -f 1,7-
) > ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
