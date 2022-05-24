#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 6. Make count tables
# ------------------------------------------------------------------------

echo 1>&2 '# Making count tables'

rm -rf ${COUNTS}
mkdir -p ${COUNTS}

# ------------------------------------------------------------------------

init_FEATURECOUNTS_ARGS

FEATURECOUNTS_ARGS+=" -t gene"
FEATURECOUNTS_ARGS+=" -g gene_id"
FEATURECOUNTS_ARGS+=" -f" # count at feature level

if [ -z "$ORIENTATION" ] ; then
    . ${STRAND}/results.sh
fi
case X"$ORIENTATION"X in
    XforwardX)
	FEATURECOUNTS_ARGS+=" -s 1" # stranded
	;;
    XreverseX)
	FEATURECOUNTS_ARGS+=" -s 2" # reverse-stranded
	;;
    X*X) echo 1>&2 cannot happen ; exit 1
esac

# ------------------------------------------------------------------------

# FIXME: cp local/annotations.gtf ${COUNTS}/annotation.gtf
fgrep $'\t'gene$'\t' ${INPUTS}/annotation.gtf \
    | ${PIPELINE}/scripts/sanitize-gtf-for-featureCounts \
	  > ${COUNTS}/annotation.gtf

(
    for i in $SAMPLES_INDICES ; do

	echo featureCounts $FEATURECOUNTS_ARGS \
    	     -a ${COUNTS}/annotation.gtf \
    	     -o ${COUNTS}/counts_${SAMPLES_NAME[$i]}.txt \
    	     ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.bam

    done
) | run_commands_from_stdin

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
