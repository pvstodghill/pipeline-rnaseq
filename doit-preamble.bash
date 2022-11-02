#! /bin/bash

# In order to help test portability, I eliminate all of my
# personalizations from the PATH, etc.
if [ "$PVSE" ] ; then
    HOWTO_CONDA_CMD="${HOWTO_CONDA_CMD:-$(type -p mamba)}"
    HOWTO_CONDA_CMD="${HOWTO_CONDA_CMD:-$(type -p conda)}"
    if [ "$HOWTO_CONDA_CMD" ] ; then
	export HOWTO_CONDA_CMD
    fi
    export PATH=/usr/local/bin:/usr/bin:/bin
    export PERL5LIB=
    export PERL_LOCAL_LIB_ROOT=
    export PERL_MB_OPT=
    export PERL_MM_OPT=
    export PYTHONPATH=
fi

# ------------------------------------------------------------------------

# yuck. ugly.

if [ -e /programs/docker/bin/docker1 ] ; then
    export HOWTO_DOCKER_CMD=/programs/docker/bin/docker1
    THREADS=32
else
    THREADS=$(nproc --all)
fi

if [ -e /programs/parallel/bin/parallel ] ; then
    PARALLEL_CMD=/programs/parallel/bin/parallel
fi

# ------------------------------------------------------------------------

# Uncomment to get packages from HOWTO
PACKAGES_FROM=howto

# Override the default number of threads (nproc --all)
#THREADS=32


# These vars are used in parameters to stubs/*, so they cannot be
# `realpath`'ed.
PIPELINE=$(dirname ${BASH_SOURCE[0]})
# v-- can be specified externally
DATA=${DATA:-data}

# ------------------------------------------------------------------------

. config.bash

# ------------------------------------------------------------------------

export HOWTO_MOUNT_DIR=$(realpath $(${PIPELINE}/howto/find-closest-ancester-dir . ${DATA} ${PIPELINE}))
export HOWTO_TMPDIR=$(realpath ${DATA})/tmp

if [ "$PACKAGES_FROM" = conda ] ; then
    if [ -z "$CONDA_EXE" ] ; then
	CONDA_EXE=$(type -p conda)
    fi
fi

case X"$PACKAGES_FROM"X in
    XcondaX)
	CONDA_PREFIX=$(dirname $(dirname $CONDA_EXE))
	. "${CONDA_PREFIX}/etc/profile.d/conda.sh"
	conda activate $CONDA_ENV || exit 1
	;;
    XX|XhowtoX|XstubsX)
	export PATH=$(realpath ${PIPELINE})/stubs:"$PATH"
	;;
    XnativeX)
	: nothing
	;;
    XX)
	echo 1>&2 "\$PACKAGES_FROM is not set"
	exit 1
	;;
    X*X)
	echo 1>&2 "\$PACKAGES_FROM is recognized: $PACKAGES_FROM"
	exit 1
	;;
    *)
	echo 1>&2 "Cannot happen"
	exit 1
esac

# ------------------------------------------------------------------------

if [ -z "$PARALLEL_CMD" ] ; then
    PARALLEL_CMD="$(type -p parallel)"
fi

# Usage: generate_commands_to_stdin | run_commands_from_stdin
function run_commands_from_stdin {
    if [ "$PARALLEL_CMD" -a "$THREADS" -gt 1 ] ; then
	eval $PARALLEL_CMD -j ${THREADS} -kv
    else
	bash -x
    fi
}

# ------------------------------------------------------------------------

set -e
set -o pipefail

export LC_ALL=C

# ------------------------------------------------------------------------

SAMPLES_INDICES=""
for i in $(seq 0 9) ; do
    if [ "${SAMPLES_R1[$i]}" ] ; then
	SAMPLES_INDICES+=" $i"
    fi
done

if [ "${SAMPLES_R2[0]}" ] ; then
    PE=1
fi

# ------------------------------------------------------------------------

function init_FEATURECOUNTS_ARGS {
    FEATURECOUNTS_ARGS=

    FEATURECOUNTS_ARGS+=" -O" # Assign reads to all their overlapping features

    #FEATURECOUNTS_ARGS+=" -M" # all multi-mapping reads reported alignments will be counted
    ##FEATURECOUNTS_ARGS+=" --fraction" # Assign fractional counts to features
    #FEATURECOUNTS_ARGS+=" --primary" # Count primary alignments only !(0x100)

    if [ "$PE" ] ; then

	V="$(featureCounts -v 2>&1 | egrep . | sed -e 's/.* v//')"
	case x"$V"x in
	    x2.0.2x|x2.0.3x)
		FEATURECOUNTS_ARGS+=" -p" # input data contains paired-end reads. (>=v2.0.2)
		FEATURECOUNTS_ARGS+=" --countReadPairs" # Count read pairs (fragments) instead of reads
		;;
	    x2.0.1x)
		FEATURECOUNTS_ARGS+=" -p" # fragments (or pairs) will be counted instead of reads (<v2.0.2)
		;;
	    *)
		echo 1>&2 Unknown featureCounts version: "$V"
		exit 1
	esac

	FEATURECOUNTS_ARGS+=" -B" # Only count read pairs that have both ends aligned.
	FEATURECOUNTS_ARGS+=" -P" # Check validity of paired-end distance
	FEATURECOUNTS_ARGS+=" -C" # Only count concordant reads
    fi

}

# ------------------------------------------------------------------------

INPUTS=${DATA}/00_inputs
FALCO1=${DATA}/01_falco
FASTP=${DATA}/02_fastp
FALCO2=${DATA}/03_falco
BOWTIE2=${DATA}/04_bowtie2
STRAND=${DATA}/05_strand
PROFILES=${DATA}/06_profiles
COUNTS=${DATA}/07_counts
DESEQ2=${DATA}/08_deseq2
STATS=${DATA}/09_stats
