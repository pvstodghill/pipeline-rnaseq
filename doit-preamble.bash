#! /bin/bash

HOWTO="./scripts/howto -q -T data/tmp -f packages.yaml"
THREADS=$(nproc --all)

export LC_ALL=C

# ------------------------------------------------------------------------

. config.bash

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

if [ "$PACKAGES_FROM" = conda ] ; then
    if [ -z "$CONDA_EXE" ] ; then
	CONDA_EXE=$(type -p conda)
    fi
fi

# In order to help test portability, I eliminate all of my
# personalizations from the PATH, etc.
if [ "$PVSE" ] ; then
    export PATH=/usr/local/bin:/usr/bin:/bin
    export PERL5LIB=
    export PERL_LOCAL_LIB_ROOT=
    export PERL_MB_OPT=
    export PERL_MM_OPT=
    export PYTHONPATH=
fi

case X"$PACKAGES_FROM"X in
    XcondaX)
	CONDA_PREFIX=$(dirname $(dirname $CONDA_EXE))
	. "${CONDA_PREFIX}/etc/profile.d/conda.sh"
	conda activate $PACKAGES_ENV

	;;
    XhowtoX|XstubsX)
	export PATH=$(dirname ${BASH_SOURCE[0]})/stubs:"$PATH"
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

PARALLEL_PATH=$(type -p parallel)

function run_parallel {
    if [ "$PARALLEL_PATH" ] ; then
	eval $PARALLEL_PATH -j ${THREADS} -kv
    else
	bash -x
    fi
}

# ------------------------------------------------------------------------

set -e
set -o pipefail

