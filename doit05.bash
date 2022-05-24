#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 5. Generating GFF files + Making profiles
# ------------------------------------------------------------------------

rm -rf ${PROFILES}
mkdir -p ${PROFILES}

echo 1>&2 '# Making profiles'

if [ "$PE" ] ; then
    opt_12=-2
else
    opt_12=-1
fi

case X"$ORIENTATION"X in
    XforwardX) opt_r= ;;
    XreverseX) opt_r=-r ;;
    X*X) echo 1>&2 cannot happen ; exit 1
esac

(
    for i in $SAMPLES_INDICES ; do

	name=${SAMPLES_NAME[i]}

	echo "samtools view -h ${BOWTIE2}/aligned_${SAMPLES_NAME[$i]}.bam \
	| ${PIPELINE}/scripts/sam2profiles $opt_12 $opt_r \
				 -e -s -n \
				 -d ${PROFILES} \
				 -t ${name}"

    done
) | run_commands_from_stdin

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
