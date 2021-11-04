#! /bin/bash

. doit-preamble.bash

BOWTIE2=data/04_bowtie2

# ------------------------------------------------------------------------
# Step 5. Generating GFF files + Making profiles
# ------------------------------------------------------------------------

PROFILES=data/05_profiles

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

	echo "samtools view -h ${BOWTIE2}/aligned_$i.bam \
	| ./scripts/sam2profiles $opt_12 $opt_r \
				 -e -s -n \
				 -d ${PROFILES} \
				 -t ${name}"

    done
) | run_commands

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
