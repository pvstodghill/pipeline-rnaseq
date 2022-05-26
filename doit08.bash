#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 8. Run DESeq2
# ------------------------------------------------------------------------

echo 1>&2 '# Prepping DESeq2...'

rm -rf ${DESEQ2}
mkdir -p ${DESEQ2}/temp

cp ${COUNTS}/annotation.gtf ${DESEQ2}/temp/regions.gtf

# ------------------------------------------------------------------------
# run_deseq2 - run a single comparison
# ------------------------------------------------------------------------

deseq2_opt_c=2.0 # fold-change "cutoff"
deseq2_opt_f=0.05 # q-value (fdr) cutoff
deseq2_opt_p=-1  # p-value cutoff

run_deseq2__usage=

run_deseq2__usage+="Usage: run_deseq2 [options] -t TAG CTL:CTL1 CTL:CTL2 TMT:TMT1 TMT:TMT2 ...\n"

run_deseq2__usage+="-c NUM - fold-change cutoff [$deseq2_opt_c]\n"
run_deseq2__usage+="-f NUM - q-value (fdr) cutoff [$deseq2_opt_f]\n"
run_deseq2__usage+="-h - this message\n"
run_deseq2__usage+="-p NUM - p-value cutoff [$deseq2_opt_p]\n"
run_deseq2__usage+="-t STR - TAG to be used for output\n"

function run_deseq2_usage {
    echo -e -n "$run_deseq2__usage" 1>&2
    exit "$@"
}    

function run_deseq2
{(

    echo 1>&2 "# run_deseq2 $*"
    
    local OPTIND opt_c opt_f opt_h opt_p opt_t

    opt_c=$deseq2_opt_c
    opt_f=$deseq2_opt_f
    opt_p=$deseq2_opt_p
    
    while getopts 'c:f:hp:t:' opt ; do
	case "$opt" in
	    h) opt_h=1 ;;
	    c) opt_c="$OPTARG" ;;
	    f) opt_f="$OPTARG" ;;
	    p) opt_p="$OPTARG" ;;
	    t) opt_t="$OPTARG" ;;
	    \?) run_deseq2_usage 1 ;;
	    *) echo 1>&2 "Can't happen" ; exit 1 ;;
	esac
    done
    shift $((OPTIND-1))

    if [ "$opt_h" ] ; then
	run_deseq2_usage
    fi

    if [ -z "$opt_t" ] ; then
	echo 1>&2 "Missing -t TAG"
	run_deseq2_usage 1
    fi
    
    if [ -z "$2" ] ; then
	echo 1>&2 "Missing args"
	run_deseq2_usage 1
    fi

    ARGS=
    for RAW_ARG in "$@" ; do
	ARGS+=" "
	# FIXME: assumes condition name does not contain ":"
	ARGS+="$(echo $RAW_ARG | sed -r -e "s|.*:(.*)|\1:${COUNTS}/counts_\1.txt|")"
    done

    ${PIPELINE}/scripts/make-counts-table-from-featurecounts \
	       ${DESEQ2}/temp/regions.gtf \
	       ${ARGS} \
	       > ${DESEQ2}/temp/counts_${opt_t}.txt

    ${PIPELINE}/scripts/prep-deseq2 -x -s ${PIPELINE}/scripts \
	       -F parametric \
	       -f $opt_f -p $opt_p \
	       -d ${DESEQ2}/temp \
	       -t _${opt_t} \
	       -c ${DESEQ2}/temp/counts_${opt_t}.txt \
	       "$@"

    Rscript ${PIPELINE}/scripts/run-deseq2 \
	    ${DESEQ2}/temp/params_${opt_t}.R

    cat ${DESEQ2}/temp/output-extended_${opt_t}.txt \
	| ${PIPELINE}/scripts/deseq-output2 -t \
		     ${DESEQ2}/temp/regions.gtf \
		     ${REFERENCE_ALIASES_TXT} \
		     > ${DESEQ2}/results_${opt_t}.txt

    cat ${DESEQ2}/temp/output_${opt_t}.txt \
	| ${PIPELINE}/scripts/deseq-output2 -g -c $opt_c \
		     ${DESEQ2}/temp/regions.gtf \
		     ${REFERENCE_ALIASES_TXT} \
		     > ${DESEQ2}/results_${opt_t}.gff

    set +e # if egrep matches nothing
    cat ${DESEQ2}/results_${opt_t}.gff \
	| egrep '; colour [23];' > ${DESEQ2}/changed_${opt_t}.gff
    set -e

)}

# ------------------------------------------------------------------------
# Run the user-specified analyses
# ------------------------------------------------------------------------

if [ -f config08.bash ] ; then
    . config08.bash
else
    ARGS=
    for i in $SAMPLES_INDICES ; do
	ARGS+=" ${SAMPLES_TREATMENT[$i]}:${SAMPLES_NAME[$i]}"
    done
    TAG="${SAMPLES_TREATMENT[0]}_${SAMPLES_TREATMENT[i]}"
    
    run_deseq2 -c 2.0 -t "$TAG" $ARGS
fi

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
