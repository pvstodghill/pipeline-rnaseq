#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# Step 7. Run DESeq2
# ------------------------------------------------------------------------

echo 1>&2 '# Running DESeq2...'

CHANGE_CUTOFF=2.0
TAG=_${SAMPLES_TREATMENT[0]}-${SAMPLES_TREATMENT[3]}

rm -rf ${DESEQ2}
mkdir -p ${DESEQ2}/temp

cp ${COUNTS}/annotation.gtf ${DESEQ2}/temp/regions.gtf

# FIXME: parameterize
${PIPELINE}/scripts/make-counts-table-from-featurecounts \
    ${DESEQ2}/temp/regions.gtf \
    ${SAMPLES_NAME[0]}:${COUNTS}/counts_${SAMPLES_NAME[0]}.txt \
    ${SAMPLES_NAME[1]}:${COUNTS}/counts_${SAMPLES_NAME[1]}.txt \
    ${SAMPLES_NAME[2]}:${COUNTS}/counts_${SAMPLES_NAME[2]}.txt \
    ${SAMPLES_NAME[3]}:${COUNTS}/counts_${SAMPLES_NAME[3]}.txt \
    ${SAMPLES_NAME[4]}:${COUNTS}/counts_${SAMPLES_NAME[4]}.txt \
    ${SAMPLES_NAME[5]}:${COUNTS}/counts_${SAMPLES_NAME[5]}.txt \
    > ${DESEQ2}/temp/counts.txt

# ------------------------------------------------------------------------

# FIXME: parameterize
${PIPELINE}/scripts/prep-deseq2 -x -s ${PIPELINE}/scripts \
		      -F parametric \
		      -d ${DESEQ2}/temp \
		      -t ${TAG} \
		      -c ${DESEQ2}/temp/counts.txt \
		      ${SAMPLES_TREATMENT[0]}:${SAMPLES_NAME[0]} \
		      ${SAMPLES_TREATMENT[1]}:${SAMPLES_NAME[1]} \
		      ${SAMPLES_TREATMENT[2]}:${SAMPLES_NAME[2]} \
		      ${SAMPLES_TREATMENT[3]}:${SAMPLES_NAME[3]} \
		      ${SAMPLES_TREATMENT[4]}:${SAMPLES_NAME[4]} \
		      ${SAMPLES_TREATMENT[5]}:${SAMPLES_NAME[5]}


# ------------------------------------------------------------------------

Rscript ${PIPELINE}/scripts/run-deseq2 \
	${DESEQ2}/temp/params${TAG}.R

# ------------------------------------------------------------------------

cat ${DESEQ2}/temp/output-extended${TAG}.txt \
    | ${PIPELINE}/scripts/deseq-output2 -t \
			      ${DESEQ2}/temp/regions.gtf \
			      ${REFERENCE_ALIASES_TXT} \
			      > ${DESEQ2}/results${TAG}.txt

cat ${DESEQ2}/temp/output${TAG}.txt \
    | ${PIPELINE}/scripts/deseq-output2 -g -c ${CHANGE_CUTOFF} \
			      ${DESEQ2}/temp/regions.gtf \
			      ${REFERENCE_ALIASES_TXT} \
			      > ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff

set +e # if egrep matches nothing
cat ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff \
    | egrep '; colour [23];' > ${DESEQ2}/changed${TAG}_${CHANGE_CUTOFF}.gff
set -e

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
