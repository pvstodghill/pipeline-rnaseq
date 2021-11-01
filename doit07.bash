#! /bin/bash

. doit-preamble.bash

COUNTS=data/06_counts

# ------------------------------------------------------------------------
# Step 7. Run DESeq2
# ------------------------------------------------------------------------

echo 1>&2 '# Running DESeq2...'

CHANGE_CUTOFF=2.0
TAG=_WT-5255

DESEQ2=data/07_deseq2

rm -rf ${DESEQ2}
mkdir -p ${DESEQ2}/temp

cp ${COUNTS}/annotation.gtf ${DESEQ2}/temp/regions.gtf

./scripts/make-counts-table-from-featurecounts \
    ${DESEQ2}/temp/regions.gtf \
    52551:${COUNTS}/counts_0.txt \
    52552:${COUNTS}/counts_1.txt \
    52553:${COUNTS}/counts_2.txt \
    WT1:${COUNTS}/counts_3.txt \
    WT2:${COUNTS}/counts_4.txt \
    WT3:${COUNTS}/counts_5.txt \
    > ${DESEQ2}/temp/counts.txt

# ------------------------------------------------------------------------

./scripts/prep-deseq2 -x -s ./scripts \
		      -F parametric \
		      -d ${DESEQ2}/temp \
		      -t ${TAG} \
		      -c ${DESEQ2}/temp/counts.txt \
 		     WT:WT1 WT:WT2 WT:WT3 \
 		     5255:52551 5255:52552 5255:52553


# ------------------------------------------------------------------------

Rscript ./scripts/run-deseq2 \
	${DESEQ2}/temp/params${TAG}.R

# ------------------------------------------------------------------------

cat ${DESEQ2}/temp/output-extended${TAG}.txt \
    | ./scripts/deseq-output2results > ${DESEQ2}/results${TAG}.txt

cat ${DESEQ2}/temp/output${TAG}.txt \
    | ./scripts/deseq-output2gff ${ACCESSION} ${CHANGE_CUTOFF} \
				 > ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff
set +e # if egrep matches nothing
cat ${DESEQ2}/results${TAG}_${CHANGE_CUTOFF}.gff \
    | egrep '; colour [23];' > ${DESEQ2}/changed${TAG}_${CHANGE_CUTOFF}.gff
set -e

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
