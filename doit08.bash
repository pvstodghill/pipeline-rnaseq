#! /bin/bash

. doit-preamble.bash

INPUTS=data/00_inputs
PROFILES=data/05_profiles
BAMGFF=${PROFILES}/tmp

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

./scripts/make-counts-table-from-gffs \
    -t -f -u \
    ${STATS}/temp/regions.gff \
    52551:${BAMGFF}/aligned_0.gff.gz \
    52552:${BAMGFF}/aligned_1.gff.gz \
    52553:${BAMGFF}/aligned_2.gff.gz \
    WT1:${BAMGFF}/aligned_3.gff.gz \
    WT2:${BAMGFF}/aligned_4.gff.gz \
    WT3:${BAMGFF}/aligned_5.gff.gz \
    > ${STATS}/stats.txt

# ------------------------------------------------------------------------
# Done.
# ------------------------------------------------------------------------

echo 1>&2 '# Done.'
