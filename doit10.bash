#! /bin/bash

. $(dirname ${BASH_SOURCE[0]})/doit-preamble.bash

# ------------------------------------------------------------------------
# software version
# ------------------------------------------------------------------------

(
    set -x
    bowtie2 --version | head -n1
    Rscript --version
    ${PIPELINE}/scripts/print-r-library-versions DESeq2
    falco --version
    fastp --version
    samtools --version | head -n1
    featureCounts -v
)


# ------------------------------------------------------------------------
# Print final git info
# ------------------------------------------------------------------------

cd ${PIPELINE}

echo ''
(
    set -x
    git status
)
echo ''
(
    set -x
    git log -n1
)
