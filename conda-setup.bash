#! /bin/bash

set -e

CONDA_PREFIX=$(dirname $(dirname $(type -p conda)))
. "${CONDA_PREFIX}/etc/profile.d/conda.sh"

PACKAGES=
#PACKAGES+=" pip"

PACKAGES+=" bowtie2"
PACKAGES+=" bioconductor-deseq2"
PACKAGES+=" falco"
PACKAGES+=" fastp"
PACKAGES+=" samtools"
PACKAGES+=" subread"

if [ "$(type -p mamba)" ] ; then
    _conda=mamba
else
    _conda=conda
fi

set -x

conda env remove -y --name pipeline-rnaseq
conda create -y --name pipeline-rnaseq
conda activate pipeline-rnaseq

$_conda install -y ${PACKAGES}

#pip install FIXME
