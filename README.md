# pipeline-rnaseq

Pipeline for analyzing _in-vitro_ RNA-Seq experiments.

Setting a [Conda](https://conda.io) environment for the pipeline,

```
# conda env remove -y --name rnaseq
conda create -y --name rnaseq
conda activate rnaseq

# v-- https://www.biostars.org/p/455593/#479908
# v-- https://www.biostars.org/p/455593/#9468633
conda config --add channels bioconda
conda config --add channels conda-forge

conda install -y bowtie2
# v-- for me 0.23.x hangs
conda install -y fastp=0.22.0
conda install -y fastqc
conda install -y bioconductor-deseq2
# v-- https://github.com/bioconda/bioconda-recipes/issues/12100#issuecomment-911569353
conda install -y "samtools>=1.10"
conda install -y subread
```
