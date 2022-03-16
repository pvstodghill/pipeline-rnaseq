# pipeline-rnaseq

Pipeline for analyzing _in-vitro_ RNA-Seq experiments.

## Installing prereqs using [Conda](https://conda.io)

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
conda install -y falco
conda install -y bioconductor-deseq2
# v-- https://github.com/bioconda/bioconda-recipes/issues/12100#issuecomment-911569353
conda install -y "samtools>=1.10"
conda install -y subread
```
## Configuring the pipeline

1. Copy `config.template.bash` to `config.bash`
2. Edit `config.bash` according to your needs and local environment.

## Running the pipeline

1. `./doit00.bash`
2. `./doit01.bash`
3. ...
