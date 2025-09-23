# pipeline-rnaseq

Pipeline for analyzing _in-vitro_ RNA-Seq experiments.

## Cloning the repo

This pipeline using Git submodules. The easiest way to clone this repo (with a recent version of `git`) is

```
git clone --recurse-submodules https://github.com/pvstodghill/pipeline-rnaseq.git
```

## Installing prereqs

One of the following:

<!-- - [Docker](https://www.docker.com/) -->
<!-- - [Singularity](https://sylabs.io/) -->
<!-- - [Apptainer](https://apptainer.org/) -->
- [Conda](https://conda.io)

You will also need,

- [Snakemake](https://snakemake.readthedocs.io/)
- [Perl](https://www.perl.org/)
- [Perl's YAML module](https://metacpan.org/dist/YAML)

## Configuring the pipeline

**Create the configuration files**

For the `example` data,

1. Copy `example/config.yaml` to  `config.yaml`.

To run the pipeline on your own data,

1. Copy `config.template.yaml` to `config.yaml`.  Edit `config.yaml` according to your needs and local environment.

## Running the pipeline

To run the pipeline using local copies of the software components:

~~~
snakemake --cores 32
~~~

Change `32` to the number of threads you want to us.

To run the pipeline using [Conda](https://conda.io) to provide software components:

~~~
snakemake --use-conda --cores 32
~~~

To run the pipeline using [Mamba](https://mamba.readthedocs.io) to provide software components:

~~~
snakemake --use-conda --conda-frontend mamba
~~~

## Software components

The following software is used by this pipeline,

- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [DESeq2](https://bioconductor.org/packages/DESeq2/)
- [Falco](https://github.com/smithlabcode/falco)
- [Fastp](https://github.com/OpenGene/fastp)
- `featureCount` from [Subread](http://subread.sourceforge.net/)
- [R](https://www.r-project.org/)
- [Samtools](https://github.com/samtools/samtools)
