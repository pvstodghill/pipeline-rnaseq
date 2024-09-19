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

- [Perl](https://www.perl.org/)
- [Perl's YAML module](https://metacpan.org/dist/YAML)

## Configuring the pipeline

**Create the configuration files**

For the `example` data,

1. Copy `example/config.mk` to  `config.mk`.

To run the pipeline on your own data,

1. Copy `config.template.mk` to `config.mk`.  Edit `config.mk` according to your needs and local environment.

**Choose how to access packages.**

<!-- Do one of the following, -->

<!-- 1. For strong reproducibility, use with Docker packages with explicit -->
<!--    versions: -->

<!--     * Edit `config.bash` -->
<!--     * Uncomment `PACKAGES_FROM=howto` -->


<!-- 1. For convenience(?), flexibility(?), use the latest versions of -->
<!--    Conda packages. -->

<!--     * Create a Conda environment with the necessary packages, perhaps -->
<!--       using `conda-setup.bash`. -->
<!--     * Edit `config.bash` -->
<!--     * Uncomment `PACKAGES_FROM=conda` -->
<!--     * Uncomment and set `CONDA_ENV=...` -->

This pipeline is currently hard coded to use [`howto`](https://github.com/pvstodghill/howto/) for accessing packages. The package manifest can be found in `packages.yaml`, which currently pulls packages using [Conda](https://conda.io).

## Running the pipeline

1. `make`

## Software components

The following software is used by this pipeline,

- [Bowtie2](https://github.com/BenLangmead/bowtie2)
- [DESeq2](https://bioconductor.org/packages/DESeq2/)
- [Falco](https://github.com/smithlabcode/falco)
- [Fastp](https://github.com/OpenGene/fastp)
- `featureCount` from [Subread](http://subread.sourceforge.net/)
- [R](https://www.r-project.org/)
- [Samtools](https://github.com/samtools/samtools)
