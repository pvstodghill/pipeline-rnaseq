# pipeline-rnaseq

Pipeline for analyzing _in-vitro_ RNA-Seq experiments.

## Cloning the repo

This pipeline using Git submodules. The easiest way to clone this repo (with a recent version of `git`) is

```
git clone --recurse-submodules https://github.com/pvstodghill/pipeline-rnaseq.git
```

## Installing prereqs

One of the following:

- [Docker](https://www.docker.com/)
- [Singularity](https://sylabs.io/)
- [Apptainer](https://apptainer.org/)
- [Conda](https://conda.io)

You will also need,

- [Perl](https://www.perl.org/)
- [Perl's YAML module](https://metacpan.org/dist/YAML)

## Configuring the pipeline

**Create the configuration files**

For the `example` data,

1. Copy `example/config.bash` to  `config.bash`.
2. Optionally, copy `example/config08.bash` to  `config08.bash`.

To run the pipeline on your own data,

1. Copy `config.template.bash` to `config.bash`.  Edit `config.bash`
   according to your needs and local environment.
2. Optionally,  `config08.template.bash` to `config08.bash` and edit.

**Choose how to access packages.**

 Do one of the following,

1. For strong reproducibility, use with Docker packages with explicit
   versions:

    * Edit `config.bash`
    * Uncomment `PACKAGES_FROM=howto`


1. For convenience(?), flexibility(?), use the latest versions of
   Conda packages.

    * Create a Conda environment with the necessary packages, perhaps
      using `conda-setup.bash`.
    * Edit `config.bash`
    * Uncomment `PACKAGES_FROM=conda`
    * Uncomment and set `CONDA_ENV=...`

## Running the pipeline

1. `./doit00.bash`
2. `./doit01.bash`
3. ...
