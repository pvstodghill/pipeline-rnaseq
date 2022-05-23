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

**Create the package manifest.**

 Do one of the following,

1. For strong reproducibility, use with Docker packages with explicit
   versions

```
cp packages.docker.yaml packages.yaml
```

1. For convenience(?), flexibility(?), use the latest versions of
   Conda packages.

```
cp packages.docker.yaml packages.yaml
```

1. Use whatever combination of package managers and versions you want.

```
$EDITOR packages.yaml
```

**Create the configuration files**

For the `example` data,

1. Copy `example/config.bash` to  `config.bash`.

To run the pipeline on your own data,

1. Copy `config.template.bash` to `config.bash`.
2. Edit `config.bash` according to your needs and local environment.

## Running the pipeline

1. `./doit00.bash`
2. `./doit01.bash`
3. ...
