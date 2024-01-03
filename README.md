# cidc_rnaseq

*Example workflow framework to assist with developing pipelines in Snakemake. This framework is currently set up to perform the following across sample files: read pre-processing (adapter trimming and quality filtering), assessment of sample contamination, read QC, reference genome alignment, and obtainment of alignment metrics.*

# Getting started

## Requirements:

The pipeline requires a computer running Linux (Ubuntu 20). This pipeline is currently set up to run on a GCP VM instance. First generate a VM instance and perform the following steps after logging into the VM instance.

Most software dependencies are managed using *conda*. To install conda, please install [miniconda3](https://conda.io/miniconda.html) and refer to installation [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
Accept the license agreement during installation, and it is recommended to allow the conda installer to prepend its path to user's .bashrc file when asked.

## Conda installation and setup:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Check if conda has successfully installed.

```
conda -h
```

If conda has installed correctly, the following output will be present.
If the below output does not generate, close and reopen the terminal.

```
$ conda
usage: conda [-h] [--no-plugins] [-V] COMMAND ...

conda is a tool for managing and deploying applications, environments and packages.

options:
  -h, --help          Show this help message and exit.
  --no-plugins        Disable all plugins that are not built into conda.
  -V, --version       Show the conda version number and exit.

commands:
  The following built-in and plugins subcommands are available.

  COMMAND
    clean             Remove unused packages and caches.
    compare           Compare packages between conda environments.
    config            Modify configuration values in .condarc.
    content-trust     See `conda content-trust --help`.
    create            Create a new conda environment from a list of specified packages.
..............
```

Next install mamba and initialize a mamba environment.
```
conda install -n base -c conda-forge mamba
mamba init
```

Add conda channels and set priorities.
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

## Pipeline setup and execution:
Clone git repository to the location where you want to run your analysis and create the conda environment that will be used to run the pipeline.
```
git clone https://github.com/NCI-CIDC/cidc_pipeline_template.git
cd cidc_pipeline_template
## Creates environment and installs dependencies automatically
conda env create -f=environment.yaml
conda activate cidc_pipeline_template
```
Pipeline configuration (config.yaml)
* This is where important workflow info is set / stored, This includes paths to reference genomes and sample metadata (sample info and gcloud location of input fastq files), as well as workflow output and source directories and options.


Test the pipeline and print jobs by performing a dryrun:
```
## Print workflow jobs from start to finish (<ncores> the number of cores for entire workflow):
snakemake --use-conda --cores <ncores> --rerun-incomplete --forceall --dryrun

## Print workflow jobs resuming to finish incomplete rules (<ncores> the number of cores for entire workflow):
snakemake --use-conda --cores <ncores> --rerun-incomplete --dryrun
```

Execute the pipeline:
```
## Execute workflow start to finish (<ncores> the number of cores for entire workflow):
snakemake --use-conda --cores <ncores> --rerun-incomplete --forceall

## Execute workflow resuming to finish incomplete rules (<ncores> the number of cores for entire workflow):
snakemake --use-conda --cores <ncores> --rerun-incomplete
```

## Please reach out to Sami Cherikh (cherikhsr@nih.gov) with any questions or to report any issues concerning this framework.
