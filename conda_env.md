# Conda/mamba is used for environment management
Some tools might have special dependency or we can't have some tools (or version) together. Conda environment saves us there to isolate things. And snakemake is quite flexible in activating and using tools from different environments. 

I use mamba, it is 10X faster than conda. Enen my bashrc file has alias like `conda='mamba'`. So, `conda` command is replaced by `mamba`. If I speifically need `conda` command, I need to run it like this:
`\conda activate ont_blast_maps`. The `\` escape the alias for it and use `conda` instead of `mamba`.

## How to install mamba
```bash
conda install -n base -c conda-forge mamba
```
`mamba` is installed inside my `base` environment.

## An example (broad) environment

I made an `environment.yaml` file to manage the environment better. I use this envoronment quite extensively.
```yml
name: ont_blast_maps
channels:
  - bioconda
  - conda-forge
  - defaults
  - nanoporetech
  - r
dependencies:
  # Workflow management
  - snakemake=7.32.4
  
  # Quality control tools
  - fastqc=0.11.9
  - nanoplot=1.40.2
  
  # Read preprocessing tools
  - porechop=0.2.4
  - cutadapt=4.4
  - nanofilt=2.8.0

  # classification tools
  - kraken2
  - krakentools
  - bracken
  - krona
  # BLAST and sequence manipulation
  - blast=2.16.0
  - seqtk=1.3
  - entrez-direct
  - taxonkit
  
  # Taxonomic classification
  - metamaps
  - nanosim
  
  # Assembly
  - flye=2.9.1
  - raven-assembler
  - megahit
  - spades
  - canu
  - quickmerge

  # Mapping tools
  - minimap2=2.24
  - samtools=1.16
  - bedtools=2.30.0
  - mummer
  # Annotation tools
  - prokka=1.14.6
  - diamond=2.1.4
  - quast  # For assembly QC
  - busco
  # Additional utilities
  - wget=1.20.3
  - python=3.9
  - biopython
  - pip
  - numpy
  - pandas
  - matplotlib
  - seaborn
  - scipy
  - scikit-learn
```
Well, we need to navigate to the folder where the `environment.yaml` file is, using `cd` command. Then I run `conda env create -f environment.yaml` to create the env with all the specified tools. So the environment name will be `ont_blast_maps`. 

## Using the env
Run `conda activate ont_blast_maps` to activate the environment. So all the specified tools in this environment will be available to use.

## If I forget the name(s) of the environment
Run `conda env list` to list all available/created environments. The output will be something like this:
```yml
base                  *  /data/users/md.rasheduzzaman/miniconda3
biobakery                /data/users/md.rasheduzzaman/miniconda3/envs/biobakery
busco                    /data/users/md.rasheduzzaman/miniconda3/envs/busco
checkm2_env              /data/users/md.rasheduzzaman/miniconda3/envs/checkm2_env
kat                      /data/users/md.rasheduzzaman/miniconda3/envs/kat
langchain_env            /data/users/md.rasheduzzaman/miniconda3/envs/langchain_env
latex-env                /data/users/md.rasheduzzaman/miniconda3/envs/latex-env
medaka                   /data/users/md.rasheduzzaman/miniconda3/envs/medaka
mol_work                 /data/users/md.rasheduzzaman/miniconda3/envs/mol_work
nanosim_env              /data/users/md.rasheduzzaman/miniconda3/envs/nanosim_env
ont_blast_maps           /data/users/md.rasheduzzaman/miniconda3/envs/ont_blast_maps
ont_meta                 /data/users/md.rasheduzzaman/miniconda3/envs/ont_meta
ont_sim                  /data/users/md.rasheduzzaman/miniconda3/envs/ont_sim
pytorch_env              /data/users/md.rasheduzzaman/miniconda3/envs/pytorch_env
```
Then activate using the name. For example, `conda activate ont_blast_maps`.
To deactivate, run `conda deactivate`. Here, we don't need the name of the environment. Since we were already inside the environment, we will just be out of it.

## Updating an env
If we make some modification (removal or addition of new tools) in the `environment.yaml` file, we can update our environment. 2 ways:
- 1. If the env is already activated:
```bash
conda env update --file environment.yaml --prune
```
N.B. The `--prune` flag tells conda/mamba to remove any packages that are currently in the environment but are no longer listed in the updated `environment.yaml` file.
- 2. If the env is not activated:
```bash
conda env update --name ont_blast_maps --file environment.yaml --prune
```

## Remove/delete an environment
Deactivate the env. It's not possible to remove an activated environment.
```bash
conda deactivate
```
Now, remove it running:
```bash
conda env remove -n <environment_name> --all
```