# To basecall from `pod5` file
The ONT machines do basecalling themselves. But sometimes it might be necessary to basecall externally because of machine error or out of space issue, for example. Here is how to do it.

## `pod5` to `fastq` file
We are going to use dorado. Dorado makes `pod5` to `bam` files. We need to use `--emit-fastq` argument to make `fastq` file in one go. Also, we need to do it on gpu node of our HPC cluster to make it faster.

But to do this, we need to download a model first. We are going to download the `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` model. Let's see how to do it.
```bash
cd /home/md.rasheduzzaman/data/molestus_Meta
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```
So, we will have the model named `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` inside our `molestus_Meta` directory. Now, we can use it.

We have a folder `data/Cx_Commensal1/test_pod5_pass` (relative to the `molestus_Meta` directors) with all the `pod5` files for a sequencing run/sample. Let's do basecalling for it.

We will use snakemake for it. The file name is `Cx_bc1.smk`, inside the `home/md.rasheduzzaman/data/molestus_Meta/dorado_bc` folder. The content is below:
```bash
#be in /home/md.rasheduzzaman/data/molestus_Meta directory
input_dir = "data/Cx_Commensal1/test_pod5_pass"
output_fastq = "data/Cx_Commensal1_1/basecalled/basecalled.fastq"
##GCCGGAGCTCTGCAGATA --primer barcode
# Define the rule for basecalling
rule all:
    input:
        output_fastq

rule basecall:
    input:
        pod5_dir = input_dir
    output:
        fastq = output_fastq
    log:
        "logs/basecall1.log"
    conda:
        "mol_work"
    threads: 32
    resources:
        mem_mb=200000
    shell:
        """
        dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v5.0.0 {input.pod5_dir} --device cuda:0 --emit-fastq > {output.fastq} 2> {log}
        """
```
Now, we can submit this job to our "gpu" node. Let's make the sbatch submission script. We named it `submit_dorado_bc1.sh`
```bash
#!/bin/bash
#SBATCH --job-name=Cx_dorado_bc1
#SBATCH --output=logs/Cx_dorado_bc1_%j.out.txt
#SBATCH --error=logs/Cx_dorado_bc1_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=200000
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

## jobid 1481800
# Run Snakemake
snakemake -j 1 --cores 32 -s dorado_bc/Cx_bc1.smk --latency-wait 120 --rerun-incomplete
```
Now, submit:
```bash
sbatch dorado_bc/submit_dorado_bc1.sh
```
It will run for a while.