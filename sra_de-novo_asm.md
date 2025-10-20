# De novo genome assembly
We need to download the data first. `sratools` is there for this purpose. We need to have this tool available for us. Better to use conda/mamba to manage this.

```yaml
name: sratools
channels:
  - bioconda
  - conda-forge
  - defaults
  - nanoporetech
  - r
dependencies:
  # Workflow management
  - sra-tools=2.11.0
```
N.B. This `sra-tools` version is old. My HPC cluster has an old glibc. So I had to use this old version. You might be able to use the newest version. That's why conda/mamba is nice to have the flexibility to use/install any version of a tool.

Run `conda env create -f environment.yaml` to create the env.
Now activate and test it:
```bash
conda activate sratools
fastq-dump --version
fastq-dump --stdout -X 2 SRR390728
```
It shows:
```bash
"fastq-dump" version 2.11.0

Read 2 spots for SRR390728
Written 2 spots for SRR390728
@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
```
So, everything is fine. It is working.

## Use it now for download
```bash
prefetch SRR35821421
#prefetch SRR34300724
```
The SRA run will be downloaded by prefetch. We need to dump it to fastq files now:
```bash
fasterq-dump --split-files SRR35821421
```
It will create two fastq files (since it was a paired read file), namely, `SRR35821421_1.fastq`, `SRR35821421_2.fastq` at the current directory. We could specify a directory using `--outdir /path/to/output`. See help for this command to know more, like this: `fasterq-dump --help`.

## Clean-up
The .sra file is not eeded after we get the fastq files. We can delete it (or the folder) safely.
```bash
rm SRR35821421/SRR35821421.sra
```

## De-novo genome assembly pipeline
I made a file named `run_denovo_asm.sh` at the same working directory. The content is below:
```bash
#!/bin/bash
# -----------------------------------------------------------------------------
# De Novo Assembly Pipeline (FastQC, fastp, SPAdes, QUAST)
# -----------------------------------------------------------------------------

set -e  # Exit immediately if a command exits with a non-zero status

# --- CONFIGURATION ---
# User Parameters
GENOME_NAME="SRR35821421_asm"
THREADS=32
MEMORY_GB=200

# Input
R1_RAW="SRR35821421_1.fastq"
R2_RAW="SRR35821421_2.fastq"

# Output directories
QC_RAW_DIR="${GENOME_NAME}_00_FastQC_Raw" 
QC_TRIM_DIR="${GENOME_NAME}_01_QC_fastp"

ASSEMBLY_DIR="${GENOME_NAME}_02_SPAdes_Assembly"
QUAST_DIR="${GENOME_NAME}_03_QUAST_Report"

# File paths
R1_CLEAN="${QC_TRIM_DIR}/SRR35821421_1_clean.fastq"
R2_CLEAN="${QC_TRIM_DIR}/SRR35821421_2_clean.fastq"
CONTIGS_FASTA="${ASSEMBLY_DIR}/contigs.fasta"

# Environment Name

ASSEMBLY_ENV_NAME="ont_blast_maps"

# Ensure all directories exist
mkdir -p "$QC_RAW_DIR" "$QC_TRIM_DIR" "$ASSEMBLY_DIR" "$QUAST_DIR"

echo "Pipeline started for ${GENOME_NAME}"
echo "Resources: ${THREADS} threads, ${MEMORY_GB} GB RAM"
echo "-----------------------------------------------------------------------------"

# --- Conda Environment Setup ---
# Load conda environment hook once
eval "$(conda shell.bash hook)"

# -----------------------------------------------------------------------------
# STEP 1: Raw Quality Check (FastQC) - NEW STEP
# -----------------------------------------------------------------------------
echo "--- Running FastQC (in ${ASSEMBLY_ENV_NAME} environment) ---"
conda activate "$ASSEMBLY_ENV_NAME"

fastqc \
  -o "$QC_RAW_DIR" \
  -t $THREADS \
  "$R1_RAW" \
  "$R2_RAW"

echo "✅ FastQC completed. Reports in ${QC_RAW_DIR}/"
echo "-----------------------------------------------------------------------------"

# -----------------------------------------------------------------------------
# STEP 2: Quality Control and Trimming (fastp)
# -----------------------------------------------------------------------------
# Environment remains activated from fastqc

echo "--- Running fastp (in ${ASSEMBLY_ENV_NAME} environment) ---"

fastp \
  -i "$R1_RAW" \
  -I "$R2_RAW" \
  -o "$R1_CLEAN" \
  -O "$R2_CLEAN" \
  -q 20 -u 50 --trim_poly_x --detect_adapter_for_pe \
  -h "${QC_TRIM_DIR}/fastp_report.html" \
  -j "${QC_TRIM_DIR}/fastp_report.json" \
  -w $THREADS

echo "✅ fastp completed. Clean reads in ${QC_TRIM_DIR}/"
echo "-----------------------------------------------------------------------------"

# -----------------------------------------------------------------------------
# STEP 3: De Novo Assembly (SPAdes)
# -----------------------------------------------------------------------------
echo "--- Running SPAdes (in ${ASSEMBLY_ENV_NAME} environment) ---"

spades.py \
  --pe1-1 "$R1_CLEAN" \
  --pe1-2 "$R2_CLEAN" \
  -o "$ASSEMBLY_DIR" \
  -t $THREADS \
  -m $MEMORY_GB

echo "✅ SPAdes completed. Contigs: ${CONTIGS_FASTA}"
echo "-----------------------------------------------------------------------------"

# -----------------------------------------------------------------------------
# STEP 4: Assembly Quality (QUAST)
# -----------------------------------------------------------------------------
echo "--- Running QUAST (in ${ASSEMBLY_ENV_NAME} environment) ---"

quast.py \
  "$CONTIGS_FASTA" \
  -o "$QUAST_DIR" \
  -t $THREADS

echo "✅ QUAST completed. Reports in ${QUAST_DIR}/"
echo "-----------------------------------------------------------------------------"

echo "--- PIPELINE FINISHED SUCCESSFULLY ---"
```

I need to make the file executable and then run it:
```bash
chmod +x run_denovo_asm.sh
srun --pty -p fat -c 32 --mem=250g -t 25:00:00 /bin/bash
bash run_denovo_asm.sh
```
I am running it inside an interactive session so that I can check what is happening in every step.

It is comparitively a small file and finishes real quick.