## kraken2 database making
Let's make a file named `kkn_bd_making.sh` like this:
```bash
#!/bin/bash
## Kraken2 DB building for virus, archaea, bacteria, and Culex host

# --- 1. Download Taxonomy ---
# This is mandatory and must be run before adding any sequences.
kraken2-build --download-taxonomy --db ../../kkn_db2 --threads 32

# --- 2. Download Standard Libraries ---
# Download RefSeq for common microbial domains.
kraken2-build --download-library viral --db ../../kkn_db2 --threads 32
kraken2-build --download-library bacteria --db ../../kkn_db2 --threads 32
kraken2-build --download-library archaea --db ../../kkn_db2 --threads 32

# --- 3. Download and Prepare Custom Host (Culex) Genome ---
# Download the specific Culex genome (GCF_016801865.2).
datasets download genome accession GCF_016801865.2 --include genome --filename GCF_016801865.2_RefSeq.zip

# Unzip the downloaded file. This creates the 'ncbi_dataset' directory structure.
unzip GCF_016801865.2_RefSeq.zip

# --- 4. Add Custom Genome to Library ---
# CAUTION: The FASTA file name inside the zip often changes.
# Using the wildcard '*' for the filename is safer to ensure it finds the correct '_genomic.fna' file.
kraken2-build --add-to-library ncbi_dataset/data/GCF_016801865.2/*_genomic.fna --db ../../kkn_db2 --threads 32

# --- 5. Build the Final Database ---
# This step creates the hash table from all combined library files. It is the most resource-intensive step.
kraken2-build --build --db ../../kkn_db2 --threads 32

# sbatch command (commented out)
# sbatch -p fat -c 32 --mem=280G --time=30:00:00 --wrap "bash kkn_bd_making.sh"
```
Now make it executable and submit it:
```bash
chmod +x kkn_bd_making.sh
sbatch -p fat -c 32 --mem=280G --time=30:00:00 --wrap "bash kkn_bd_making.sh"
```
We will have our kraken2 database ready.