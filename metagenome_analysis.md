# Metagenomic Classification from Illumina Novaseq (PE reads), 150bp
We need a `config.yaml` file to store file paths, parameters, etc. This configuration file is the easiest on eto understand. So, we always should make a good config file. It helps! Snakemake file (rules) are somewhat cryptic.

The `config.yaml` file content:
```bash
# Configuration file for RNA metagenomic analysis pipeline
# NovaSeq PE 150 data

# Sample name (without path and without *1/*2 suffix)
sample_name: "NG-A5666_7Usu_Cx_libLAQ4872"

# Number of threads to use
threads: 16

# Kraken2 database path (REQUIRED - update this path)
kraken_db: "../../k2_db2/"

# Kraken2 confidence threshold (0.0 to 1.0)
# 0.0 = no confidence threshold (default)
# 0.1 = more conservative classification
kraken_confidence: 0.05

# Assembly parameters for RNA
min_contig_length: 200  # Shorter for RNA transcripts

# Memory for rnaSPAdes (in GB)
spades_memory: 100

# Bracken parameters
bracken_read_length: 400
bracken_threshold: 10

# fastp trimming parameters
min_length_after_trim: 50
qualified_quality_phred: 15
cut_mean_quality: 20
```
Now, we need the sankefile with all the rules to do our analysis. Each rule means the step we take in our step-by-step analysis. The filename is `Snakefile` without any extension. We could use with extension though. In that case, we could name it like `metagenome.smk` or similar. The `.smk` ending/suffix means it is a snakemake file. But snakemake command takes `Snakefile` from our current/working directory by default. So, if we are having only one snakemake file, naming it to `Snakefile` is easier. But we can specify even this file name as well. You know what is a default value and how to override it, right? 
```bash
# Snakemake pipeline for RNA metagenomic analysis
# NovaSeq PE 150 data processing
# Two approaches: 1) Kraken2 on reads, 2) Assembly -> Kraken2 on contigs 3) Krona plot for visualization

import os
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Get sample info from config
SAMPLE = config.get("sample_name", "NG-A5666_7Usu_Cx_libLAQ4872")
THREADS = config.get("threads", 16)
KRAKEN_DB = config["kraken_db"]
KRAKEN_CONFIDENCE = config.get("kraken_confidence", 0.05)
MIN_CONTIG_LENGTH = config.get("min_contig_length", 200)  # Shorter for RNA

# Define data directory and input files
DATA_DIR = "../../data/NGS"

# Define all outputs with expand for the sample
rule all:
    input:
        # QC reports
        expand("results/fastqc/raw/{sample}_multiqc_report.html", sample=[SAMPLE]),
        expand("results/fastqc/trimmed/{sample}_multiqc_report.html", sample=[SAMPLE]),
        
        # Host removal stats (if needed for RNA)
        expand("results/trimmed/{sample}_trimming_stats.txt", sample=[SAMPLE]),
        
        # Approach 1: Kraken2 on reads
        expand("results/kraken2/reads/{sample}_kraken_report.txt", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_kraken_output.txt", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_bracken_species.txt", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_bracken_genus.txt", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_classified_1.fq", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_classified_2.fq", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_unclassified_1.fq", sample=[SAMPLE]),
        expand("results/kraken2/reads/{sample}_unclassified_2.fq", sample=[SAMPLE]),

        # Approach 2: Assembly -> Kraken2
        expand("results/assembly/{sample}_transcripts.fasta", sample=[SAMPLE]),
        expand("results/assembly/quast/{sample}_report.html", sample=[SAMPLE]),
        expand("results/kraken2/contigs/{sample}_kraken_report.txt", sample=[SAMPLE]),
        expand("results/kraken2/contigs/{sample}_kraken_output.txt", sample=[SAMPLE]),
        
        # Visualization
        expand("results/kraken2/reads/{sample}_krona.html", sample=[SAMPLE]),
        expand("results/kraken2/contigs/{sample}_krona.html", sample=[SAMPLE]),
        
        # Summary statistics
        expand("results/summary/{sample}_pipeline_stats.txt", sample=[SAMPLE])

# Rule 1: FastQC on raw reads
rule fastqc_raw:
    input:
        r1=f"{DATA_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{DATA_DIR}/{{sample}}_2.fastq.gz"
    output:
        html1="results/fastqc/raw/{sample}_1_fastqc.html",
        html2="results/fastqc/raw/{sample}_2_fastqc.html",
        zip1="results/fastqc/raw/{sample}_1_fastqc.zip",
        zip2="results/fastqc/raw/{sample}_2_fastqc.zip"
    params:
        outdir="results/fastqc/raw"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}
        """

# Rule 2: MultiQC for raw reads
rule multiqc_raw:
    input:
        expand("results/fastqc/raw/{{sample}}_{read}_fastqc.zip", read=["1", "2"])
    output:
        "results/fastqc/raw/{sample}_multiqc_report.html"
    params:
        outdir="results/fastqc/raw"
    threads: 2
    shell:
        """
        multiqc -f -o {params.outdir} -n {wildcards.sample}_multiqc_report.html {params.outdir}/*{wildcards.sample}*
        """

# Rule 3: Trim reads with fastp (better for NovaSeq and RNA-seq)
rule fastp_trim:
    input:
        r1=f"{DATA_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{DATA_DIR}/{{sample}}_2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2="results/trimmed/{sample}_2_trimmed.fastq.gz",
        json="results/trimmed/{sample}_fastp.json",
        html="results/trimmed/{sample}_fastp.html"
    params:
        # NovaSeq specific parameters
        qualified_quality_phred=15,
        length_required=50,  # Minimum length after trimming
        cut_window_size=4,
        cut_mean_quality=20,
        # Poly-X tail trimming for RNA-seq
        trim_poly_x="--trim_poly_x",
        # NovaSeq adapter sequences
        adapter_sequence="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
        adapter_sequence_r2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    threads: THREADS
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --length_required {params.length_required} \
            --cut_front --cut_tail \
            --cut_window_size {params.cut_window_size} \
            --cut_mean_quality {params.cut_mean_quality} \
            --adapter_sequence {params.adapter_sequence} \
            --adapter_sequence_r2 {params.adapter_sequence_r2} \
            {params.trim_poly_x} \
            --correction \
            --overrepresentation_analysis \
            --thread {threads} \
            2> {log}
        """

# Rule 4: Extract trimming statistics
rule trimming_stats:
    input:
        "results/trimmed/{sample}_fastp.json"
    output:
        "results/trimmed/{sample}_trimming_stats.txt"
    run:
        import json
        with open(input[0], 'r') as f:
            data = json.load(f)
        
        with open(output[0], 'w') as out:
            out.write(f'=== Trimming Statistics for {wildcards.sample} ===\n\n')
            out.write('Before filtering:\n')
            out.write(f"  Total reads: {data['summary']['before_filtering']['total_reads']:,}\n")
            out.write(f"  Total bases: {data['summary']['before_filtering']['total_bases']:,}\n")
            out.write(f"  Q20 rate: {data['summary']['before_filtering']['q20_rate']:.3f}\n")
            out.write(f"  Q30 rate: {data['summary']['before_filtering']['q30_rate']:.3f}\n\n")
            
            out.write('After filtering:\n')
            out.write(f"  Total reads: {data['summary']['after_filtering']['total_reads']:,}\n")
            out.write(f"  Total bases: {data['summary']['after_filtering']['total_bases']:,}\n")
            out.write(f"  Q20 rate: {data['summary']['after_filtering']['q20_rate']:.3f}\n")
            out.write(f"  Q30 rate: {data['summary']['after_filtering']['q30_rate']:.3f}\n\n")
            
            out.write(f"Reads passed filter: {data['filtering_result']['passed_filter_reads']:,}\n")
            out.write(f"Reads with low quality: {data['filtering_result']['low_quality_reads']:,}\n")
            out.write(f"Reads too short: {data['filtering_result']['too_short_reads']:,}\n")

# Rule 5: FastQC on trimmed reads
rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_1_trimmed.fastq.gz",
        "results/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        "results/fastqc/trimmed/{sample}_1_trimmed_fastqc.html",
        "results/fastqc/trimmed/{sample}_2_trimmed_fastqc.html",
        "results/fastqc/trimmed/{sample}_1_trimmed_fastqc.zip",
        "results/fastqc/trimmed/{sample}_2_trimmed_fastqc.zip"
    params:
        outdir="results/fastqc/trimmed"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input}
        """

# Rule 6: MultiQC for trimmed reads
rule multiqc_trimmed:
    input:
        expand("results/fastqc/trimmed/{{sample}}_{read}_trimmed_fastqc.zip", read=["1", "2"]),
        "results/trimmed/{sample}_fastp.json"
    output:
        "results/fastqc/trimmed/{sample}_multiqc_report.html"
    params:
        outdir="results/fastqc/trimmed",
        indir="results/trimmed"
    threads: 2
    shell:
        """
        multiqc -f -o {params.outdir} -n {wildcards.sample}_multiqc_report.html \
            {params.outdir}/*{wildcards.sample}* {params.indir}/*{wildcards.sample}*
        """

########### APPROACH 1: Kraken2 on reads ###############
rule kraken2_reads:
    input:
        r1="results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2="results/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        report="results/kraken2/reads/{sample}_kraken_report.txt",
        output="results/kraken2/reads/{sample}_kraken_output.txt",
        classified_1="results/kraken2/reads/{sample}_classified_1.fq",
        classified_2="results/kraken2/reads/{sample}_classified_2.fq",
        unclassified_1="results/kraken2/reads/{sample}_unclassified_1.fq",
        unclassified_2="results/kraken2/reads/{sample}_unclassified_2.fq"
    params:
        db=KRAKEN_DB,
        confidence=KRAKEN_CONFIDENCE,
        classified_base="results/kraken2/reads/{sample}_classified#.fq",
        unclassified_base="results/kraken2/reads/{sample}_unclassified#.fq"
    threads: THREADS
    log:
        "logs/kraken2/reads_{sample}.log"
    shell:
        """
        mkdir -p results/kraken2/reads
        
        kraken2 --db {params.db} \
            --threads {threads} \
            --confidence {params.confidence} \
            --report {output.report} \
            --output {output.output} \
            --classified-out {params.classified_base} \
            --unclassified-out {params.unclassified_base} \
            --paired \
            --gzip-compressed \
            {input.r1} {input.r2} 2> {log}
        """

# Rule 8: Bracken for abundance estimation at species level
rule bracken_species:
    input:
        report="results/kraken2/reads/{sample}_kraken_report.txt"
    output:
        "results/kraken2/reads/{sample}_bracken_species.txt"
    params:
        db=KRAKEN_DB,
        level="S",  # Species level
        threshold=10,
        read_len=400
    log:
        "logs/bracken/{sample}_species.log"
    shell:
        """
        bracken -d {params.db} \
            -i {input.report} \
            -o {output} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} 2> {log}
        """

# Rule 9: Bracken for abundance estimation at genus level
rule bracken_genus:
    input:
        report="results/kraken2/reads/{sample}_kraken_report.txt"
    output:
        "results/kraken2/reads/{sample}_bracken_genus.txt"
    params:
        db=KRAKEN_DB,
        level="G",  # genus level
        threshold=10,
        read_len=400
    log:
        "logs/bracken/{sample}_genus.log"
    shell:
        """
        bracken -d {params.db} \
            -i {input.report} \
            -o {output} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} 2> {log}
        """

############### APPROACH 2: Assembly-based analysis ###############
# Rule 10: Assembly with rnaSPAdes
rule rnaspades_assembly:
    input:
        r1="results/trimmed/{sample}_1_trimmed.fastq.gz",
        r2="results/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        transcripts="results/assembly/{sample}_transcripts.fasta",
        log="results/assembly/{sample}_spades.log"
    params:
        outdir="results/assembly/{sample}_spades_out",
        memory=100  # GB of RAM
    threads: THREADS
    log:
        "logs/assembly/{sample}_rnaspades.log"
    shell:
        """
        # Clean up any previous runs
        rm -rf {params.outdir}
        
        # Run rnaSPAdes
        rnaspades.py -1 {input.r1} -2 {input.r2} \
            -t {threads} \
            -m {params.memory} \
            -o {params.outdir} \
            --phred-offset 33 \
            2> {log}
        
        # Copy the transcripts
        if [ -f {params.outdir}/transcripts.fasta ]; then
            cp {params.outdir}/transcripts.fasta {output.transcripts}
        else
            # Fallback to hard_filtered_transcripts if main file doesn't exist
            cp {params.outdir}/hard_filtered_transcripts.fasta {output.transcripts}
        fi
        
        cp {params.outdir}/spades.log {output.log}
        """

# Rule 11: Assembly QC with QUAST
rule quast:
    input:
        "results/assembly/{sample}_transcripts.fasta"
    output:
        "results/assembly/quast/{sample}_report.html"
    params:
        outdir="results/assembly/quast/{sample}",
        min_contig=MIN_CONTIG_LENGTH
    threads: THREADS
    log:
        "logs/quast/{sample}.log"
    shell:
        """
        quast.py {input} \
            -o {params.outdir} \
            -t {threads} \
            --min-contig {params.min_contig} \
            --no-plots \
            2> {log}
        
        # Copy report to expected location
        cp {params.outdir}/report.html {output}
        """

# Rule 12: Kraken2 on assembled contigs
rule kraken2_contigs:
    input:
        contigs="results/assembly/{sample}_transcripts.fasta"
    output:
        report="results/kraken2/contigs/{sample}_kraken_report.txt",
        output="results/kraken2/contigs/{sample}_kraken_output.txt"
    params:
        db=KRAKEN_DB,
        confidence=KRAKEN_CONFIDENCE
    threads: THREADS
    log:
        "logs/kraken2/{sample}_contigs.log"
    shell:
        """
        mkdir -p results/kraken2/contigs
        
        kraken2 --db {params.db} \
            --threads {threads} \
            --confidence {params.confidence} \
            --report {output.report} \
            --output {output.output} \
            --use-names \
            {input.contigs} 2> {log}
        """
# Rule 13: Krona visualization for reads
rule krona_reads:
    input:
        "results/kraken2/reads/{sample}_kraken_report.txt"
    output:
        "results/kraken2/reads/{sample}_krona.html"
    log:
        "logs/krona/{sample}_reads.log"
    shell:
        """
        # Convert Kraken report to Krona format
        python -c "
import sys
with open('{input}', 'r') as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 6:
            count = parts[2]
            tax_path = parts[5].strip().replace(' ', '_')
            if tax_path and count != '0':
                print(f'{{count}}\\t{{tax_path}}')
        " > results/kraken2/reads/{wildcards.sample}_krona.txt
        
        ktImportText results/kraken2/reads/{wildcards.sample}_krona.txt -o {output} 2> {log}
        rm results/kraken2/reads/{wildcards.sample}_krona.txt
        """

# Rule 14: Krona visualization for contigs
rule krona_contigs:
    input:
        "results/kraken2/contigs/{sample}_kraken_report.txt"
    output:
        "results/kraken2/contigs/{sample}_krona.html"
    log:
        "logs/krona/{sample}_contigs.log"
    shell:
        """
        # Convert Kraken report to Krona format
        python -c "
import sys
with open('{input}', 'r') as f:
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 6:
            count = parts[2]
            tax_path = parts[5].strip().replace(' ', '_')
            if tax_path and count != '0':
                print(f'{{count}}\\t{{tax_path}}')
        " > results/kraken2/contigs/{wildcards.sample}_krona.txt
        
        ktImportText results/kraken2/contigs/{wildcards.sample}_krona.txt -o {output} 2> {log}
        rm results/kraken2/contigs/{wildcards.sample}_krona.txt
        """

# Rule 15: Generate summary statistics
rule summary_stats:
    input:
        raw_r1=f"{DATA_DIR}/{{sample}}_1.fastq.gz",
        trimmed_json="results/trimmed/{sample}_fastp.json",
        kraken_reads="results/kraken2/reads/{sample}_kraken_report.txt",
        kraken_contigs="results/kraken2/contigs/{sample}_kraken_report.txt",
        assembly="results/assembly/{sample}_transcripts.fasta"
    output:
        "results/summary/{sample}_pipeline_stats.txt"
    run:
        import json
        import gzip
        from pathlib import Path
        
        # Read fastp stats
        with open(input.trimmed_json, 'r') as f:
            fastp_data = json.load(f)
        
        # Count assembled transcripts
        transcript_count = 0
        total_length = 0
        with open(input.assembly, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    transcript_count += 1
                else:
                    total_length += len(line.strip())
        
        # Parse Kraken reports
        def parse_kraken_report(file_path):
            classified = 0
            unclassified = 0
            with open(file_path, 'r') as f:
                for line in f:
                    if 'unclassified' in line.lower():
                        parts = line.strip().split('\t')
                        unclassified = int(parts[2])
                    elif 'root' in line.lower() and '\t1\t' in line:
                        parts = line.strip().split('\t')
                        classified = int(parts[2])
            return classified, unclassified
        
        reads_classified, reads_unclassified = parse_kraken_report(input.kraken_reads)
        
        # Write summary
        with open(output[0], 'w') as out:
            out.write('=' * 60 + '\n')
            out.write('METAGENOMIC ANALYSIS PIPELINE SUMMARY\n')
            out.write(f'Sample: {wildcards.sample}\n')
            out.write('=' * 60 + '\n\n')
            
            out.write('INPUT DATA:\n')
            out.write('-' * 40 + '\n')
            out.write('Library Type: SsTS RNA\n')
            out.write('Technology: Illumina NovaSeq\n')
            out.write('Sequence mode: NovaSeq PE 150\n')
            out.write(f"Raw reads: {fastp_data['summary']['before_filtering']['total_reads']:,}\n")
            out.write(f"Raw bases: {fastp_data['summary']['before_filtering']['total_bases']:,}\n\n")
            
            out.write('QUALITY CONTROL:\n')
            out.write('-' * 40 + '\n')
            out.write(f"Reads after filtering: {fastp_data['summary']['after_filtering']['total_reads']:,}\n")
            out.write(f"Bases after filtering: {fastp_data['summary']['after_filtering']['total_bases']:,}\n")
            out.write(f"Q30 rate: {fastp_data['summary']['after_filtering']['q30_rate']:.1%}\n\n")
            
            out.write('ASSEMBLY STATISTICS:\n')
            out.write('-' * 40 + '\n')
            out.write(f'Assembled transcripts: {transcript_count:,}\n')
            out.write(f'Total assembly length: {total_length:,} bp\n')
            if transcript_count > 0:
                out.write(f'Average transcript length: {total_length/transcript_count:.1f} bp\n\n')
            
            out.write('TAXONOMIC CLASSIFICATION:\n')
            out.write('-' * 40 + '\n')
            out.write('Approach 1 (Read-based):\n')
            total_reads = reads_classified + reads_unclassified
            if total_reads > 0:
                out.write(f'  Classified reads: {reads_classified:,} ({reads_classified/total_reads:.1%})\n')
                out.write(f'  Unclassified reads: {reads_unclassified:,} ({reads_unclassified/total_reads:.1%})\n\n')
            
            out.write('Approach 2 (Assembly-based):\n')
            out.write(f'  See results/kraken2/contigs/{wildcards.sample}_kraken_report.txt\n\n')
            
            out.write('=' * 60 + '\n')

# Optional: Clean intermediate files
rule clean:
    shell:
        """
        echo "Cleaning intermediate files..."
        rm -rf results/assembly/*/spades_out
        rm -f results/kraken2/reads/*_classified*.fq
        rm -f results/kraken2/reads/*_unclassified*.fq
        echo "Clean complete!"
        """

# Rule to print workflow information
rule info:
    run:
        print("==========================================")
        print("RNA METAGENOMIC ANALYSIS PIPELINE")
        print("==========================================")
        print(f"Sample: {SAMPLE}")
        print("Input files:")
        print(f"  - {DATA_DIR}/{SAMPLE}_1.fastq.gz")
        print(f"  - {DATA_DIR}/{SAMPLE}_2.fastq.gz")
        print(f"Kraken2 database: {KRAKEN_DB}")
        print(f"Threads: {THREADS}")
        print("==========================================")
```
Now, we can run the snakemake pipeline like this:
```bash
snakemake \
    --cores 36 \
    --rerun-incomplete \
    --latency-wait 60 \
    --verbose \
    --stats snakemake_stats.txt \
    2>&1 | tee snakemake_log.txt
```
But never really run it on the head node, HPC Cluster and the administrator(s) won't be happy. Take an interactive session (using `srun` command) or submit a sbatch jub (using SLURM scheduler). 

Run it wrapping like this:
```bash
sbatch -p fat -c 36 --mem=180G --time=25:00:00 --wrap "snakemake \
    --cores 36 \
    --rerun-incomplete \
    --latency-wait 60 \
    --verbose \
    --stats snakemake_stats.txt \
    2>&1 | tee novaseq_metagenome_log.txt
```
Some HPC Cluster could have PBS instead of SLURM. Know yours. Mine is SLURM.

Instead of wrapping, we can write a proper submission script. Here goes the submission script, named `submit_pipeline.sh`:
```bash
#!/bin/bash
#SBATCH --job-name=rna_metagenomics
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180G
#SBATCH --time=12:00:00

# conda activate metagenomics_env

# Print job information
echo "Job started on $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# Set up temporary directory for SPAdes (uses local scratch for speed)
export TMPDIR=/scratch/$USER/$SLURM_JOB_ID
mkdir -p $TMPDIR

# Run Snakemake with cluster configuration
snakemake \
    --cores 36 \
    --rerun-incomplete \
    --latency-wait 60 \
    --verbose \
    --stats snakemake_stats.txt \
    2>&1 | tee novaseq_metagenome_log.txt

# Check exit status
if [ $? -eq 0 ]; then
    echo "Pipeline completed successfully!"
else
    echo "Pipeline failed with exit code $?"
fi

# Clean up temporary directory
rm -rf $TMPDIR

echo "Job finished on $(date)"
```
You see, I did not mention the snakefile. Snakemake got it by default. 

Usage of a sbatch/submission script is better. We can have error and log file for the whole run and it helps debugging.

By the way, we just wrote the script. Need to submit now. Write `submit_pipeline.sh` to submit the job. You will see the job submission confirmation message with a job id. 

We can use the job id to check status. Run `squeue <job_id>`
When the job is finished, we can check how much resources it used by running `sacct <job_id>`. So, we can decide later how much resource to ask and for how long.
