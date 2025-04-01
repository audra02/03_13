import glob
import os

# Input path
input_path = "/home/shared/raw_data/ATACseq_example"

# Get the list of sample names
files_R1 = glob.glob(os.path.join(input_path, "*_R1.fastq.gz"))
samples = sorted(set(os.path.basename(f).replace("_R1.fastq.gz", "") for f in files_R1))

print(samples)

# Define the all rule: The final target files
rule all:
    input:
        expand("results/fastqc/{sample}_{suffix}_fastqc.html", sample=samples, suffix=["R1", "R2", "R1_trimmed", "R2_trimmed"]),
        expand("results/trimmed/{sample}_R1_trimmed.fq.gz", sample=samples),
        expand("results/trimmed/{sample}_R2_trimmed.fq.gz", sample=samples),
        expand("results/mapped/{sample}.bam", sample=samples),
        "NC_000022.fasta",  # Ensure chr22 is downloaded before indexing
        expand("results/mapped/{sample}_filtered.bam", sample=samples),
        expand("results/mapped/{sample}_no_duplicates.bam", sample=samples),  # Added no duplicates BAM
        expand("results/mapped/{sample}_duplication_metrics.txt", sample=samples),  # Added duplication metrics
        expand("results/metrics/{sample}_insert_size_metrics.txt", sample=samples),  # Added insert size metrics
        expand("results/metrics/{sample}_insert_size_histogram.pdf", sample=samples)  # Added insert size histogram

rule fastqc:
    input:
        fq1 = lambda wildcards: f"{input_path}/{wildcards.sample}_R1.fastq.gz",
        fq2 = lambda wildcards: f"{input_path}/{wildcards.sample}_R2.fastq.gz",
        fq3 = "results/trimmed/{sample}_R1_trimmed.fq.gz",
        fq4 = "results/trimmed/{sample}_R2_trimmed.fq.gz"  
    output:
        "results/fastqc/{sample}_R1_fastqc.html",
        "results/fastqc/{sample}_R2_fastqc.html",
        "results/fastqc/{sample}_R1_trimmed_fastqc.html",
        "results/fastqc/{sample}_R2_trimmed_fastqc.html",
        "results/fastqc/{sample}_R1_fastqc.zip",
        "results/fastqc/{sample}_R2_fastqc.zip",
        "results/fastqc/{sample}_R1_trimmed_fastqc.zip",
        "results/fastqc/{sample}_R2_trimmed_fastqc.zip"
    conda: "envs/preprocess_rnaseq.yaml"
    threads: 4
    shell:
        """
        fastqc -t {threads} --outdir results/fastqc {input.fq1} {input.fq2} {input.fq3} {input.fq4}
        """


rule cutadapt:
    input:
        fq1 = lambda wildcards: f"{input_path}/{wildcards.sample}_R1.fastq.gz",
        fq2 = lambda wildcards: f"{input_path}/{wildcards.sample}_R2.fastq.gz",
    output:
        fq1 = "results/trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = "results/trimmed/{sample}_R2_trimmed.fq.gz"
    conda: "envs/preprocess_rnaseq.yaml"
    threads: 4
    shell:
        """
        cutadapt -j {threads} \
        -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        -m 20 -q 20  -o {output.fq1} -p {output.fq2} \
        {input.fq1} {input.fq2}
        """

rule fetch_chr22:
    output:
        "NC_000022.fasta"
    shell:
        """
        efetch -db nucleotide -id NC_000022 -format fasta > {output}
        """

rule bowtie2_index:
    input:
        "NC_000022.fasta"
    output:
        expand("chr22_index.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"])
    shell:
        """
        bowtie2-build {input} chr22_index
        """
        
rule bowtie2_map:
    input:
        index = expand("chr22_index.{ext}", ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]),
        fq1 = "results/trimmed/{sample}_R1_trimmed.fq.gz",
        fq2 = "results/trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        "results/mapped/{sample}.sam"
    threads: 8
    shell:
        """
        bowtie2 --very-sensitive --end-to-end --dovetail \
        --maxins 1000 -x chr22_index -1 {input.fq1} -2 {input.fq2} \
        -S {output}
        """

rule sam_to_bam:
    input:
        "results/mapped/{sample}.sam"
    output:
        "results/mapped/{sample}.bam",
        "results/mapped/{sample}.bam.bai"
    shell:
        """
        samtools view -@ 8 -bS {input} | samtools sort -@ 8 -o {output[0]}
        samtools index {output[0]}
        """

rule filter_bam:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}_filtered.bam",
        "results/mapped/{sample}_filtering_stats.txt"
    threads: 4
    shell:
        """
        original_count=$(samtools view -c {input})
        bamtools filter -in {input} -out {output[0]} -mapQuality ">=30" -isProperPair true
        filtered_count=$(samtools view -c {output[0]})
        echo "Original Reads: $original_count" > {output[1]}
        echo "Filtered Reads: $filtered_count" >> {output[1]}
        echo "Removed Reads: $(($original_count - $filtered_count))" >> {output[1]}
        """
rule picard_mark_duplicates:
    input:
        bam = "results/mapped/{sample}_filtered.bam"
    output:
        bam_no_duplicates = "results/mapped/{sample}_no_duplicates.bam",
        bam_metrics = "results/mapped/{sample}_duplication_metrics.txt"
    conda: "envs/preprocess_rnaseq.yaml"
    shell:
        """
        picard MarkDuplicates \
        I={input.bam} \
        O={output.bam_no_duplicates} \
        M={output.bam_metrics} \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT
        """

rule collect_insert_size_metrics:
    input:
        bam = "results/mapped/{sample}_no_duplicates.bam"
    output:
        insert_size_metrics = "results/metrics/{sample}_insert_size_metrics.txt",
        insert_size_histogram = "results/metrics/{sample}_insert_size_histogram.pdf"
    shell:
        """
        picard CollectInsertSizeMetrics \
        I={input.bam} \
        O={output.insert_size_metrics} \
        H={output.insert_size_histogram} \
        M=0.5
        """


