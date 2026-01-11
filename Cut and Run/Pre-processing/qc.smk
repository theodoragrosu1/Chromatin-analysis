#qc.smk

##sequencing QC
rule fastqc:
    input:
        "fastq/{sample}_{read}.fq.gz"
    output:
        "fastq/fastqc/{sample}_{read}_fastqc.html",
        "fastq/fastqc/{sample}_{read}_fastqc.zip"
    shell:
        "fastqc -o fastq/fastqc/ -f fastq {input}"

rule multiqc:
    input:
        expand("fastq/fastqc/{sample}_{read}_fastqc.html", sample = SAMPLE_NAMES, read = ["R1", "R2"]),
        expand("fastq/fastqc/{sample}_{read}_fastqc.zip", sample = SAMPLE_NAMES, read = ["R1", "R2"])
    output:
        "multiqc_report.html"
    shell:
        "multiqc ."


rule count_reads_fastq_before_trimming: 
    input:
        input_files = "fastq/{sample}_{read}.fq.gz",
        check_files_fastqc = rules.fastqc.output
    output:
        "count_reads/before_trimming/{sample}_{read}_reads.txt"
    shell:
        "echo $(zcat {input.input_files} |wc -l)/4|bc > {output}"



rule trimming:
	input:
		read1 = "fastq/{sample}_R1.fq.gz",
		read2 = "fastq/{sample}_R2.fq.gz"
	output:
		paired_read1 = "fastq_trimmed/{sample}_R1_paired.fq.gz",
		unpaired_read1 = "fastq_trimmed/{sample}_R1_unpaired.fq.gz",
		paired_read2 = "fastq_trimmed/{sample}_R2_paired.fq.gz",
		unpaired_read2 = "fastq_trimmed/{sample}_R2_unpaired.fq.gz"
	shell:
		"java -jar ./Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 " #fastp
		"{input.read1} {input.read2} "
		"{output.paired_read1} {output.unpaired_read1} {output.paired_read2} {output.unpaired_read2} "
		"ILLUMINACLIP:./Trimmomatic-0.39/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:15:4:1:true MINLEN:20"

rule fastqc_post_trimming:
    input:
        input_files = "fastq_trimmed/{sample}_{read}_paired.fq.gz",
        check_files_present1 = rules.trimming.output.paired_read1, #these two lines check if the trimmed files are present before running the rule 
        check_files_present2 = rules.trimming.output.paired_read2
    output:
        "fastq_trimmed/fastqc_post_trimming/{sample}_{read}_paired_fastqc.html",
        "fastq_trimmed/fastqc_post_trimming/{sample}_{read}_paired_fastqc.zip"
    shell:
        "fastqc -o fastq_trimmed/fastqc_post_trimming/ -f fastq {input.input_files}"#

rule count_reads_fastq_after_trimming: 
	input:
		input_files = "fastq_trimmed/{sample}_{read}_paired.fq.gz",
		check_files_fastqc = rules.fastqc_post_trimming.output
	output:
		"count_reads/before_alignment/{sample}_{read}.txt"
	shell:
		"echo $(zcat {input.input_files} |wc -l)/4|bc > {output}"


rule multiqc_post_trim:
    input:
        expand("fastq_trimmed/fastqc_post_trimming/{sample}_{read}_paired_fastqc.html", sample = SAMPLE_NAMES, read = ["R1", "R2"]),
        expand("fastq_trimmed/fastqc_post_trimming/{sample}_{read}_paired_fastqc.zip", sample = SAMPLE_NAMES, read = ["R1", "R2"])
    output:
        "multiqc_report_post_trimming/multiqc_report.html"
    shell:
        "multiqc . -o multiqc_report_post_trimming"
