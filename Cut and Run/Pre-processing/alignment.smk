#aligment.smk

rule align_reads:
	input:
		read1 = "fastq_trimmed/{sample}_R1_paired.fq.gz",
		read2 = "fastq_trimmed/{sample}_R2_paired.fq.gz",
	output:
		alignment = temp("sams/{sample}.sam"), #removes temp sam after the .bams are created
		stats = "stats_{sample}_bowtie2.txt"
	threads: 8
	params:
		genome_idx = "./GRCh38_genome" #path to genome sample
	shell:
		"bowtie2 --end-to-end --dovetail -I 10 -X 700 --no-mixed --no-discordant --threads {threads} -x {params.genome_idx} -1 {input.read1} -2 {input.read2} "
		"-S {output.alignment} &> {output.stats}"



rule sam_to_bam:
	input:
		rules.align_reads.output.alignment
	output:
		temp("bams/{sample}.bam")
	shell:
		"samtools view -bS {input} > {output}"


rule count_reads_after_alignment: 
	input:
		rules.sam_to_bam.output
	output:
		"count_reads/after_alignment/{sample}_reads_after_alignment.txt"
	shell:
		"bamtools stats -in {input} > {output}"

rule sort_bam:
    input:
        bam_files = rules.sam_to_bam.output
    output:
        temp("sorted_bams/{sample}/{sample}_sorted.bam")
    shell:
        "java -jar ./picard.jar SortSam I={input.bam_files} "
        "O={output} "
        "SORT_ORDER=coordinate "
        "VALIDATION_STRINGENCY=LENIENT"

rule add_read_groups:
    input:
        bam = rules.sort_bam.output
    output:
        bam_rg = "sorted_bams/{sample}/{sample}_sorted_rg.bam"
    params:
        rgid = lambda wildcards: wildcards.sample,
        rglb = "lib1",
        rgpl = "illumina",
        rgpu = "unit1",
        rgsm = lambda wildcards: wildcards.sample
    log:
        "logs/add_read_groups_{sample}.log"
    shell:
        """
        java -jar /mnt/Data/Bond/tgrosu/cutrun/picard.jar AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.bam_rg} \
            RGID={params.rgid} \
            RGLB={params.rglb} \
            RGPL={params.rgpl} \
            RGPU={params.rgpu} \
            RGSM={params.rgsm} \
            VALIDATION_STRINGENCY=LENIENT \
            &> {log}
        """


rule mark_duplicates:
	input:
		rules.add_read_groups.output
	output:
		dupes = "dupes/{sample}/{sample}_marked_duplicates.bam",
		metrics = "dupes/{sample}_dup_metrics.txt"
	shell:
		"java -jar ./picard.jar MarkDuplicates I={input} "
		"O={output.dupes} "
		"M={output.metrics} "
		"VALIDATION_STRINGENCY=LENIENT"


rule index_bams_dupes:
	input:
		rules.mark_duplicates.output.dupes
	output:
		"dupes/{sample}/{sample}_marked_duplicates.bam.bai"
	shell:
		"samtools index {input}"
