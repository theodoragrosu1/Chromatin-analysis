#filtering.smk

rule filter_duplicates:
	input:
		input_files = "dupes/{sample}/{sample}_marked_duplicates.bam",
		index_files = "dupes/{sample}/{sample}_marked_duplicates.bam.bai"
	output:
		"deduped/{sample}/{sample}_deduped.bam"
	shell:
		"samtools view -F 1280 -f 2 -b {input.input_files} > {output}" #remove multimappers -F 1280 (multimappers & pcr/optical duplicates) ###-F 3328 remove supplementary alignment as well



#rule remove_blacklisted_reg:
#	input: 
#		rules.filter_duplicates.output
#	output:
#		"deduped/{sample}/{sample}_deduped_blacklist_filtered.bam"
#	shell:
#		"bedtools intersect -v -abam {input} -b ./hg38-blacklist.v2.bed > {output}"


rule index_bams:
	input:
		rules.filter_duplicates.output
	output:
		"deduped/{sample}/{sample}_deduped.bam.bai"
	shell:
		"samtools index {input}"


rule count_reads_after_removing_dupes: 
	input:
		input_files = rules.filter_duplicates.output,
		wait_for_finish = rules.index_bams.output
	output:
		"count_reads/after_alignment/{sample}_reads_deduped_bam.txt"
	shell:
		"bamtools stats -in {input.input_files} > {output}"
