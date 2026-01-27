#!/bin/bash
#make sure you are in the correct folder before starting this


#checks quality of fastq files
fastqc ./*/*.gz -t 4 #number of threads it's using, can be changed
echo "fastqc finished"

#trimmomatic step: 
master_directory="./raw_data"

#Loop through each pair of input files
for file_1 in ./*/*_1.fq.gz; do
    file_2="${file_1/_1/_2}"
    file_name=$(basename "$file_1")
    output_directory=$(echo "$file_name" | cut -d "_" -f 1-2)
    sample_name="${file_name%_1.fq.gz}"
    
    
    # RuntrimmomaticPE
    TrimmomaticPE -phred33 -threads 8 "$file_1" "$file_2" \
        "$master_directory/$output_directory/$sample_name"_1_paired.fq.gz \
        "$master_directory/$output_directory/$sample_name"_1_unpaired.fq.gz \
        "$master_directory/$output_directory/$sample_name"_2_paired.fq.gz \
        "$master_directory/$output_directory/$sample_name"_2_unpaired.fq.gz \
        ILLUMINACLIP:NexteraPE-PE.fa:1:30:4:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 #adapt as required - these should be good for ATAC-Seq with nextera adapters
done


echo "Trimmomatic finished"


#remove unnecesary files
rm ./*/*_L1_1.fq.gz
rm ./*/*_L1_2.fq.gz
rm ./*/*_L2_1.fq.gz
rm ./*/*_L2_2.fq.gz

#checks quality of trimmed adapters fastq files
fastqc ./*/*_paired.fq.gz -t 4 
echo "fastqc part2 finished"

#Alignment with bowtie

#bowtie2-build ./references/hg38.fa GRCh38
mkdir ./aligned

for file_1 in ./*/*_1_paired.fq.gz; do
    file_2="${file_1/_1_paired/_2_paired}"
    file_name=$(basename "$file_1")
    output_name="${file_name%_1_paired.fq.gz}"
    output_directory="./aligned"

    bowtie2 -p 4 -X 2000 --very-sensitive -t -x ./references/genome -1 "$file_1" -2 "$file_2" | samtools view -bS - | samtools sort -n > "$output_directory/$output_name".bam

done

echo "bowtie2 alignment finished"

mkdir ./peaks
peaks_dir="./peaks"

cd ./aligned


#calling peaks using Genrich
for bam_file in *.bam; do 
    file_name=$(basename "$bam_file")
    output_name="${file_name%.bam}"
    ./Genrich/Genrich -e chrM -r -j -q 0.05 -t "$bam_file" -o "$peaks_dir/$output_name".narrowpeak -k "$peaks_dir/$output_name".bedGraph #genrich installed in this directory
done

echo "peak calling finished"

cd ../peaks

#getting additional files (bedGraph)
for bedgraph_file in *.bedGraph; do 
    file_name=$(basename "$bedgraph_file")
    output_name="${file_name%.bedGraph}"
    cat "$bedgraph_file" | tail -n +3 | cut -f 1-4 | sort -k1,1 -k2,2n > "$output_name"_cut.bedGraph
done

echo "bedgraph manipulation finished"


#create index for the genome file
samtools faidx ./references/genome.fa

#get chromosome sizes file
cut -f 1,2 ./references/genome.fa.fai > ./references/chrom.sizes


#activate virtual environment that has installed the bedGraphToBigWig function
conda activate bigwig_processing

#getting bigwig wiles
for cut_bedgraph_file in *_cut.bedGraph; do 
    file_name=$(basename "$cut_bedgraph_file")
    output_name="${file_name%.bedGraph}"
    bedGraphToBigWig "$cut_bedgraph_file" ./references/chrom.sizes "$output_name".bw 

done

echo "getting bigwig files finished - yay all done ^_^"