#!/bin/bash

# Parent directory to search for subdirectories with fastq files
echo "starting"
parent_dir1="./RawData_october"
parent_dir2="./RawData_november"

# Loop through each subdirectory in the parent directory
for dir in "$parent_dir1"/*; do
  
  # Get the name of the subdirectory
  echo $(basename "$dir")
  subdir=$(basename "$dir")

  echo "$parent_dir1/$subdir"
  echo "$parent_dir2/$subdir"

  # Merge all files with the suffix "_1.fq.gz"
  cat "$parent_dir1/$subdir"/*_1.fq.gz "$parent_dir2/$subdir"/*_1.fq.gz > ./merged/"${subdir}_R1.fq.gz"

  # Merge all files with the suffix "_2.fq.gz"
  cat "$parent_dir1/$subdir"/*_2.fq.gz "$parent_dir2/$subdir"/*_2.fq.gz > ./merged/"${subdir}_R2.fq.gz"

  echo "merged files created"
done