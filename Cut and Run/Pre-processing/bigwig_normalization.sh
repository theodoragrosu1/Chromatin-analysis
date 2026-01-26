#!/bin/bash

#this script normalizes all the merged bigwigs (generated with bamCov from merged bam files) to their respective IgG for log fold change visualization 

bigwigCompare -b1 WT_me3.bigwig -b2 WT_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/WT_me3_log2.bigwig

bigwigCompare -b1 WT_ac.bigwig -b2 WT_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/WT_ac_log2.bigwig

bigwigCompare -b1 WT_ub.bigwig -b2 WT_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/WT_ub_log2.bigwig

bigwigCompare -b1 A7_me3.bigwig -b2 A7_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/A7_me3_log2.bigwig

bigwigCompare -b1 A7_ac.bigwig -b2 A7_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/A7_ac_log2.bigwig

bigwigCompare -b1 A7_ub.bigwig -b2 A7_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/A7_ub_log2.bigwig

bigwigCompare -b1 C9_me3.bigwig -b2 C9_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/C9_me3_log2.bigwig

bigwigCompare -b1 C9_ac.bigwig -b2 C9_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/C9_ac_log2.bigwig

bigwigCompare -b1 C9_ub.bigwig -b2 C9_IgG.bigwig -bl ../public_datasets/genome_files/hg38-blacklist.v2.bed -p max/2 -o ./normalized_to_IgG/C9_ub_log2.bigwig