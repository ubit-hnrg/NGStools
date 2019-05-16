#!/bin/bash
bamfile=$1
gen=$2;
sample=$3
zcat /home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_canonical_transcript_GRCh37.tsv.gz |grep -P $gen'\t'|head -n1 |cut -f1,2,3 > $gen.bed
intersectBed -a $bamfile -b $gen.bed > $sample'_'$gen.bam
samtools index $sample'_'$gen.bam $sample'_'$gen.bai