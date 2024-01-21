#!/bin/bash
bamfile=$1
bedfile=$2;
/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/intersectBed -a $bamfile -b $bedfile > $bedfile.bam
samtools index $bedfile.bam $bedfile.bai
