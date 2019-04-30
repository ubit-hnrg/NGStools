#!/bin/bash
bamfile=$1
bedfile=$2;
intersectBed -a $bamfile -b $bedfile > $bedfile.bam
samtools index $bedfile.bam $bedfile.bai