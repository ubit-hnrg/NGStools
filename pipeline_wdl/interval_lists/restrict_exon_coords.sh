#!/bin/bash
bedfile=$1
exoncoords='/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_GRCh37.tsv'
intersectBed -wa -a $exoncoords -b $bedfile |sort -k1,1 -k2,2n -V |uniq 
