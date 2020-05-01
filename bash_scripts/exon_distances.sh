#!/bin/bash
set -e

vcf=$1
vcf_name = $(basename $vcf .final_annot.HNRG_freq.vcf )
vcf_pwd = $(pwd $vcf)
exon_file=$2


grep "^#" $vcf > head_vcf.vcf
grep -v "#" $vcf| sort -k1,1 -k2,2n >> head_vcf.vcf

bedtools closest -a head_vcf.vcf -b $exon_file -D ref | cut -f 1,2,12-17  > $vcf_pwd/$vcf_name'.exon_distance.tsv'

rm head_vcf.vcf


