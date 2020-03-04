#!/bin/bash


#comparar 2 vcfs con isec 
set -e
vcf1=$1
vcf2=$2
output=$3

bcftools isec -p $output $1 $2
