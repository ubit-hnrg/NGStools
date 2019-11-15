#!/bin/bash

bams=$1
bed=$2
fasta=$3
out=$4

# bed = /home/usuario/bam/intervalo_b37_padded_100.bed
# fasta = /home/usuario/reference/hs37d5.fa
###step1
sudo Rscript ReadInBams.R --bams $bams --bed $bed --fasta $fasta --out $out/resultado1 \ 

##step2
Rscript IdentifyFailures.R --Rdata $out/resultado1.RData --mincorr .98 --mincov 100 --out $out/resultado2 \

###step 3
Rscript makeCNVcalls.R --Rdata $out/resultado1.RData --transProb 0.01 --out $out/resultado3 --plot None \


&& while read line; do name=$(basename $line .bam); head -n1 $out/resultado3_all.txt > $out/$name.txt && grep "$name" $out/resultado3_all.txt >> $out/$name.txt; done < $bams