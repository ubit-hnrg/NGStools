#!/bin/bash

for f in locate *MAPQ0.bam;
do
basef= $(basename f)
echo $basef
done;
#tools/variantCallers/gatk-4.0.7.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R bundle/hgref/human_g1k_v37_decoy.fasta -I samples/exomasbitgenia/GENETIX/68641LSNDA_Result/68641LSNDA_MAPQ0.bam -O outputariel/test.g.vcf -ERC GVCF -L outputariel/querySNP.list
