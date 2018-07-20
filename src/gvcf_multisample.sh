#!/bin/bash

# this script generate a  multisample vcf file by analyzing gvcf files of each sample   
# the output is called      
# test_output_genotypeGVCFs.vcf

# this must be runned calling a docker with 

# para montar el docker de gatk sobre un directorio existente

sudo docker run -v ~/data_bundle_wdl/WDLdata:/gatk/my_data/ -it broadinstitute/gatk




#father
gatk --java-options "-Xmx4g" HaplotypeCaller -R my_data/ref/ref.fasta -I my_data/bams/father.bam -O ./my_data/father.g.vcf \
-ERC GVCF 

#mother
gatk --java-options "-Xmx4g" HaplotypeCaller -R my_data/ref/ref.fasta -I my_data/bams/mother.bam -O ./my_data/mother.g.vcf \
-ERC GVCF 

#son
gatk --java-options "-Xmx4g" HaplotypeCaller -R my_data/ref/ref.fasta -I my_data/bams/son.bam -O ./my_data/son.g.vcf \
-ERC GVCF


############### ES MUUUY LENTOOOO!!! ####
#gatk CombineGVCFs \
#   -R my_data/ref/ref.fasta \
#   --variant my_data/son.g.vcf \
#   --variant my_data/father.g.vcf \
#   --variant my_data/mother.g.vcf \
#   -O ./my_data/trio.g.vcf


#esto va como piña muy rapido y mismo output
gatk GenomicsDBImport \
    -V my_data/mother.g.vcf \
    -V my_data/father.g.vcf \
    -V my_data/son.g.vcf \
    --genomicsdb-workspace-path my_data/trio_database \
    -L 20

# el output de esto es una bosta, necesitamos para poder leerlo hacer:

 gatk SelectVariants \
    -R my_data/ref/ref.fasta \
    -V gendb://my_data/trio_database \
    -O my_data/combined.g.vcf   

###################
## esto anda, pero no entiendo la notacion de gennotipos

 gatk GenotypeGVCFs \
    -R my_data/ref/ref.fasta \
    -V gendb://my_data/trio_database \
    -O my_data/test_output_genotypeGVCFs.vcf
