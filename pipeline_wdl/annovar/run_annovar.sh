#!/bin/bash
vcf_in=$1
out_prefix=$2
dbpath=/home/hnrg/HNRG-pipeline-V0.1/dbs/hg19_annovar/
annovar=/home/bitgenia/dbs/annovar/table_annovar.pl
perl $perl $vcf_in $dbpath -vcfinput  -buildver hg19 -remove -out $out_prefix -protocol refGene,avsnp150,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_exome,gnomad_genome,clinvar_20180603,intervar_20180118,dbscsnv11,dbnsfp35a,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g -nastring . -otherinfo
