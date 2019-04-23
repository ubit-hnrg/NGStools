#!/bin/bash
vcf_in=$1
out_prefix=$2

dbpath=/home/hnrg/HNRG-pipeline-V0.1/dbs/hg19_annovar/
annovar=/home/bitgenia/dbs/annovar/table_annovar.pl
perl $annovar $vcf_in $dbpath -vcfinput  -buildver hg19 -remove -out $out_prefix -protocol dbnsfp35a -operation f -nastring .
