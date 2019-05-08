#!/bin/bash
# convert multianno into annovar vcf db (sorted)
vcf=/home/hnrg/resultsHNRG/1711242/1711242.final_annot.vcf
base=$(basename $vcf)
DIR=$(readlink -f $(dirname {$file}))
input=/home/hnrg/resultsHNRG/1711242/1711242_TSO_renamed_one_sample.hg19_multianno.txt
out=17111242_intervarDB.vcf
python NGStools/pipeline_wdl/anotacion_funcional/create_InterVarDB.py -i=$input -o $out

# index
bgzip $out
tabix $out.gz


#java -Xmx4G -jar /home/bitgenia/samples/HNRG-pipeline-V0.1/tools/SnpSift.jar annotate -v -info InterVarEvidence,InterVarVeredict $out.gz $vcf > $vcf'intervar.vcf'
/home/hnrg/HNRG-pipeline-V0.1/tools/vcfanno_linux64 -p 4 /home/hnrg/NGStools/pipeline_wdl/anotacion_funcional/intervar_anno_fromVCF.tom $vcf > $DIR/$BASE'_intervar.vcf'
