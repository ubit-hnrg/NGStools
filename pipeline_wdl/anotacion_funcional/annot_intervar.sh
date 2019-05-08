#!/bin/bash
# convert multianno into annovar vcf db (sorted)
vcfinput=$1
multianno=$2
vcfDB='intervar_sample_BD.vcf'


base=$(basename $vcfinput)
BASE="${BASE%.*}"
DIR=$(readlink -f $(dirname {$file}))
echo $DIR/$base'_intervar.vcf'

# multianno to vcf db file
python NGStools/pipeline_wdl/anotacion_funcional/create_InterVarDB.py -i=$multianno -o $vcfDB

# index
bgzip $vcfDB
tabix $vcfDB.gz

sed -e "s|__vcfDB__|$vcfDB.gz|g" /home/hnrg/NGStools/pipeline_wdl/anotacion_funcional/intervar_vcfanno_template_fromVCF.tom > config_vcfanno.tom
#java -Xmx4G -jar /home/bitgenia/samples/HNRG-pipeline-V0.1/tools/SnpSift.jar annotate -v -info InterVarEvidence,InterVarVeredict $vcfDB.gz $vcf > $vcf'intervar.vcf'
/home/hnrg/HNRG-pipeline-V0.1/tools/vcfanno_linux64 -p 4 config_vcfanno.tom $vcfinput > $DIR/$BASE'_intervar.vcf'
