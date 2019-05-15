toolpath=/home/hnrg/HNRG-pipeline-V0.1/tools
ngstools_path=/home/hnrg/NGStools
runID="TSO20192803"
joint_vcf=/home/hnrg/resultsHNRG/$runID/TSO20192803.vcf

tso_vcf=/home/hnrg/resultsHNRG/$runID/$runID'_TSO.vcf'
prefix=$runID'_TSO';

tso_padded=/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/TruSight_One_v1_padded_100_GRCh37.bed

#output=one_smple_vcf

declare -a samples=("1711242" "1802140" "1804642" "1805817" "1809341" "1809720" "1810255" "1901364" "1901981" "EB761" "EB790" "EB802")


# 1rst step restrict VCF to regions in bedfile. 
cat $joint_vcf | java -jar $toolpath/SnpSift.jar intervals $tso_padded > $tso_vcf

# 2nd step prepare the Joint VCF file to be renamed and splited into multiple ones. 
#tso_renamed=$runID'_TSO_renamed.vcf'
#cp $tso_vcf $tso_renamed


#############################
#iterate across samples. 
#################################

for i in "${samples[@]}";
# rename numeric id samples
do echo $i;sed -i -e"s/$i/ID$i/g" $tso_vcf;

# prepare id  ## This is only for supportting samplenames starting with numbers. (SnpSift compatibility)
if [[ $i =~ ^1.* ]]
then
    id='ID'$i
else
    id=$i
fi
echo $id;

#deffine input/outpus
faceted_one_sample_vcf=/home/hnrg/resultsHNRG/$runID/$i/$i'_TSO_faceted_one_sample_SnpSift.vcf';
one_smple_vcf=/home/hnrg/resultsHNRG/$runID/$i/$i'_TSO_one_sample_SnpSift.vcf';
one_sample_prefix_path=/home/hnrg/resultsHNRG/$runID/$i/$i'_TSO_faceted_one_sample_SnpSift';
output=/home/hnrg/resultsHNRG/$runID/$i/$i'_TSO_one_sample_SnpSift.vcf';

one_sample_prefix_basename=$i'_TSO_one_sample_SnpSift'
multianno_multisample_tsv=/home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv


# split vcf according to $i sample. This file contain all samples but only those relevant for $i one
cat /home/hnrg/resultsHNRG/$runID/tso_vcf | java -jar /home/hnrg/HNRG-pipeline-V0.1/tools/SnpSift.jar filter "(GEN[$id].GT!='./.')&(GEN[$id].GT != '0/0')" > $faceted_one_sample_vcf
#these steps remove the remaining samples of the vcf.
cat <(grep '^##' $faceted_one_sample_vcf) <(grep -v '^##' $faceted_one_sample_vcf| csvcut -t -c '#CHROM',POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,$id | csvformat -T) > $one_sample_vcf
rm $faceted_one_sample_vcf


#run annovar over this file.
$ngstools_path/pipeline_wdl/annovar/run_annovar.sh $one_smple_vcf $one_sample_prefix_path

# prepare mjultianno.tsv from multianno.txt for deliver (do not forget postprocess InterVar)
#1) le saco las columnas 163-165. OJO CON ESTO QUE 
#2) modifico el header agregando las columnas que faltan del vcf (remplazo x otherinfo)
#3) joineo al archivo multinanno las columnas de genotipo de las restantes muestras de  la corrida.

    #columnas a cortar (localizando Otherinfo column y las 2 siguientes)
    nl0=$(head -n1 $one_sample_prefix_path.hg19_multianno.txt|tr '\t' '\n'|nl|grep 'Otherinfo'|cut -f1)
    nl1=$((nl0 + 1))
    nl2=$((nl0 + 2))

    # meto header (dejando el campo 'Otherinfo' que despues va a aser remplazado por las columnas del vcf original)
    head -n1 $one_sample_prefix_path.hg19_multianno.txt > $one_sample_prefix_path.hg19_multianno.tsv
    tail -n+2 $one_sample_prefix_path.hg19_multianno.txt|cut -f$nl0,$nl1,$nl2 --complement >>  $one_sample_prefix.hg19_multianno.tsv;
    vcf_header=$(grep '#CH' $one_sample_vcf);

    #remplazo columnas
    sed -i "s/Otherinfo/$vcf_header/g" $one_sample_prefix.hg19_multianno.tsv;

    #join one multianno tsv file AND joint genotyped vcf. This script (join_vcf.py) also postprocess Intervar columns.
    python /home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py --multianno_tsv=$one_sample_prefix_path.hg19_multianno.tsv --vcf_multisample=$tso_vcf --output=$multianno_multisample_tsv
    #change dots by tabs.
    sed -i -e "s|\.	|	|g" $multianno_multisample_tsv



    # take annovar.multianno.vcf and annotate this for bplaform usage .
    # Intervar postprocessing is written inside wdl annotation workflow

    path='/home/hnrg/resultsHNRG/'$runID/$i
    inputvcf=$one_sample_prefix_path.hg19_multianno.vcf
    mkdir -p inputs
    sed -e "s|__input__|$inputvcf|g" $ngstools_path/pipeline_wdl/anotacion_funcional/template_hnrg-anotacion_funcional_sin_CADD.json > ./inputs/$i'_hnrg-anotacion_funcional_sin_CADD.json'
    sed -i -e "s|__samplename__|$i|g" ./inputs/$i'_hnrg-anotacion_funcional_sin_CADD.json'
    sed -i -e "s|__path_save__|$path|g" ./inputs/$i'_hnrg-anotacion_funcional_sin_CADD.json'
    java -Dconfig.file=/home/hnrg/HNRG-pipeline-V0.1/src/wdl/cromwell.conf -jar $toolpath/cromwell-37.jar run -i ./inputs/$i'_hnrg-anotacion_funcional_sin_CADD.json' $ngstools_path/pipeline_wdl/anotacion_funcional/anotaciones_hnrg.wdl 



done
