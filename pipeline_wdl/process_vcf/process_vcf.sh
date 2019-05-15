toolpath=/home/bitgenia/samples/HNRG-pipeline-V0.1/tools
joint_vcf=/home/hnrg/resultsHNRG/TSO20190328/TSO20192803.vcf
tso_vcf=/home/hnrg/resultsHNRG/TSO20190328/TSO20190328_TSO.vcf
tso_padded=/home/bitgenia/samples/HNRG-pipeline-V0.1/libraries/GRCh37/TruSight_One_v1_padded_100_GRCh37.bed

# 1rst step restrict VCF to regions in bedfile. 
cat $joint_vcf | java -jar $toolpath/SnpSift.jar intervals $tso_padded > $tso_vcf


# 2nd step rename and split Joint VCF into multiple ones. 
tso_renamed=TSO20190328_TSO_renamed.vcf
cp TSO20190328_TSO.vcf $tso_renamed


# rename numeric id samples
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981;
do echo $i;sed -i -e"s/$i/ID$i/g" $tso_renamed;
done


## use annovar for splitting multisample files. 
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
do
if [[ $i =~ ^1.* ]]
then
    id='ID'$i
else
    id=$i
fi
echo $id;
prefix='TSO20190328_TSO_renamed';
/home/bitgenia/dbs/annovar/convert2annovar.pl -format vcf4 TSO20190328_TSO_renamed.vcf -outfile $prefix -allsample -include -comment;
grep -P '^#' $prefix.$id.avinput > $prefix.$id.vcf;
grep -v -P '^#' $prefix.$id.avinput | cut -f 6- >>$prefix.$id.vcf;
mv $prefix.$id.vcf /home/hnrg/resultsHNRG/$i/;
sample_vcf=/home/hnrg/resultsHNRG/$i/$prefix.$id.vcf;
out_prefix=/home/hnrg/resultsHNRG/$i/$i'_TSO_renamed_one_sample';
~/NGStools/pipeline_wdl/annovar/run_annovar.sh $sample_vcf $out_prefix;
done

###   modifico el archivo multianno de annovar:
#1) le saco las columnas 163-165. OJO CON ESTO QUE 
#2) modifico el header agregando las columnas que faltan del vcf (remplazo x otherinfo)
#3) joineo al archivo multinanno las columnas de genotipo de las restantes muestras de  la corrida. 
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
    do
    if [[ $i =~ ^1.* ]]
    then
        id='ID'$i
    else
        id=$i
    fi
    echo $id;
    out_prefix=/home/hnrg/resultsHNRG/$i/$i'_TSO_renamed_one_sample';

    #columnas a cortar (localizando Otherinfo column y las 2 siguientes)
    nl0=$(head -n1 $out_prefix.hg19_multianno.txt|tr '\t' '\n'|nl|grep 'Otherinfo'|cut -f1)
    nl1=$((nl0 + 1))
    nl2=$((nl0 + 2))

    # meto header (dejando el campo 'Otherinfo' que despues va a aser remplazado por las columnas del vcf original)
    head -n1 $out_prefix.hg19_multianno.txt > $out_prefix.hg19_multianno.tsv
    tail -n+2 $out_prefix.hg19_multianno.txt|cut -f$nl0,$nl1,$nl2 --complement >>  $out_prefix.hg19_multianno.tsv;

    vcf_header=$(grep '#CH' /home/hnrg/resultsHNRG/$i/TSO20190328_TSO_renamed.$id.vcf);
    #remplazo columnas
    sed -i "s/Otherinfo/$vcf_header/g" $out_prefix.hg19_multianno.tsv;
    python /home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py --multianno_tsv=$out_prefix.hg19_multianno.tsv --vcf_multisample=/home/hnrg/resultsHNRG/TSO20190328/TSO20190328_TSO_renamed.vcf --output=/home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv
    sed -i -e "s|\. |   |g" /home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv
    done

#este vcf se puede usar para anotar. 
#El header debe skipear tres columnas de other info 163,164,y 165. 
#Luego las restantes (166 hasta el final) son las columnas del vcf de entrada.
# la ventaja de hacerlo así es que naturalmente quedan sólo los alternativos correspondientes a la muestra que se está mirando. 




## genero listas de genes
python ~/NGStools/pipeline_wdl/process_vcf/parse-sheets.py -i=/home/hnrg/metadataHNRG/TSO-20190328-bioinfo.xlsx -o=/home/hnrg/metadataHNRG/

for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
    do
    head -n1 /home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv > /home/hnrg/resultsHNRG/$i/$i'_geneList'.multiano_multisample.tsv
    grep /home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv -f /home/hnrg/metadataHNRG/$i/$i'_genelist.txt' >> /home/hnrg/resultsHNRG/$i/$i'_geneList'.multiano_multisample.tsv
    sed -i -e "s|\. |   |g" /home/hnrg/resultsHNRG/$i/$i'_geneList'.multiano_multisample.tsv
done


### merge results into a xlsx file
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
    do echo $i;
    /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py /home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv:Variants \
    /home/hnrg/resultsHNRG/$i/$i'_geneList'.multiano_multisample.tsv:GeneListVariants /home/hnrg/resultsHNRG/$i/$i'_coverage_statistics_by_exon.tsv':ExonCoverage \
    /home/hnrg/resultsHNRG/$i/$i.variants.xlsx
    done


## ONLY FOR CHECK. 
#check if gene List in our reference gene dabtase
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
    do echo $i;
    grep /home/hnrg/metadataHNRG/$i/$i'_genelist.txt' -v -f refGene.list >/home/hnrg/metadataHNRG/$i/$i'NotIn_refGene.txt';
    done




## use snpsift for splitting multisample files. 
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
do
if [[ $i =~ ^1.* ]]
then
    id='ID'$i
else
    id=$i
fi
echo $id;
prefix='TSO20190328_TSO_renamed';
output=/home/hnrg/resultsHNRG/$i/$i'_TSO_renamed_one_sample_SnpSift.vcf';

cat /home/hnrg/resultsHNRG/TSO20190328/TSO20190328_TSO_renamed.vcf | java -jar ~/HNRG-pipeline-V0.1/tools/SnpSift.jar filter "(GEN[$id].GT!='./.')&(GEN[$id].GT != '0/0')" > $output
done


for i in 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
    do echo $i;
    path=/home/ariel/Dropbox/UBIT/TSO/TSO-20190328/ResultsVCFs/$i
    mkdir $path
    mv /home/ariel/tmp/$i*xlsx $path/
    done



###############  FACETEO EL VCF ANOTADO PARA BPLAT ###############################################
## use snpsift for splitting multisample files. 

# rename numeric id samples
renamed_vcf=TSO20190328_TSO_renamed.vcf

for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
do
if [[ $i =~ ^1.* ]]
then
    id='ID'$i
else
    id=$i
fi
echo $id;
output=/home/hnrg/resultsHNRG/$i/$i'_TSO20190328_TSO_renamed_one_sample_SnpSift.vcf';
cat $renamed_vcf | java -jar ~/HNRG-pipeline-V0.1/tools/SnpSift.jar filter "(GEN[$id].GT!='./.')&(GEN[$id].GT != '0/0')" > aux.vcf
cat <(grep '^##' aux.vcf) <(grep -v '^##' aux.vcf| csvcut -t -c '#CHROM',POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,$id | csvformat -T) > $output
rm aux.vcf
done




for i in 1805817 EB802;
do
inputvcf=/home/hnrg/resultsHNRG/$i/$i'_TSO20190328_TSO_renamed_one_sample_SnpSift.vcf';
annotat
sed -e "s|__input__|$vcf|g" template_hnrg-anotacion_funcional_sin_CADD.json > $i'_hnrg-anotacion_funcional_sin_CADD.json'


##### NOT RUN
for i in 1711242 1802140 1804642 1805817 1809341 1809720 1810255 1901364 1901981 EB761 EB790 EB802;
    do
    if [[ $i =~ ^1.* ]]
    then
        id='ID'$i
    else
        id=$i
    fi
    echo $id;
    
    sed -i -e "s|\. |   |g" /home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv
    sed -i -e "s|\. |   |g" /home/hnrg/resultsHNRG/$i/$i'_geneList'.multiano_multisample.tsv
    done