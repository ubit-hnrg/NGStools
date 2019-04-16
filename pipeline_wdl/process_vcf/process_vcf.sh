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

    head -n1 $out_prefix.hg19_multianno.txt > $out_prefix.hg19_multianno.tsv
    tail -n+2 $out_prefix.hg19_multianno.txt|cut -f163,164,165 --complement >>  $out_prefix.hg19_multianno.tsv;
    vcf_header=$(grep '#CH' /home/hnrg/resultsHNRG/$i/TSO20190328_TSO_renamed.$id.vcf);
    sed -i "s/Otherinfo/$vcf_header/g" $out_prefix.hg19_multianno.tsv;
    python /home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py --multianno_tsv=$out_prefix.hg19_multianno.tsv --vcf_multisample=/home/hnrg/resultsHNRG/TSO20190328/TSO20190328_TSO_renamed.vcf --output=/home/hnrg/resultsHNRG/$i/$i.multiano_multisample.tsv
    done


#este vcf se puede usar para anotar. 
#El header debe skipear tres columnas de other info 163,164,y 165. 
#Luego las restantes (166 hasta el final) son las columnas del vcf de entrada.
# la ventaja de hacerlo así es que naturalmente quedan sólo los alternativos correspondientes a la muestra que se está mirando. 

