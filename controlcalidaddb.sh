#!/bin/bash


####programa en bash para control de calidad en bases de datos
#### uit dic 2018

### patron a buscar
lista=("SNPEFF_CODON_CHANGE" "RS" "HW" \
 "GWASCAT_" "HOM" "HET" "MAF" \
 "PhastCons" "CLNSIG" "PGKB_" "AC_Adj")

 ###de que step del pipeline de anotacion
db=("paso_1" "dbSNP151_step3" "1000G_step4" \
 "GWASCAT_step6" "VARTYPE_step8_HOM" "VARTYPE_step8_HET" "EVS_step9" \
 "PhastCons_step10" "CLINVAR_step12" "PHARMGKB_step13" "EXAC_finalstep")


echo 'Base de datos' $1 $2 > reporte.txt
for ((i=0;i< ${#lista[@]}; i++))
do 
echo ${db[i]} > db1.txt
awk '{if(/'${lista[$i]}'/) {++si }else {++no}} END {print "aparece", si; print "noAparece", no; print "Nºvariantes", NR}' $1 > col1.txt
awk '{if(/'${lista[$i]}'/) {++si }else {++no}} END {print "aparece", si; print "noAparece", no; print "Nºvariantes", NR}' $2 | cut -f2 -d " " > col2.txt

paste -d'\t' col1.txt col2.txt > ${db[$i]}.txt

cat ${db[$i]}.txt | head -n1 | cut -d " " -f2 > cabecera.txt
paste -d'\t' db1.txt cabecera.txt >> reporte.txt

done

echo 'reporte generado al comparar VCFs' >> reporte.txt
 rm col1.txt col2.txt cabecera.txt db1.txt
