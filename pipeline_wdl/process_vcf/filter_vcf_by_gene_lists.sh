#!/bin/bash
set -e
runID=$1
declare -a samples=("CC1807815" "CH1803670" "CH1805038" "CO1886835" "CT1802289" "CT1900424" "EB456" "EB665" "EB666" "ECM1803301" "ECM1803562" "IDP1902830")
#runID="TSO20190328"
toolpath=/home/hnrg/HNRG-pipeline-V0.1/tools
ngstools_path=/home/hnrg/NGStools
date=$(echo $runID|sed 's|TSO||')


parseshets=/home/hnrg/NGStools/pipeline_wdl/process_vcf/parse-sheets.py
make_excel=/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py
exon_coverage=/home/hnrg/resultsHNRG/$runID/$i/$i'_coverage_statistics_by_exon.tsv'

#############################
#iterate across samples. 
#################################


## task1: genero listas de genes 
python $parseshets -i=/home/hnrg/metadataHNRG/$runID/TSO$date-bioinfo.xlsx -o=/home/hnrg/metadataHNRG/$runID/


#scatter
for i in "${samples[@]}";
do 
    echo $i;
    # prepare id  ## This is only for supportting samplenames starting with numbers. (SnpSift compatibility)
    if [[ $i =~ ^[0-9].* ]]
    then
        id='ID'$i
    else
        id=$i
    fi
    echo $id;
    
    filtered=/home/hnrg/resultsHNRG/$runID/$i/$i'_geneList'.multiano_multisample.tsv
    variants=/home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv
    output_xlsx=/home/hnrg/resultsHNRG/$runID/$i/$i.variants.xlsx
    metadata_path=/home/hnrg/metadataHNRG/$runID/$i/$i


    head -n1 /home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv > $filtered
    grep /home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv -f $metadata_path'_genelist.txt' >> $filtered
    sed -i -e "s|\. |   |g" $filtered
    
    ### merge results into a xlsx file
    $make_excel $variants:Variants \
    $filtered:GeneListVariants $exon_coverage:ExonCoverage $output_xlsx


    #check if any gene do not mach with our list
    grep $metadata_path'_genelist.txt' -v -f refGene.list > $metadata_path'NotIn_refGene.txt';
done
