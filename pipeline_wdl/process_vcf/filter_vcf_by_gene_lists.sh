#!/bin/bash
set -e
runID=$1
declare -a samples=("CH1803670" "CH1805038" "CO1886835" "CT1802289" "CT1900424" "EB456" "EB665" "EB666" "ECM1803301" "ECM1803562" "IDP1902830")

#runID="TSO20190328"
toolpath=/home/hnrg/HNRG-pipeline-V0.1/tools
ngstools_path=/home/hnrg/NGStools
date=$(echo $runID|sed 's|TSO||')


parseshets=/home/hnrg/NGStools/pipeline_wdl/process_vcf/parse-sheets.py
make_excel=/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py
refGenelist=/home/hnrg/metadataHNRG/refGene.list
#exon_coverage=/home/hnrg/resultsHNRG/$runID/$i/$i'_coverage_statistics_by_exon.tsv'
#############################
#iterate across samples. 
#################################


## task1: genero listas de genes 
python $parseshets -i=/home/hnrg/metadataHNRG/$runID/TSO-$date'-bioinfo.xlsx' -o=/home/hnrg/metadataHNRG/$runID/


#scatter
for i in "${samples[@]}";
    do echo $i; 
    # prepare id  ## This is only for supportting samplenames starting with numbers. (SnpSift compatibility)
    if [[ $i =~ ^[0-9].* ]]
    then
        id='ID'$i
    else
        id=$i
    fi
    
    filtered=/home/hnrg/resultsHNRG/$runID/$i/$i'_geneList'.multiano_multisample.tsv
    variants=/home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv
    output_xlsx=/home/hnrg/resultsHNRG/$runID/$i/$i.variants.xlsx
    out_bad_genes=/home/hnrg/resultsHNRG/$runID/$i/$i'_NotIn_refGene.txt'
    metadata_path=/home/hnrg/metadataHNRG/$runID/$i/$i
    exon_coverage=$(readlink -f /home/hnrg/executionsHNRG/$runID/cromwell-executions/main_workflow/7b733377-eb59-4ea3-81e7-13f9594d2d5b/call-quality_control/qual_control.quality_control/ffc1d0a9-f582-4a47-a47e-9e46961817cf/call-bam_depth/*/execution/$i'_coverage_statistics_by_exon.tsv')
    

    head -n1 /home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv > $filtered
    grep /home/hnrg/resultsHNRG/$runID/$i/$i.multiano_multisample.tsv -f $metadata_path'_genelist.txt' >> $filtered
    sed -i -e "s|\. |   |g" $filtered
    
    ### merge results into a xlsx file
    $make_excel $variants:Variants \
    $filtered:GeneListVariants $exon_coverage:ExonCoverage $output_xlsx


    #check if any gene do not mach with our list
    grep $metadata_path'_genelist.txt' -v -f $refGenelist > $out_bad_genes;


done
