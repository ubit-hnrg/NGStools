#!/bin/bash
repo_path='/home/hnrg/HNRG-pipeline-V0.1'
filename=$(basename "$1")
basename=${filename%%.*}
merging_tolerance=$2   ## second argument refere to merge tolerance length
padding=100
out_basename=$basename'_padded_'$padding
chrlen='/home/hnrg/NGStools/pipeline_wdl/interval_lists/chromosome_lengths_hg19.txt'


# padding bedfile
bedtools slop  -i $1 -b $padding -g $chrlen > $out_basename'.bed'

# merge it 
sort -k1,1 -k2,2n $out_basename'.bed' | mergeBed -d $merging_tolerance > $out_basename'_merged_'$merging_tolerance'.bed'
# this line if necessary only in hg19 nomenclature
#sort -k1,1 -k2,2n $out_basename'.bed' | mergeBed -d $merging_tolerance | awk -F":" 'sub(/^chr/, "")' > $out_basename'_merged_'$merging_tolerance'.bed'

sudo java -jar $repo_path/tools/gatk-package-4.0.8.1-local.jar BedToIntervalList -I=$out_basename'_merged_'$merging_tolerance'.bed' -O=$out_basename'_merged_'$merging_tolerance'_preprocessing.interval_list' -SD=$repo_path/references/hs37d5/hs37d5.dict

