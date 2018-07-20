#!/bin/bash

#Get names
casename=$1
sample=$2

echo $casename
echo $sample

#Paths
casepath=/home/bitgenia/samples/data/$casename
#samplepath=$casepath/$sample
samplepath=$casepath
resultpath=$samplepath/result
finalpath=$samplepath/result/$sample\_Result

echo "Remove fastq"
rm $samplepath/*.fastq

echo "Compress logs"
tar -C $resultpath -czf $resultpath/$sample\_logs.tar.gz $sample\_logs $sample\_time_logs
rm -r $resultpath/$sample\_logs $resultpath/$sample\_time_logs

mkdir $finalpath/backup

echo "Compress vcfs"
tar -czf $finalpath/$sample\_vcfs.tar.gz -C $finalpath $sample\_variant_calling.vcf $sample\_final_annot.vcf

echo "Backup important files: realigned.bam, vcfs, fastqc, metrics"
mv $finalpath/*.txt $finalpath/*.zip $finalpath/*realigned.ba* $finalpath/*.tar.gz $finalpath/backup/

echo "Remove temporary files"
rm $finalpath/*
mv $finalpath/backup/* $finalpath
rm -r $finalpath/backup

