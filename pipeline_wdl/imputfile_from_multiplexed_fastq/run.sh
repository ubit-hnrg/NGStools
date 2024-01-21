#/bin/bash

#####################################################################################
# run as follow: 
#./run.sh "/home/hnrg/samplesHNRG/TSO20190425/TSO_FONARSEC-125446323/FASTQ_*/*/*.gz"

#### "DO NOT FORGET THE QUOTES!!!!"

#####################################################################################


path_to_samples=$1
tso=$(echo $path_to_samples |cut -f5 -d'/')
echo 'samples for run '$tso

path=/home/hnrg/executionsHNRG/$tso/inputs/;
mkdir -p $path;
readlink -f $path_to_samples | ./generate_fasq_tsv.sh;
mv ./samples.txt  $path$tso'_samples.txt'
