#/bin/bash
tso=TSO20190425;
path=/home/hnrg/executionsHNRG/$tso/inputs/;
mkdir -p $path;
readlink -f /home/hnrg/samplesHNRG/$tso/TSO_FONARSEC-125446323/FASTQ_*/*/*.gz|sort | ./generate_fasq_tsv.sh;
mv ./samples.txt  $path$tso'_samples.txt'
