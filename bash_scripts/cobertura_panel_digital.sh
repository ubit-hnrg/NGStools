#!/bin/bash
set -e

bam=$1
intervalo_restricted=$2
lista_genes=$3

out=$4

sample_name=$(basename $bam .bam)
panel_name=$(basename $3 .genes)

python /home/hnrg/NGStools/python_scripts/panel_virtual.py -l $3 -ic $2 -o intervalo_panel_digital.bed

/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -a intervalo_panel_digital -b $1 -hist -sorted > $sample_name.$panel_digital.hist.aux1
echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
grep -v '^all' $sample_name.$panel_digital.hist.aux1 > $sample_name.$panel_digital.hist.aux2
cat header.txt $sample_name.$panel_digital.hist.aux2 > $sample_name.$panel_digital.hist.hist
rm $sample_name.$panel_digital.hist.aux1 $sample_name.$panel_digital.hist.aux2 header.txt

