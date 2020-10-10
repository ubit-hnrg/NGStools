#!/bin/bash
set -e

bam=$1
intervalo_restricted=$2
lista_genes=$3

#out=$4

sample_name=$(basename $bam .bam)
panel_name=$(basename $lista_genes .genes)

python /home/hnrg/NGStools/python_scripts/panel_virtual.py -l $lista_genes -ic $intervalo_restricted -o intervalo_panel_digital.bed

/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -abam $bam  -b intervalo_panel_digital.bed -hist -sorted > $sample_name.$panel_digital.hist.aux1
echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
grep -v '^all' $sample_name.$panel_digital.hist.aux1 > $sample_name.$panel_digital.hist.aux2
cat header.txt $sample_name.$panel_digital.hist.aux2 > $sample_name.$panel_digital.hist
rm $sample_name.$panel_digital.hist.aux1 $sample_name.$panel_digital.hist.aux2 header.txt


# make tsv coverage report by exon
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/local_coverage_report_ENS_intersect_Library.py -i=$sample_name.$panel_digital.hist -o $sample_name_$panel_digital_coverage.tsv -s=$sample_name
       



