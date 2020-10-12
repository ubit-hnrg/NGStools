#!/bin/bash
set -e

bam=$1
intervalo_restricted=$2
lista_genes=$3

#out=$4

sample_name=$(basename $bam .bam)
panel_name=$(basename $lista_genes .genes)

python /home/hnrg/NGStools/python_scripts/panel_virtual.py -l $lista_genes -ic $intervalo_restricted -o intervalo_panel_digital.bed

#### quito duplicados
##/home/hnrg/HNRG-pipeline-V0.1/tools/samtools-1.9
samtools view -F1024 $bam> bam_nodups.bam
samtools stats bam_nodups.bam -t intervalo_panel_digital.bed  > $sample_name'_samtools_nodup.stats'


/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -b bam_nodups.bam  -a intervalo_panel_digital.bed -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.genome -hist -sorted > $sample_name.hist.aux1
echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
grep -v '^all' $sample_name.hist.aux1 > $sample_name.hist.aux2
cat header.txt $sample_name.hist.aux2 > $sample_name.hist
rm $sample_name.hist.aux1 $sample_name.hist.aux2 header.txt bam_nodups.bam


# make tsv coverage report by exon
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/local_coverage_report_pd_20X.py -i=$sample_name.hist -o $sample_name'_digital_panel_coverage.tsv' -pn=$panel_name
       



