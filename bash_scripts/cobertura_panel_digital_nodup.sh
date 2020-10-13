#!/bin/bash
set -e

bam=$1
intervalo_restricted=$2
lista_genes=$3

#out=$4

sample_name=$(basename $bam .bam)
panel_name=$(basename $lista_genes .genes)

python /home/hnrg/NGStools/python_scripts/panel_virtual.py -l $lista_genes -ic $intervalo_restricted -o intervalo_panel_digital.bed -gne $panel_name'_no_encontrado.tsv' -ge $panel_name'_incluidos.tsv'
sed -i 's/\t*$//' intervalo_panel_digital.bed > bed_ok.bed
#### quito duplicados
/home/hnrg/HNRG-pipeline-V0.1/tools/samtools-1.9/samtools view -bh -F1024 $bam  > bam_nodups.bam
/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -b bam_nodups.bam  -a bed_ok.bed -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.genome -hist -sorted > $sample_name.hist.aux1

####con dup
#/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -b $bam  -a intervalo_panel_digital.bed -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.genome -hist -sorted > $sample_name.hist.aux1

echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
grep -v '^all' $sample_name.hist.aux1 > $sample_name.hist.aux2
cat header.txt $sample_name.hist.aux2 > $sample_name.hist
rm $sample_name.hist.aux1 $sample_name.hist.aux2 header.txt bam_nodups.bam

#rm $sample_name.hist.aux1 $sample_name.hist.aux2 header.txt 


#make tsv coverage report by exon
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/global_coverage_report_digital_panel.py -i=$sample_name.hist -o $panel_name'_digital_panel_coverage.tsv' -pn=$panel_name
       

#rm $sample_name.hist

