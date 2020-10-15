#!/bin/bash
set -e

bam=$1
intervalo_restricted=$2
lista_genes=$3

#out=$4

sample_name=$(basename $bam .bam)
panel_name=$(basename $lista_genes .genes)

####for intervalo de captura raw.
#sort -k1,1V -k2,2n $4 > intervalo_sorted.bed

python /home/hnrg/NGStools/python_scripts/panel_virtual.py -l $lista_genes -ic $intervalo_restricted -o intervalo_panel_digital.bed -gne $panel_name'_no_encontrado.tsv' -ge $panel_name'_incluidos.tsv'

#### quito duplicados
/home/hnrg/HNRG-pipeline-V0.1/tools/samtools-1.9/samtools view -h -F1024 $bam -u > bam_nodups.sam
/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -b bam_nodups.sam  -a intervalo_panel_digital.bed -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.genome -hist -sorted > $sample_name'.hist.aux'

#####fines de testing: cobertura en la libreria del experimento: NO DUPS
#/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -b bam_nodups.sam  -a intervalo_sorted.bed -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.genome -hist -sorted > $sample_name'.intervalo_captura.aux'
#samtools stats bam_nodups.sam -t intervalo_sorted.bed > $sample_name'_samtools_nodups.stats'



#####fines de testing: cobertura en la libreria del experimento: DUPS
#/home/hnrg/HNRG-pipeline-V0.1/tools/bedtools2/bin/coverageBed -b $bam  -a intervalo_sorted.bed -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.genome -hist -sorted > $sample_name'.intervalo_captura_dup.aux'
#samtools stats $bam -t intervalo_sorted.bed > $sample_name'_samtools_dups.stats'

#####para ensembl
echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
grep -v '^all' $sample_name'.hist.aux' > $sample_name'.hist.aux2'
cat header.txt $sample_name'.hist.aux2' > $sample_name'.hist'
rm $sample_name'.hist.aux2' header.txt bam_nodups.sam 

#histograma global del bam nodup restringido a toda la librería ####all, no lo uso. este va en el excel
#grep '^all' $sample_name'.hist.aux' > global_nodup.hist
#echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global_nodup.header.txt
#cat global_nodup.header.txt global_nodup.hist > $sample_name'_global_nodup.hist'
#rm global_nodup.header.txt global_nodup.hist $sample_name'.hist.aux'

# ########################################################### para intervalo_lib.
# echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header_lib.txt
# grep -v '^all' $sample_name'.intervalo_captura.aux' > $sample_name'.IC.aux2'
# cat header_lib.txt $sample_name'.IC.aux2' > $sample_name'_lib_IC.hist' 
# rm $sample_name'.IC.aux2' header_lib.txt

# #histograma global del bam nodup restringido a toda la librería
# grep '^all' $sample_name'.intervalo_captura.aux' > global_all.hist
# echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global_nodup.header.txt
# cat global_nodup.header.txt global_all.hist > $sample_name'_IC_all.hist'
# rm global_nodup.header.txt global_all.hist

# #####intervalo sin duplicados $sample_name'.intervalo_captura_dup.aux'
# echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header_lib_dup.txt
# grep -v '^all' $sample_name'.intervalo_captura_dup.aux' > $sample_name'.IC_dup.aux2'
# cat header_lib_dup.txt $sample_name'.IC_dup.aux2' > $sample_name'_lib_IC_dup.hist' 
# rm $sample_name'.IC_dup.aux2' header_lib_dup.txt

# #histograma global del bam nodup restringido a toda la librería
# grep '^all' $sample_name'.intervalo_captura_dup.aux' > global_all_dup.hist
# echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global_dup.header.txt
# cat global_dup.header.txt global_all_dup.hist > $sample_name'_IC_all_dup.hist'
# rm global_dup.header.txt global_all_dup.hist


#make tsv coverage report by exon
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/global_coverage_report_digital_panel.py -i=$sample_name'.hist' -o $panel_name'_digital_panel_coverage.tsv' -pn=$panel_name
       
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py $panel_name'_digital_panel_coverage.tsv':$panel_name $panel_name'_no_encontrado.tsv':genes_NO_encontrados $panel_name'_incluidos.tsv':genes_cubiertos  $panel_name'_digital_coverage_report.xlsx'


rm $sample_name.hist intervalo_panel_digital.bed $panel_name'_digital_panel_coverage.tsv'

