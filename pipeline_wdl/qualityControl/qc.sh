#!/bin/bash
set -e

input_bam=/data/resultsHNRG/*/CC1707556/CC1707556.bam
TSO_bed=/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/TruSight_One_v1.1_GRCh37.bed
ensembl=/home/hnrg/HNRG-pipeline-V0.1/libraries/intervalos/ensembl_canonicos_GRCh37_0based.tsv 
name=$(basename $input_bam)
toolpath=/home/hnrg/HNRG-pipeline-V0.1/tools

# esto reporta la cobertura en cada intervalo de captura y hace un histograma global también con el keyword "all"
$toolpath/bedtools2/bin/coverageBed -a $TSO_bed -b $input_bam  -hist > $name.hist.aux
echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
cat header.txt $name.hist.aux > $name.hist; 
rm $name.hist.aux header.txt

#histograma global del bam restringido a toda la librería
grep '^all' $name.hist > global.hist
echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global.header.txt
cat global.header.txt global.hist > $name.global.hist
rm global.header.txt global.hist

# make tsv report
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/global_coverage_report_inLibrary.py -i=$name.global.hist -o $name'_experiment_global_report.tsv' -s $name




#### EXONES ####
#histograma restringido a cada exon de ensembl que está en la librería de captura
$toolpath/bedtools2/bin/intersectBed -a $ensembl -b $TSO_bed > ensembl_in_library.tsv
exon_coords=./ensembl_in_library.tsv

$toolpath/bedtools2/bin/coverageBed -a $exon_coords -b $input_bam  -hist > $name.ENS.hist.aux1
echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
grep -v '^all' $name.ENS.hist.aux1 > $name.ENS.hist.aux2
cat header.txt $name.ENS.hist.aux2 > $name.ENS.hist;
rm $name.ENS.hist.aux1 $name.ENS.hist.aux2 rm header.txt

# make tsv coverage report by exon
python /home/hnrg/NGStools/pipeline_wdl/qualityControl/local_coverage_report_ENS_intersect_Library.py -i=$name.ENS.hist -o $name'_ENS_local_report.tsv' -s=$name



####################################################
################# NO VA MAS ########################

## ok for reduce bam
#Exon_coords=/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_GRCh37.tsv

##$toolpath/bedtools2/bin/intersectBed -a $input_bam -b $Exon_coords > $name'_exonTSO_reduced.bam'

#Bam coverage. No está estrictamente limitado al intervalo porque contiene grandes zonas con cobertura CERO 
#que están presentes por una mínima interseccion con el intervalo buscado)
##$toolpath/bedtools2/bin/genomeCoverageBed -ibam $name'_exonTSO_reduced.bam' -bga > $name'_exon.coverage'

# Left outer join EXON_coords and BAM coverage. Intermediate file. Incorpora a cada intervalo la covertura, (pero sigue preservando zonas con mínima intersección que traen ruido)
##$toolpath/bedtools2/bin/intersectBed -a $Exon_coords -b $name'_exon.coverage' -loj > $name'_loj.txt'


#esto limita la covertura estrictamente al intervalo paddeado.
# Es decicr, esta linea elimina las grandes zonas del bam con cobertura zero que tenían unas pocas bases de interseccion con el intervalo buscado
# primero: reordeno las columnas de loj.txt, para que las coordenadas start end no sean las del exon, sino las de la cobertura
##awk -F"\t" '{print $8"\t"$9"\t"$10"\t"$11"\t"$5"\t"$4"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3}' ${name}_loj.txt > ${name}_loj_sorted_cols.tsv

#$toolpath/bedtools2/bin/intersectBed -a $name'_loj_sorted_cols.tsv' -b $Exon_coords  > $name'_loj_exon_filtered.coverage'
#para ensembl
#echo -e 'chr\tstart\tend\tcount\tgeneSymbol\tENSEMBL_ID\texon_number\tstrand\texon_chr\texon_start\texon_end' > header.tsv
#cat header.tsv $name'_loj_exon_filtered.coverage' > $name'_exon_filtered_coverage.tsv'
#rm ${name}_loj_exon_filtered.coverage ${name}_loj_sorted_cols.tsv header.tsv ${name}_exon.coverage ${name}_loj.txt

#####coverage statistics cambio la forma del input... 
###usamos este con ENS
#/home/hnrg/NGStools/pipeline_wdl/qualityControl/coverage_statistics_v1.0_ENS.py -i $name'_exon_filtered_coverage.tsv' -g $name'_global_coverage_statistics.tsv' -e $name'_coverage_statistics_by_exon.tsv' -s $name

#rm ${name}_exonTSO_reduced.bam ${name}_exon_filtered_coverage.tsv

####################################################
################# NO VA MAS ########################


# ACÁ RETOMARÍAS con samtools

#####task del summary_metrics de samtools
#input_orig_bam = 
#nameorig = basename(input_orig_bam, ".bam")
