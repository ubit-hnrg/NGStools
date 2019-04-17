####paso1
task bam_depth {

File input_bam
File Exon_coords
#String sampleID

###herramientas
#String gatk_jar
String toolpath

command <<<

#!/bin/bash
set -e
set -o pipefail
name=$(basename ${input_bam} .bam)

###input
#BAM='/home/hnrg/resultsHNRG/'$sampleID'/'$sampleID'.bam'  ###comentarrrrrrrrrrrrr


  
###input
#EXON_coords=/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_in_TSO_V1_padded.bed ###comentarrrrrrrrrrrrr


#toolpath=/home/hnrg/HNRG-pipeline-V0.1/tools
# ok for reduce bam
${toolpath}bedtools2/bin/intersectBed -a ${input_bam} -b ${Exon_coords} > $name'_exonTSO_reduced.bam'

#Bam coverage. No está estrictamente limitado al intervalo porque contiene grandes zonas con cobertura CERO 
#que están presentes por una mínima interseccion con el intervalo buscado)
${toolpath}bedtools2/bin/genomeCoverageBed -ibam $name'_exonTSO_reduced.bam' -bga > $name'_exon.coverage'

# Left outer join EXON_coords and BAM coverage. Intermediate file. Incorpora a cada intervalo la covertura, (pero sigue preservando zonas con mínima intersección que traen ruido)
${toolpath}bedtools2/bin/intersectBed -a ${Exon_coords} -b $name'_exon.coverage' -loj > $name'_loj.txt'


#esto limita la covertura estrictamente al intervalo paddeado.
# Es decicr, esta linea elimina las grandes zonas del bam con cobertura zero que tenían unas pocas bases de interseccion con el intervalo buscado
# primero: reordeno las columnas de loj.txt, para que las coordenadas start end no sean las del exon, sino las de la cobertura
awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$5"\t"$1"\t"$2"\t"$3}' $name'_loj.txt' > $name'_loj_sorted_cols.tsv'
${toolpath}bedtools2/bin/intersectBed -a $name'_loj_sorted_cols.tsv' -b ${Exon_coords}  > $name'_loj_exon_filtered.coverage'

echo -e 'chr\tstart\tend\tcount\tgeneSymbol\texon_number\texon_chr\texon_start\texon_end' > header.tsv
cat header.tsv $name'_loj_exon_filtered.coverage' > $name'_exon_filtered_coverage.tsv'
rm $name'_loj_exon_filtered.coverage' $name'_loj_sorted_cols.tsv' header.tsv $name'_exon.coverage' $name'_loj.txt'


${toolpath}/qualityControl/coverage_statistics.py $name'_exon_filtered_coverage.tsv'

>>>


output {

#File reduced_bam = $name'_exonTSO_reduced.bam'
#File exon_coverage = $name'_exon_filtered_coverage.tsv'
File glob_cov_stats = "$name'_coverage_statistics.tsv'"
File cov_stats_by_name = "$name'_coverage_statistics_by_exon.tsv'"
String sample_Name = "name"

}

}


#####task del summary_metrics de samtools
task samtools_stat{

###herramientas
#String gatk_jar
String toolpath
File TSO_bed #./TruSight_One_v1_padded_100_GRCh37.bed
File input_orig_bam
#

command {

##paso1
name=$(basename ${input_orig_bam} .bam)

${toolpath}samtools stats ${input_orig_bam}  > $name'_orig_samtools.stats'

#paso2
##input es el bam original  y el intervalo de captura
${toolpath}samtools stats ${input_orig_bam} -t ${TSO_bed} > $name'_TSO_samtools.stats'

}
output {

#File samtools_stats = 
#File samtools_reduced_bam = $name'_samtools_reduced.stats'
File samtools_original_bam = "$name'_orig_samtools.stats'"
File samtools_TSO_bam = "$name'_TSO_samtools.stats'"

#File 

}

}

#README_run_merge_coverage_global_reports.sh
#task merge coverage_global_reports
# /home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py coverage_global_statistics.files TSO20190328_coverage_statistics.tsv

task coverage_global_reports {
####inputs del paso1 
File coverage_global_stats
File coverage_stats 


command{
python << CODE 
#!/usr/bin/python
import sys
import pandas as pd
import os

reports = ${coverage_global_stats}
reporte_salida = ${coverage_stats}
merged_report = []

reports = pd.read_csv(reports).values.flatten()
for rep in reports:
    print rep
    sample = os.path.basename(rep).split('_')[0]
    print sample
    df = pd.read_table(rep,index_col=[0])
    df.columns=[sample]
    merged_report.append(df)

merged_report = pd.concat(merged_report,axis = 1)
merged_report.index.name=None
merged_report.to_csv(reporte_salida,sep ='\t')
CODE

}

output {
File merged_report = "reporte_salida"

}

}



##/home/hnrg/NGStools/pipeline_wdl/qualityControl/samtools_stats_report.py

task reports_file {

String sampleID
#File sample
File samtools_global_report
File samtools_report
String toolpath

#String path_salida 

command {
${toolpath}/qualityControl/samtools_stats_report.py -g=${samtools_global_report} -l=${samtools_report} -o=${sampleID}_samtools_report.tsv

}

output {
 
File output_report = "output_rep" 

}

}

workflow quality_control {

File analysis_readybam 
String toolpath
File exon_coords
File tso_bed

call bam_depth {
input: 

input_bam = analysis_readybam,
Exon_coords = exon_coords,
toolpath = toolpath

}


call samtools_stat {
input:
toolpath = toolpath,
TSO_bed = tso_bed, #./TruSight_One_v1_padded_100_GRCh37.bed
input_orig_bam = analysis_readybam


}


call coverage_global_reports {

input: 
####inputs del paso1 
coverage_global_stats = bam_depth.glob_cov_stats,
coverage_stats = bam_depth.cov_stats_by_name

} 

call reports_file {

input: 
sampleID = bam_depth.sample_Name,
samtools_global_report = samtools_stat.samtools_original_bam,
samtools_report = samtools_stat.samtools_TSO_bam,
toolpath = toolpath

}

}