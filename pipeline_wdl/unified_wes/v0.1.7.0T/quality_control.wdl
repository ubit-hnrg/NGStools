###fastp

task fastp_qual {
Array[File] inputs_json_report
String Tso_name

command <<<
./estadistica_fastp.py -i ${inputs_json_report} -o ${Tso_name}_fastp_report.tsv
>>>


output {
File fastp_stats = "${Tso_name}_fastp_report.tsv"

}
}


####paso1 para calidad de bams
task bam_depth {

File input_bam
File Exon_coords
#String sampleID

###herramientas
String name = basename(input_bam ,".bam")
String toolpath

command <<<

#!/bin/bash
set -e
set -o pipefail
#name=$(basename ${input_bam} .bam)
###input
#BAM='/home/hnrg/resultsHNRG/'$sampleID'/'$sampleID'.bam'  ###comentarrrrrrrrrrrrr

###input
#EXON_coords=/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/exon_coords_in_TSO_V1_padded.bed ###comentarrrrrrrrrrrrr

#toolpath=/home/hnrg/HNRG-pipeline-V0.1/tools
# ok for reduce bam
${toolpath}bedtools2/bin/intersectBed -a ${input_bam} -b ${Exon_coords} > ${name}'_exonTSO_reduced.bam'

#Bam coverage. No está estrictamente limitado al intervalo porque contiene grandes zonas con cobertura CERO 
#que están presentes por una mínima interseccion con el intervalo buscado)
${toolpath}bedtools2/bin/genomeCoverageBed -ibam ${name}'_exonTSO_reduced.bam' -bga > ${name}'_exon.coverage'

# Left outer join EXON_coords and BAM coverage. Intermediate file. Incorpora a cada intervalo la covertura, (pero sigue preservando zonas con mínima intersección que traen ruido)
${toolpath}bedtools2/bin/intersectBed -a ${Exon_coords} -b ${name}'_exon.coverage' -loj > ${name}'_loj.txt'


#esto limita la covertura estrictamente al intervalo paddeado.
# Es decicr, esta linea elimina las grandes zonas del bam con cobertura zero que tenían unas pocas bases de interseccion con el intervalo buscado
# primero: reordeno las columnas de loj.txt, para que las coordenadas start end no sean las del exon, sino las de la cobertura
awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$5"\t"$1"\t"$2"\t"$3}' ${name}'_loj.txt' > ${name}'_loj_sorted_cols.tsv'
${toolpath}bedtools2/bin/intersectBed -a ${name}'_loj_sorted_cols.tsv' -b ${Exon_coords}  > ${name}'_loj_exon_filtered.coverage'

echo -e 'chr\tstart\tend\tcount\tgeneSymbol\texon_number\texon_chr\texon_start\texon_end' > header.tsv
cat header.tsv ${name}'_loj_exon_filtered.coverage' > ${name}'_exon_filtered_coverage.tsv'
rm ${name}'_loj_exon_filtered.coverage' ${name}'_loj_sorted_cols.tsv' header.tsv ${name}'_exon.coverage' ${name}'_loj.txt'

#####coverage statistics cambio la forma del input... 
/home/hnrg/NGStools/pipeline_wdl/qualityControl/coverage_statistics_v1.0.py -i ${name}'_exon_filtered_coverage.tsv' -g ${name}_global_coverage_statistics.tsv -e ${name}_coverage_statistics_by_exon.tsv -s ${name}


rm ${name}'_exonTSO_reduced.bam' ${name}'_exon_filtered_coverage.tsv'
>>>


output {

#File reduced_bam = $name'_exonTSO_reduced.bam'
#File exon_coverage = $name'_exon_filtered_coverage.tsv'


#### agregar un writelines para generar un archivo .files con el 

File cov_stats_by_exon = "${name}_coverage_statistics_by_exon.tsv"
File glob_cov_stats = "${name}_global_coverage_statistics.tsv"
String sample_Name = "${name}"

}

}

#task stat_files {
#Array[File] files_in
#String sample_name 


#command <<<
#mv ${write_lines(files_in)}  ${sample_name}.txt

#>>>

#output {
#File path_stat_files = "${sample_name}.txt"

#}

#}


#####task del summary_metrics de samtools
task samtools_stat{

###herramientas
#String gatk_jar

String toolpath
File TSO_bed #./TruSight_One_v1_padded_100_GRCh37.bed
File input_orig_bam
String name = basename(input_orig_bam ,".bam")

command {

##paso1
###$name=$(basename ${input_orig_bam} .bam)

${toolpath}samtools stats ${input_orig_bam}  > ${name}'_orig_samtools.stats'

#paso2
##input es el bam original  y el intervalo de captura
${toolpath}samtools stats ${input_orig_bam} -t ${TSO_bed} > ${name}'_TSO_samtools.stats'

}
output {

#File samtools_stats = 
#File samtools_reduced_bam = $name'_samtools_reduced.stats'
File samtools_stat_original_bam = "${name}_orig_samtools.stats"
File samtools_stat_TSO_bam = "${name}_TSO_samtools.stats"

#File 

}

}





##/home/hnrg/NGStools/pipeline_wdl/qualityControl/samtools_stats_report.py

task samtools_reports_file {

String sampleID
#File sample
File samtools_global_report
File samtools_library_report
String toolpath

#String path_salida 

command {
/home/hnrg/NGStools/pipeline_wdl/qualityControl/samtools_stats_report.py -g=${samtools_global_report} -l=${samtools_library_report} -o=${sampleID}_samtools_report.tsv

}

output {
 
File output_global_report = "${sampleID}_samtools_report.tsv" 

}

}


#README_run_merge_coverage_global_reports.sh
#task merge coverage_global_reports
# /home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py coveerage_global_statistics.files TSO20190328_coverage_statistics.tsv


#########este script es alimentado por un archivo que tiene  nombre_reporte.tsv por cada linea. y se llama en uno de los ultimos pasos 

#####
task merge_coverage_global_reports {
####inputs del paso1 
Array[File] coverage_global_files
String TSO_name
#String toolpath
#File coverage_stats 


command{
/home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py -i ${coverage_global_files} -o ${TSO_name}.merged_global_report
}

output {
File merged_glob_report = "${TSO_name}.merged_global_report"

}

}

task merge_samtools_reports {
####inputs del paso1 
Array[File] samtools_reports_files
String TSO_name
#String toolpath
#File coverage_stats 


command{
/home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py -i ${samtools_reports_files} -o ${TSO_name}.merged_st_report
}

output {
File merged_st_report = "${TSO_name}.merged_st_report"

}

}

workflow quality_control {

Array[File] analysis_readybam 
String toolpath
File exon_coords
File tso_bed
Array[File] fastp_json_files
String Tso_name 

call fastp_qual {
input:
inputs_json_report = fastp_json_files,
Tso_name = Tso_name

}


######################scatter por los bams... analysis_readybam

scatter (bams_ready in analysis_readybam)  {

call bam_depth {
input: 
input_bam = bams_ready,
Exon_coords = exon_coords,
toolpath = toolpath

}

#call stat_files {
#    input:
# files_in = ["${bam_depth.glob_cov_stats}","${bam_depth.cov_stats_by_exon}"], 
# sample_name = bam_depth.sample_Name
#
#}

#File glob_cov_stats = "${name}_coverage_statistics.tsv"
#File Files_global_stats = "${name}.files"
#File cov_stats_by_name = "${name}_coverage_statistics_by_exon.tsv"
#String sample_Name = "${name}"


call samtools_stat {
input:
toolpath = toolpath,
TSO_bed = tso_bed, #./TruSight_One_v1_padded_100_GRCh37.bed
input_orig_bam = analysis_readybam

}


call samtools_reports_file {

input: 
sampleID = bam_depth.sample_Name,
samtools_global_report = samtools_stat.samtools_stat_original_bam,
samtools_library_report = samtools_stat.samtools_stat_TSO_bam,
toolpath = toolpath

}

}

Array[File] bams_stat_depth_global_coverage_stats = ["${bam_depth.glob_cov_stats}"]


####### esto mergea archivos de distintas muestras
call merge_coverage_global_reports {

input:  
#toolpath = toolpath,
coverage_global_files = bams_stat_depth_global_coverage_stats,
#sample_Name = bam_depth.sample_Name
#coverage_stats = bam_depth.cov_stats_by_name

} 

#Array[File] stat_alineamiento = ["${samtools_reports_file.output_global_report}"]







output {
Array[File] depth_global_cov_stats = bam_depth.glob_cov_stats ###estadistica del alineamiento...
Array[File] by_exon_depth = bam_depth.cov_stats_by_exon
#File coverage_merged_report = merge_coverage_global_reports.merged_report
Array[File] reporte_final = samtools_reports_file.output_global_report ### archivo para mergear... estadistica en la libreria del experimento

Array[File] Samt_bam_stat = samtools_stat.samtools_stat_original_bam 
Array[File] Samt_TSO_stat = samtools_stat.samtools_stat_TSO_bam


}

}