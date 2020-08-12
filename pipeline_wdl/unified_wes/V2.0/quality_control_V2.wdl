###fastp

task fastp_qual {
  File inputs_json_report
  String report_name = basename(inputs_json_report, ".txt")

  #${sep=' -I ' input_bqsr_reports}
  command <<<
  /home/hnrg/NGStools/pipeline_wdl/qualityControl/estadistica_fastp.py -i ${inputs_json_report} -o ${report_name}_fastp_report.tsv
  >>>

  output {
    File fastp_stats = "${report_name}_fastp_report.tsv"

  }
}


####paso1 para calidad de bams
task bam_depth {

  File input_bam
  File Exon_coords
  

  ###herramientas
  String name = basename(input_bam, ".bam")
  String toolpath
  String pipeline_version

  command <<<

  #!/bin/bash
  set -e
  set -o pipefail
  
  #reduce bam
  ${toolpath}bedtools2/bin/intersectBed -a ${input_bam} -b ${Exon_coords} > ${name}_exonTSO_reduced.bam

  ###ya que estamos, prediccion de sexo:
  /home/hnrg/NGStools/python_scripts/bam_sex_xy.py -b $input_bam > ${name}_sex.txt

  #Bam coverage. No está estrictamente limitado al intervalo porque contiene grandes zonas con cobertura CERO 
  #que están presentes por una mínima interseccion con el intervalo buscado)
  ${toolpath}bedtools2/bin/genomeCoverageBed -ibam ${name}_exonTSO_reduced.bam -bga > ${name}_exon.coverage

  # Left outer join EXON_coords and BAM coverage. Intermediate file. Incorpora a cada intervalo la covertura, (pero sigue preservando zonas con mínima intersección que traen ruido)
  ${toolpath}bedtools2/bin/intersectBed -a ${Exon_coords} -b ${name}_exon.coverage -loj > ${name}_loj.txt


  #esto limita la covertura estrictamente al intervalo paddeado.
  # Es decicr, esta linea elimina las grandes zonas del bam con cobertura zero que tenían unas pocas bases de interseccion con el intervalo buscado
  # primero: reordeno las columnas de loj.txt, para que las coordenadas start end no sean las del exon, sino las de la cobertura
  #awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$5"\t"$1"\t"$2"\t"$3}' ${name}_loj.txt > ${name}_loj_sorted_cols.tsv

  ### para el intervalo de ENSEMBL usar el de abajo y no el de arriba:
  awk -F"\t" '{print $8"\t"$9"\t"$10"\t"$11"\t"$5"\t"$4"\t"$6"\t"$7"\t"$1"\t"$2"\t"$3}' ${name}_loj.txt > ${name}_loj_sorted_cols.tsv

  ${toolpath}bedtools2/bin/intersectBed -a ${name}_loj_sorted_cols.tsv -b ${Exon_coords}  > ${name}_loj_exon_filtered.coverage

  echo -e 'chr\tstart\tend\tcount\tgeneSymbol\tENSEMBL_ID\texon_number\tstrand\texon_chr\texon_start\texon_end' > header.tsv
  cat header.tsv ${name}_loj_exon_filtered.coverage > ${name}_exon_filtered_coverage.tsv
  rm ${name}_loj_exon_filtered.coverage ${name}_loj_sorted_cols.tsv header.tsv ${name}_exon.coverage ${name}_loj.txt

  #####coverage statistics cambio la forma del input... 
  #/home/hnrg/NGStools/pipeline_wdl/qualityControl/coverage_statistics_v1.0.py -i ${name}_exon_filtered_coverage.tsv -g ${name}_global_coverage_statistics.tsv -e ${name}_coverage_statistics_by_exon.tsv -s ${name}

  ###usamos este con ENS
  /home/hnrg/NGStools/pipeline_wdl/qualityControl/coverage_statistics_v1.0_ENS.py -i ${name}_exon_filtered_coverage.tsv -g ${name}_global_coverage_statistics_${pipeline_version}.tsv -e ${name}_coverage_statistics_by_exon_${pipeline_version}.tsv -s ${name}

  rm ${name}_exonTSO_reduced.bam ${name}_exon_filtered_coverage.tsv
  >>>


  output {
    File cov_stats_by_exon = "${name}_coverage_statistics_by_exon_${pipeline_version}.tsv"
    File glob_cov_stats = "${name}_global_coverage_statistics_${pipeline_version}.tsv"
    File sex_prediction = "${name}_sex.txt"
    String sample_Name = "${name}"
  }

}



#####task del summary_metrics de samtools
task samtools_stat{

  String toolpath
  File TSO_bed #./TruSight_One_v1_padded_100_GRCh37.bed
  File input_orig_bam
  String name = basename(input_orig_bam, ".bam")

  command {

  ##paso1
  ###$name=$(basename ${input_orig_bam} .bam)

  ${toolpath}samtools stats ${input_orig_bam}  > ${name}_orig_samtools.stats

  #paso2
  ##input es el bam original  y el intervalo de captura
  ${toolpath}samtools stats ${input_orig_bam} -t ${TSO_bed} > ${name}_TSO_samtools.stats

  }
  output {

    #File samtools_stats = 
    #File samtools_reduced_bam = $name'_samtools_reduced.stats'
    File samtools_stat_original_bam = "${name}_orig_samtools.stats"
    File samtools_stat_TSO_bam = "${name}_TSO_samtools.stats"

  }

}





##/home/hnrg/NGStools/pipeline_wdl/qualityControl/samtools_stats_report.py

task samtools_reports_file {

  String sampleID
  File samtools_library_report
  String toolpath
  String pipeline_version


  command {
    /home/hnrg/NGStools/pipeline_wdl/qualityControl/samtools_stats_report_v1.0.py -l=${samtools_library_report} -o=${sampleID}_samtools_report_${pipeline_version}.tsv
  }

  output {
 
    File output_global_report = "${sampleID}_samtools_report_${pipeline_version}.tsv" 

  }

}


#README_run_merge_reports.sh
#task merge coverage_global_reports
# /home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py coveerage_global_statistics.files TSO20190328_coverage_statistics.tsv


#########este script es alimentado por un archivo que tiene  nombre_reporte.tsv por cada linea. y se llama en uno de los ultimos pasos 

#####
task merge_reports {
####inputs del paso1 
File files_to_merge
String TSO_name
#String toolpath
#File coverage_stats 

#gvcfs = ['${sep="','" input_gvcfs}']

command<<<
/home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py -i ${files_to_merge} -o ${TSO_name}.merged_report
>>>

output {
File merged_report = "${TSO_name}.merged_report"

}

}

#task merge_samtools_reports {
####inputs del paso1 
#Array[File] samtools_reports_files
#String TSO_name
#String toolpath
#File coverage_stats 


#command{
#/home/hnrg/NGStools/pipeline_wdl/qualityControl/merge_sample_reports.py -i ${samtools_reports_files} -o ${TSO_name}.merged_st_report
#}

#output {
#File merged_st_report = "${TSO_name}.merged_st_report"

#}

#}

task CreateFoFN {
  # Command parameters
  Array[File] array_of_files
  String fofn_name
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list \
   
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }
}

task make_excel {
  String pipeline_version
  String Tso_name
  File tabla1
  String pestana1
  File tabla2
  String pestana2
  File tabla3
  String pestana3

  command{
    /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${tabla1}:${pestana1} ${tabla2}:${pestana2} ${tabla3}:${pestana3} ${Tso_name}_qual_report_${pipeline_version}.xlsx
 
  }

  output {
    File reporte_excel = "${Tso_name}_qual_report_${pipeline_version}.xlsx"

  }

}


task symlink_important_files {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}



workflow quality_control_V2 {

  Array[File]+ analysis_readybam 
  String toolpath
  File exon_coords
  #File tso_bed
  Array[File] fastp_json_files
  String Tso_name
  Array[String] path_save
  String pipeline_v
  #Map[String,String] bams_N_reads

  scatter (fastp in fastp_json_files){
    call fastp_qual {
      input:
      inputs_json_report = fastp
    }
  }
 

  ######################scatter por los bams... analysis_readybam

  scatter (bams_ready in analysis_readybam) {

    call bam_depth {
      input: 
      input_bam = bams_ready,
      Exon_coords = exon_coords,
      toolpath = toolpath,
      pipeline_version = pipeline_v
    }

  }


  #Array[File] bams_sex_prediction = bam_depth.sex_prediction
  #Array[File] bams_stat_depth_global_coverage_stats = bam_depth.glob_cov_stats

  Array[File] stat_alineamiento #samtools_reporte_final from bam2gvcf


 #Create a file with a list of the generated bam_depth.glob_cov_stats
  # call CreateFoFN {
  #   input:
  #     array_of_files = bams_stat_depth_global_coverage_stats,
  #     fofn_name = Tso_name
     
  # }

 #Create a file with a list of the generated output_global_report
#   call CreateFoFN as CreateFoFN_samtools{
#     input:
#       array_of_files = stat_alineamiento,
#       fofn_name = Tso_name
     
#   }

#  call CreateFoFN as CreateFoFN_fastp{
#     input:
#       array_of_files = fastp_stats,
#       fofn_name = Tso_name
     
#   }

  ####### esto mergea archivos de distintas muestras
  # call merge_reports {

  #   input:  
  #   files_to_merge = CreateFoFN.fofn_list,
  #   TSO_name = Tso_name

  # } 

  # call merge_reports as merge_samtools_reports{

  #   input:  
  #     files_to_merge = CreateFoFN_samtools.fofn_list,
  #     TSO_name = Tso_name
  # } 


  #call merge_reports as merge_fastp_reports{

    # input:  
    #   files_to_merge = CreateFoFN_fastp.fofn_list,
    #   TSO_name = Tso_name
    # } 



  # call make_excel { 
  #   input:
  #   Tso_name = Tso_name,
  #   pipeline_version = pipeline_v,
  #   tabla1 = merge_fastp_reports.merged_report,
  #   pestana1 = "Filtrado",
  #   tabla2 = merge_samtools_reports.merged_report, 
  #   pestana2 = "Alineamiento",
  #   tabla3 = merge_reports.merged_report,
  #   pestana3 = "Profundidad-en-libreria"

  # }

  # Array[File] reportes_salidas = ["${make_excel.reporte_excel}"]
  # Array[Pair[String,File]] samples_x_files = cross (path_save, reportes_salidas)
  # scatter (pairs in samples_x_files) {
  #   call symlink_important_files {
  #       input:
  #       output_to_save = pairs.right,
  #       path_save = pairs.left
  #   }
  # }  



output {
Array[File] depth_global_cov_stats = bam_depth.glob_cov_stats ###estadistica del alineamiento...
Array[File] by_exon_depth = bam_depth.cov_stats_by_exon

####pdf_report
Array[File] bams_sex_prediction = bam_depth.sex_prediction
  #Array[File] bams_stat_depth_global_coverage_stats = bam_depth.glob_cov_stats
  #Array[File] stat_alineamiento #samtools_reporte_final from bam2gvcf

###test  borrar dps de pdf report.
  #Array[File] bams_stat_depth_global_coverage_stats_out = bam_depth.glob_cov_stats
 # Array[File] bams_sex_prediction_out = bam_depth.sex_prediction
Array[File] fastp_rep_out = fastp_qual.fastp_stats

#File coverage_merged_report = merge_reports.merged_report
#Array[File] reporte_final = samtools_reports_file.output_global_report ### archivo para mergear... estadistica en la libreria del experimento
#File excel_qual_report = make_excel.reporte_excel
#Array[File] Samt_bam_stat = samtools_stat.samtools_stat_original_bam 
#Array[File] Samt_TSO_stat = samtools_stat.samtools_stat_TSO_bam

}

}

