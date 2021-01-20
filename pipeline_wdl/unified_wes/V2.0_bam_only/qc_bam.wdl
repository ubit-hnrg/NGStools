####quality control. New release

###fastp

task fastp_qual {
  File inputs_json_report
  String report_name = basename(inputs_json_report, ".txt")
  String ngs_toolpath

  #${sep=' -I ' input_bqsr_reports}
  command <<<
  ${ngs_toolpath}/pipeline_wdl/qualityControl/estadistica_fastp_V2.py -i ${inputs_json_report} -o ${report_name}_fastp_report.tsv -bb ${report_name}_N_bases_before.txt -ba ${report_name}_N_bases_after.txt -ra ${report_name}_N_reads_after.txt

  >>>

  output {
    File fastp_stats = "${report_name}_fastp_report.tsv"
    File bases_after = "${report_name}_N_bases_after.txt"
    File bases_before = "${report_name}_N_bases_before.txt"
    File reads_after = "${report_name}_N_reads_after.txt"

  }
}


task cobertura {
    
        File intervalo_captura
        File input_bam
        File input_bam_index
        File ensembl2intervalo_captura
        String sample_name = basename( input_bam,'.bam')
        String toolpath
        String ngs_toolpath
        String path_save
        String pipeline_version
        String? sorted 
    
    
    command {   
        #!/bin/bash
        set -e
        set -o pipefail

        #sort del intervalo de caputura para que funcione bien el coverageBED.
        
        #### COBERTURA  ##################################
        #### EXONES     ##################################
        #histograma restringido a cada exon de ensembl que está en la librería de captura
        
        
        samtools view -uf 0x2 ${input_bam} > proper_pair.bam
        ${toolpath}bedtools2/bin/coverageBed -b proper_pair.bam -a ${ensembl2intervalo_captura} -hist ${sorted} -g /home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5_sizes.genome > ${sample_name}.ENS_${pipeline_version}.hist.aux1
        echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        grep -v '^all' ${sample_name}.ENS_${pipeline_version}.hist.aux1 > ${sample_name}.ENS_${pipeline_version}.hist.aux2
        cat header.txt ${sample_name}.ENS_${pipeline_version}.hist.aux2 > ${sample_name}.ENS_${pipeline_version}.hist
        rm ${sample_name}.ENS_${pipeline_version}.hist.aux1 ${sample_name}.ENS_${pipeline_version}.hist.aux2 header.txt proper_pair.bam
        ##
        ##regiones no cubiertas en el intervalo de captura. -bga reporta la profunidad in bedgraph format. reporta las regiones con 0 cobertura. 
        ## por lo que dps se puede filtrar lo no cubierto.-
        
        ##intersect con el TSO




        cp -L ${sample_name}.ENS_${pipeline_version}.hist ${path_save}
       

           
    }
    #${sample_name}_${pipeline_version}.txt

    output {
        File histo_exon = "${sample_name}.ENS_${pipeline_version}.hist"

    }

}
#   task sex_pred {

#     File input_bam
#     File input_bam_index
#     String sample_name = basename( input_bam,'.bam')
#     String ngs_toolpath
#     String path_save 
       
#     command <<<
#       #!/bin/bash
#       set -e
#       set -o pipefail

#       ####prediccion de sexo para reporte pdf
#       ###ya que estamos, prediccion de sexo:
#       ${ngs_toolpath}/python_scripts/bam_sex_xy.py -b ${input_bam} > ${sample_name}_sex_${pipeline_version}.txt
#       cp -L ${sample_name}_sex_${pipeline_version}.txt ${path_save} 
      
#     >>>

#   output {
#     File sex_prediction = "${sample_name}_sex_${pipeline_version}.txt"
#     }
#     }


task samtools_reports_file {

  String sampleID
  String N_total_reads
  String N_bases_before
  String N_bases_after ##from fastp_report
  File samtools_library_report
  File samtools_dup
  String path_save
  String ngs_toolpath
  String pipeline_version

  command {
  ${ngs_toolpath}/pipeline_wdl/qualityControl/samtools_stats_report_V2.py -N=${N_total_reads} -l=${samtools_library_report} -d ${samtools_dup} -ba ${N_bases_after} -bb ${N_bases_before} -o=${sampleID}_samtools_report.tsv
  
  cp -L ${sampleID}_samtools_report.tsv ${path_save}

  }

  output {
 
  File output_global_report = "${sampleID}_samtools_report.tsv" 

  }

}


task make_tsv_reports {
    
        File by_exon_cov  
        #File global_cov
        #File global_cov_nodups
        String ngs_toolpath
        String sample_name
        String path_save
        String pipeline_version
    
# make global tsv report
        #python ${ngs_toolpath}/pipeline_wdl/qualityControl/global_coverage_report_inLibrary.py -i=${global_cov} -o ${sample_name}_experiment_global_report.tsv -op ${sample_name}.distributions.eps -s ${sample_name}
#local_coverage_report_ENS_intersect_Library.py
    command {

        #!/bin/bash
        set -e

      
        # make tsv coverage report by exon
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/local_coverage_report_ENS_intersect_Library_full_bam.py -i=${by_exon_cov} -o ${sample_name}_ENS_local_report.tsv -s=${sample_name}
       
        cp -L  ${sample_name}_ENS_local_report.tsv ${path_save}
    }
 
    output {
        File hist_by_exon = "${sample_name}_ENS_local_report.tsv" 
        #File hist_global = "${sample_name}_experiment_global_report.tsv"
        #File distributions_plot = "${sample_name}.distributions.eps"
        #File hist_global_nodups = "${sample_name}_experiment_nodups_global_report.tsv"
        #File distributions_plot_nodups = "${sample_name}_nodups_distributions.eps"

    }

}

task merge_reports {

####inputs del paso1 
File files_to_merge
String experiment_name
String ngs_toolpath

command<<<
${ngs_toolpath}/pipeline_wdl/qualityControl/merge_sample_reports.py -i ${files_to_merge} -o ${experiment_name}.merged_report
>>>

output {
File merged_report = "${experiment_name}.merged_report"

}

}

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

task intersect_bam_tso {
String sampleID
String path_save
File reporte_cob_tsv
File intervalo_TSO
String ngs_toolpath
String toolpath

command {
${toolpath}bedtools2/bin/intersectBed -a ${reporte_cob_tsv} -b ${intervalo_TSO} -v > ${sampleID}_offtarget_raw.tsv

python ${ngs_toolpath}/python_scripts/coverage_report_filter.py -r ${sampleID}_offtarget_raw.tsv -l ${intervalo_TSO} -o1 ${sampleID}_exones_offtarget_tso.tsv -o2 ${sampleID}_offtarget_notin_tso.tsv

cp -L ${sampleID}_offtarget_raw.tsv ${path_save}
cp -L ${sampleID}_exones_offtarget_tso.tsv ${path_save}
cp -L ${sampleID}_offtarget_notin_tso.tsv ${path_save}

}

}

task make_excel {
    String experiment_name
    File tabla1
    String pestana1
    File tabla2
    String pestana2
    File tabla3
    String pestana3
    #File tabla4
    #String pestana4
    String ngs_toolpath
    String pipeline_version
    ##${tabla3}:${pestana3}
    command{
    ${ngs_toolpath}/pipeline_wdl/qualityControl/make_excel_report.py ${tabla1}:${pestana1} ${tabla2}:${pestana2} ${tabla3}:${pestana3} ${experiment_name}_qual_report_${pipeline_version}.xlsx
 
  }

  output {
    File reporte_excel = "${experiment_name}_qual_report_${pipeline_version}.xlsx"

  }

}

task symlink_important_files {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}



workflow qual_control{

Array[File] analysis_readybam
Array[File] analysis_readybam_index
String toolpath
String ngs_toolpath
#Array[File]+ fastp_json_files
Array[String] path_save
String experiment_path
String pipeline_v
String experiment_name
File exon_coords
File intervalo_captura
#Array[String] bams_N_reads

####inputs from bam2gvcf.reporte_final
#Array[File] stat_alineamiento

# scatter (fastp in fastp_json_files){
#     call fastp_qual {
#       input:
#       inputs_json_report = fastp,
#       ngs_toolpath = ngs_toolpath
#     }
#   }


#  Array[File] total_reads_fastq = fastp_qual.reads_after ###ahora es sobre N_bases
#  Array[File] N_bases_after_filtering_fastq_fastq = fastp_qual.bases_after
#  Array[File]  N_bases_before_filtering_fastq_fastq = fastp_qual.bases_before
######################scatter por los bams... analysis_readybam

  #scatter (bams_ready in analysis_readybam) {
   scatter (idx in range(length(analysis_readybam))){
    call cobertura {
      input: 
      input_bam = analysis_readybam[idx],#bams_ready,
      input_bam_index = analysis_readybam_index[idx],
      pipeline_version = pipeline_v,
      intervalo_captura = intervalo_captura,
      ensembl2intervalo_captura = exon_coords,
      toolpath = toolpath,
      ngs_toolpath = ngs_toolpath,
      path_save = path_save[idx]
    }
 

    ###crear tsv 
  call make_tsv_reports {
      input:
      by_exon_cov = cobertura.histo_exon,
      #global_cov = cobertura.histo_global,
      ngs_toolpath = ngs_toolpath,
      sample_name = basename(analysis_readybam[idx], '.bam'),
      path_save = path_save[idx],
      #global_cov_nodups = cobertura.histo_global_nodup,
      pipeline_version = pipeline_v
       
      }

    call intersect_bam_tso {
      input: 
      sampleID = basename(analysis_readybam[idx], '.bam'),
      path_save = path_save[idx] ,
      reporte_cob_tsv = make_tsv_reports.hist_by_exon,
      intervalo_TSO = intervalo_captura,
      ngs_toolpath = ngs_toolpath,
      toolpath = toolpath
 }
}


#   call merge_reports as merge_fastp_reports{

#     input:  
#       files_to_merge = CreateFoFN_fastp.fofn_list,
#       experiment_name = experiment_name,
#       ngs_toolpath = ngs_toolpath
#     } 

#  call make_excel { 
#     input:
#     experiment_name = experiment_name,
#     tabla1 = merge_fastp_reports.merged_report,
#     pestana1 = "Filtrado",
#     tabla2 = merge_samtools_reports.merged_report, 
#     pestana2 = "Alineamiento",
#     #tabla3 = merge_reports.merged_report,
#     #pestana3 = "Profundidad-en-libreria",
#     tabla3= merge_nodups_report.merged_report,
#     pestana3 = "Profundidad-en-libreria",
#     ngs_toolpath = ngs_toolpath,
#     pipeline_version = pipeline_v

#   }

################solo guardo un excel de calidad x experimento en carpeta /data/resulthsHNRG/experiment_name/excel_qual_xlsx

# call symlink_important_files as save_excel_qual {
#          input:
#          output_to_save = make_excel.reporte_excel,
#          path_save = experiment_path
#      }

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
        
    #Array[File] depth_global_cov_stats = cobertura.histo_global ###estadistica del alinamiento...
    #Array[File] depth_global_cov_stats = cobertura.samtools_stat_nodup_experiment_bam#histo_global_nodup
    #Array[File] bams_stat_depth_global_coverage_stats = make_tsv_reports.hist_global_nodups
    ##plot distribution
    #Array[File] plot_distribution = make_tsv_reports.distributions_plot_nodups
#tsv_results
    Array[File] tsv_exon = make_tsv_reports.hist_by_exon
    #Array[File] nocubierto = cobertura.no_cubierto_intervalo


    #Array[File] by_exon_depth = cobertura.histo_exon
    #File excel_qual_report = make_excel.reporte_excel
    ####for pdf_report
    #Array[File] bams_sex_prediction = cobertura.sex_prediction

    #Array[File] fastp_rep_out = fastp_qual.fastp_stats
    
    ##fram samtools_stats
    #Array[File] reporte_final_alineamiento = samtools_reports_file.output_global_report ### archivo para mergear... estadistica en la libreria del experimento
    #Array[File] samtools_stat_report_from_reduced_bam = cobertura.samtools_stat_experiment_bam

 
 }

}
