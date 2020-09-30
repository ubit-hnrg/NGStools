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
        sort -k1,1V -k2,2n ${intervalo_captura} > intervalo_sorted.bed
        

        # # esto reporta la cobertura en cada intervalo de captura y hace un histograma global también con el keyword "all"
        # #${toolpath}bedtools2/bin/sort ${input_bam} -m 1G | 
        # ${toolpath}bedtools2/bin/coverageBed -a intervalo_sorted.bed -b ${input_bam} -hist ${sorted} > ${sample_name}.hist.aux
        # #${toolpath}bedtools2/bin/coverageBed -a ${intervalo_captura} -b ${input_bam} -sorted -hist > ${sample_name}.hist.aux
        # echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        # cat header.txt ${sample_name}.hist.aux > ${sample_name}.hist 
        # rm ${sample_name}.hist.aux header.txt

        # #histograma global del bam restringido a toda la librería
        # grep '^all' ${sample_name}.hist > global.hist
        # echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global.header.txt
        # cat global.header.txt global.hist > ${sample_name}.global.hist
        # rm global.header.txt global.hist
        
        ###sex prediction
        ${ngs_toolpath}/python_scripts/bam_sex_xy.py -b ${input_bam} > ${sample_name}_sex_${pipeline_version}.txt

        ###septiembre,20: se agrega eliminar duplicados.
        ${toolpath}samtools view -F1024 -u ${input_bam} > bam_nodups.sam
         ${toolpath}samtools stats bam_nodups.sam -t intervalo_sorted.bed > ${sample_name}_nodups.stats_${pipeline_version}.txt_${pipeline_version}.txt_${pipeline_version}.txt_${pipeline_version}.txt_${pipeline_version}.txt

        ####global_hist for no dups_bams
        ${toolpath}bedtools2/bin/coverageBed -a intervalo_sorted.bed -b bam_nodups.sam -hist ${sorted} > ${sample_name}_nodup.hist.aux
        ##${toolpath}bedtools2/bin/coverageBed -a ${intervalo_captura} -b ${input_bam} -sorted -hist > ${sample_name}.hist.aux
        echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header_nodup.txt
        cat header_nodup.txt ${sample_name}_nodup.hist.aux > ${sample_name}_nodup.hist 
        rm ${sample_name}_nodup.hist.aux header_nodup.txt bam_nodups.sam

        #histograma global del bam nodup restringido a toda la librería
        grep '^all' ${sample_name}_nodup.hist > global_nodup.hist
        echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global_nodup.header.txt
        cat global_nodup.header.txt global_nodup.hist > ${sample_name}_nodups.stats_${pipeline_version}.txt
        rm global_nodup.header.txt global_nodup.hist

        
         ####samtools stat
        ${toolpath}samtools stats ${input_bam} -t ${intervalo_captura} > ${sample_name}_samtools_${pipeline_version}.stats
         /usr/local/bin/plot-bamstats ${sample_name}_samtools_${pipeline_version}.stats -p ${path_save}samtools_plots/${sample_name}

         #### COBERTURA  ##################################
        #### EXONES     ##################################
        #histograma restringido a cada exon de ensembl que está en la librería de captura
        
        ##${toolpath}bedtools2/bin/sort ${input_bam} -m 1G | 
        ${toolpath}bedtools2/bin/coverageBed -a ${ensembl2intervalo_captura} -b ${input_bam} -hist ${sorted} > ${sample_name}.ENS_${pipeline_version}.hist.aux1
        echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        grep -v '^all' ${sample_name}.ENS_${pipeline_version}.hist.aux1 > ${sample_name}.ENS_${pipeline_version}.hist.aux2
        cat header.txt ${sample_name}.ENS_${pipeline_version}.hist.aux2 > ${sample_name}.ENS_${pipeline_version}.hist
        rm ${sample_name}.ENS_${pipeline_version}.hist.aux1 ${sample_name}.ENS_${pipeline_version}.hist.aux2 header.txt
        ##
        ##regiones no cubiertas en el intervalo de captura. -bga reporta la profunidad in bedgraph format. reporta las regiones con 0 cobertura. 
        ## por lo que dps se puede filtrar lo no cubierto.-
        bedtools genomecov -ibam ${input_bam} -bga | awk '$4==0'| bedtools intersect -a intervalo_sorted.bed -b - > ${sample_name}.no_cubierto_intervalo_${pipeline_version}.tsv

        
        #cp -L ${sample_name}.global.hist ${path_save}
        cp -L ${sample_name}_samtools_${pipeline_version}.stats ${path_save}
        cp -L ${sample_name}.ENS_${pipeline_version}.hist ${path_save}
        cp -L ${sample_name}_sex_${pipeline_version}.txt ${path_save} 
        cp -L ${sample_name}_nodups.stats_${pipeline_version}.txt ${path_save}
        cp -L ${sample_name}.no_cubierto_intervalo_${pipeline_version}.tsv ${path_save}
       

           
    }

    output {
        #File histo_global ="${sample_name}.global.hist"
        File histo_global_nodup = "${sample_name}_nodups.stats_${pipeline_version}.txt"
        File samtools_stat_experiment_bam = "${sample_name}_samtools_${pipeline_version}.stats"
        File histo_exon = "${sample_name}.ENS_${pipeline_version}.hist"
        File sex_prediction = "${sample_name}_sex_${pipeline_version}.txt"
        File nodups = "{sample_name}_nodups.stats_${pipeline_version}.txt"
        File no_cubierto_intervalo = "${sample_name}.no_cubierto_intervalo_${pipeline_version}.tsv"


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
  File pipeline_version

  command {
  ${ngs_toolpath}/pipeline_wdl/qualityControl/samtools_stats_report_V2.py -N=${N_total_reads}  -l=${samtools_library_report} -d ${samtools_dup} -ba ${N_bases_after} -bb ${N_bases_before} -o=${sampleID}_samtools_report.tsv
  
  cp -L ${sampleID}_samtools_report.tsv ${path_save}

  }

  output {
 
  File output_global_report = "${sampleID}_samtools_report.tsv" 

  }

}


task make_tsv_reports {
    
        File by_exon_cov  
        #File global_cov
        File global_cov_nodups
        String ngs_toolpath
        String sample_name
        String path_save
        File pipeline_version
    
# make global tsv report
        #python ${ngs_toolpath}/pipeline_wdl/qualityControl/global_coverage_report_inLibrary.py -i=${global_cov} -o ${sample_name}_experiment_global_report.tsv -op ${sample_name}.distributions.eps -s ${sample_name}

    command {

        #!/bin/bash
        set -e

        
        # make global_nodups tsv report
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/global_coverage_report_inLibrary.py -i=${global_cov_nodups} -o ${sample_name}_experiment_nodups_global_report.tsv -op ${sample_name}.nodups.distributions.eps -s ${sample_name}

        # make tsv coverage report by exon
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/local_coverage_report_ENS_intersect_Library.py -i=${by_exon_cov} -o ${sample_name}_ENS_local_report.tsv -s=${sample_name}
       
        cp -L  ${sample_name}_ENS_local_report.tsv ${sample_name}_experiment_nodups_global_report.tsv ${sample_name}.nodups.distributions.eps ${path_save}
####${sample_name}.distributions.eps ${sample_name}_experiment_global_report.tsv
    }
 
    output {
        File hist_by_exon = "${sample_name}_ENS_local_report.tsv" 
        #File hist_global = "${sample_name}_experiment_global_report.tsv"
        #File distributions_plot = "${sample_name}.distributions.eps"
        File hist_global_nodups = "${sample_name}_experiment_nodups_global_report.tsv"
        File distributions_plot_nodups = "${sample_name}.nodups.distributions.eps"

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
    File pipeline_version
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
Array[File]+ fastp_json_files
Array[String] path_save
String experiment_path
String pipeline_v
String experiment_name
File exon_coords
File intervalo_captura
#Array[String] bams_N_reads

####inputs from bam2gvcf.reporte_final
#Array[File] stat_alineamiento

scatter (fastp in fastp_json_files){
    call fastp_qual {
      input:
      inputs_json_report = fastp,
      ngs_toolpath = ngs_toolpath
    }
  }


 Array[File] total_reads_fastq = fastp_qual.reads_after ###ahora es sobre N_bases
 Array[File] N_bases_after_filtering_fastq_fastq = fastp_qual.bases_after
 Array[File]  N_bases_before_filtering_fastq_fastq = fastp_qual.bases_before
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
 
       
  call samtools_reports_file {

  input: 
  sampleID = basename(analysis_readybam[idx], '.bam'),#base_file_name,
  # total_reads_fastq = bams_N_reads[idx], ###ahora es sobre N_bases
  # N_bases_after_filtering_fastq_fastq = fastp_qual.bases_after,
  # N_bases_before_filtering_fastq = fastp_qual.bases_before, 
  #N_total_reads = read_string(N_total_reads), ###ahora es sobre N_bases
  #N_bases_after = read_string(N_bases_after),
  #N_bases_before = read_string(N_bases_before),
  samtools_dup = cobertura.nodups,
  #ensembl2intervalo_captura = ensembl2intervalo_captura,
  N_total_reads = read_string(total_reads_fastq[idx]),
  N_bases_before = read_string(N_bases_before_filtering_fastq_fastq[idx]),
  N_bases_after =  read_string(N_bases_after_filtering_fastq_fastq[idx]),


  #samtools_global_report = samtools_stat.samtools_stat_original_bam,
  samtools_library_report = cobertura.samtools_stat_experiment_bam,
  ngs_toolpath = ngs_toolpath,
  path_save = path_save[idx]


  }

    ###crear tsv 
    call make_tsv_reports {
        input:
        by_exon_cov =   cobertura.histo_exon,
        #global_cov = cobertura.histo_global,
        ngs_toolpath = ngs_toolpath,
        sample_name = basename(analysis_readybam[idx], '.bam'),
        path_save = path_save[idx],
        global_cov_nodups = cobertura.histo_global_nodup,
        pipeline_version = pipeline_v


    }
 }


 #reportes global nuevo from ari script
 #Array[File] bams_stat_depth_global_coverage_stats = make_tsv_reports.hist_global

 #report por exon... 
 #Array[File] by_exon_report = make_tsv_reports.hist_by_exon

 ###reportes calidad fastp
 Array[File] fastp_rep = fastp_qual.fastp_stats
 
 ###File hist_by_exon = "${sample_name}_ENS_local_report.tsv"  ##este es output y va al main para armar el excel x exones.


 #Create a file with a list of the generated histo glob_stats for merge in excel report
  # call CreateFoFN {
  #   input:
  #     array_of_files = make_tsv_reports.hist_global,#bams_stat_depth_global_coverage_stats,
  #     fofn_name = experiment_name
     
  # }

  call CreateFoFN as CreateFoFN_no_dups {
    input:
      array_of_files = make_tsv_reports.hist_global_nodups,#bams_stat_depth_global_coverage_stats,
      fofn_name = experiment_name
     
  }

 #Create a file with a list of the generated output_global_report
  call CreateFoFN as CreateFoFN_samtools{
    input:
      array_of_files = samtools_reports_file.output_global_report,#stat_alineamiento,
      fofn_name = experiment_name
     
  }

 call CreateFoFN as CreateFoFN_fastp{
    input:
      array_of_files = fastp_rep,
      fofn_name = experiment_name
     
  }


 ####### esto mergea archivos de distintas muestras
  # call merge_reports {

  #   input:  
  #   files_to_merge = CreateFoFN.fofn_list,
  #   experiment_name = experiment_name,
  #   ngs_toolpath = ngs_toolpath

  # } 

   call merge_reports as merge_nodups_report{

    input:  
    files_to_merge = CreateFoFN_no_dups.fofn_list,
    experiment_name = experiment_name,
    ngs_toolpath = ngs_toolpath

  } 


  call merge_reports as merge_samtools_reports{

    input:  
      files_to_merge = CreateFoFN_samtools.fofn_list,
      ngs_toolpath = ngs_toolpath,
      experiment_name = experiment_name
  } 


  call merge_reports as merge_fastp_reports{

    input:  
      files_to_merge = CreateFoFN_fastp.fofn_list,
      experiment_name = experiment_name,
      ngs_toolpath = ngs_toolpath
    } 

 call make_excel { 
    input:
    experiment_name = experiment_name,
    tabla1 = merge_fastp_reports.merged_report,
    pestana1 = "Filtrado",
    tabla2 = merge_samtools_reports.merged_report, 
    pestana2 = "Alineamiento",
    #tabla3 = merge_reports.merged_report,
    #pestana3 = "Profundidad-en-libreria",
    tabla3= merge_nodups_report.merged_report,
    pestana3 = "Profundidad-en-libreria_1024",
    ngs_toolpath = ngs_toolpath

  }

################solo guardo un excel de calidad x experimento en carpeta /data/resulthsHNRG/experiment_name/excel_qual_xlsx

call symlink_important_files as save_excel_qual {
         input:
         output_to_save = make_excel.reporte_excel,
         path_save = experiment_path
     }

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
    Array[File] depth_global_cov_stats = cobertura.histo_global_nodup
    Array[File] bams_stat_depth_global_coverage_stats = make_tsv_reports.hist_global_nodups
    ##plot distribution
    Array[File] plot_distribution = make_tsv_reports.distributions_plot_nodups
#tsv_results
    Array[File] tsv_exon = make_tsv_reports.hist_by_exon
    Array[File] nocubierto = cobertura.no_cubierto_intervalo


    #Array[File] by_exon_depth = cobertura.histo_exon
    File excel_qual_report = make_excel.reporte_excel
    ####for pdf_report
    Array[File] bams_sex_prediction = cobertura.sex_prediction

    Array[File] fastp_rep_out = fastp_qual.fastp_stats
    
    ##fram samtools_stats
    Array[File] reporte_final = samtools_reports_file.output_global_report ### archivo para mergear... estadistica en la libreria del experimento
    Array[File] samtools_stat_report_from_reduced_bam = cobertura.samtools_stat_experiment_bam

 
 }

}
