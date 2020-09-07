####quality control. New release

###fastp

task fastp_qual {
  File inputs_json_report
  String report_name = basename(inputs_json_report, ".txt")

  #${sep=' -I ' input_bqsr_reports}
  command <<<
  /home/hnrg/NGStools/pipeline_wdl/qualityControl/estadistica_fastp_V2.py -i ${inputs_json_report} -o ${report_name}_fastp_report.tsv
  >>>

  output {
    File fastp_stats = "${report_name}_fastp_report.tsv"

  }
}



task histo_cob {
    
#input_bam=/data/resultsHNRG/*/CC1707556/CC1707556.bam
        File intervalo_captura
        File input_bam
        Array[File] input_bam_index
        File ensembl2intervalo_captura
        String sample_name = basename( input_bam,'.bam')
        String toolpath
        String ngs_toolpath 
    
    
    command {   
        #!/bin/bash
        set -e
        set -o pipefail

        # esto reporta la cobertura en cada intervalo de captura y hace un histograma global también con el keyword "all"
        ${toolpath}/bedtools2/bin/coverageBed -a ${intervalo_captura} -b ${input_bam}  -hist > ${sample_name}.hist.aux
        echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        cat header.txt ${sample_name}.hist.aux > ${sample_name}.hist 
        rm ${sample_name}.hist.aux header.txt

        #histograma global del bam restringido a toda la librería
        grep '^all' ${sample_name}.hist > global.hist
        echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global.header.txt
        cat global.header.txt global.hist > ${sample_name}.global.hist
        rm global.header.txt global.hist

        ####prediccion de sexo para reporte pdf
        ###ya que estamos, prediccion de sexo:
        ${ngs_toolpath}/python_scripts/bam_sex_xy.py -b ${input_bam} > ${sample_name}_sex.txt


        #### COBERTURA  ##################################
        #### EXONES     ##################################
        #histograma restringido a cada exon de ensembl que está en la librería de captura
        
        ${toolpath}/bedtools2/bin/coverageBed -a ${ensembl2intervalo_captura} -b ${input_bam}  -hist > ${sample_name}.ENS.hist.aux1
        echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        grep -v '^all' ${sample_name}.ENS.hist.aux1 > ${sample_name}.ENS.hist.aux2
        cat header.txt ${sample_name}.ENS.hist.aux2 > ${sample_name}.ENS.hist
        rm ${sample_name}.ENS.hist.aux1 ${sample_name}.ENS.hist.aux2 header.txt
        ##


        

        
    }

    output {
        File histo_exon = "${sample_name}.ENS.hist"
        File histo_global ="${sample_name}.global.hist"
        File sex_prediction = "${sample_name}_sex.txt"

    }

}


task make_tsv_reports {
    
        File by_exon_cov  
        File global_cov
        String ngs_toolpath
        String sample_name
        #String path_save
    

    command {

        #!/bin/bash
        set -e

        # make global tsv report
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/global_coverage_report_inLibrary.py -i=${global_cov} -o ${sample_name}_experiment_global_report.tsv -op ${sample_name}.distributions.eps -s ${sample_name}

        # make tsv coverage report by exon
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/local_coverage_report_ENS_intersect_Library.py -i=${by_exon_cov} -o ${sample_name}_ENS_local_report.tsv -s=${sample_name}


    }

    output {
        File hist_by_exon = "${sample_name}_ENS_local_report.tsv" 
        File hist_global = "${sample_name}_experiment_global_report.tsv"
        File distributions_plot = "${sample_name}.distributions.eps"

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
  String ngs_toolpath

  command{
   ${ngs_toolpath}/pipeline_wdl/qualityControl/make_excel_report.py ${tabla1}:${pestana1} ${tabla2}:${pestana2} ${tabla3}:${pestana3} ${experiment_name}_qual_report.xlsx
 
  }

  output {
    File reporte_excel = "${experiment_name}_qual_report.xlsx"

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
String pipeline_v
String experiment_name
File exon_coords
File intervalo_captura

####inputs from bam2gvcf.reporte_final
Array[File] stat_alineamiento

scatter (fastp in fastp_json_files){
    call fastp_qual {
      input:
      inputs_json_report = fastp
    }
  }

######################scatter por los bams... analysis_readybam

  scatter (bams_ready in analysis_readybam) {

    call histo_cob {
      input: 
      input_bam = bams_ready,
      input_bam_index = analysis_readybam_index,
      #pipeline_version = pipeline_v,
      intervalo_captura = intervalo_captura,
      ensembl2intervalo_captura = exon_coords,
      toolpath = toolpath,
      ngs_toolpath = ngs_toolpath
     
    }

    ###crear tsv 
    call make_tsv_reports {
        input:
        by_exon_cov =   histo_cob.histo_exon,
        global_cov = histo_cob.histo_global,
        ngs_toolpath = ngs_toolpath,
        sample_name = basename(bams_ready, '.bam')
    }
 }


 #reportes global nuevo from ari script
 #Array[File] bams_stat_depth_global_coverage_stats = make_tsv_reports.hist_global

 #report por exon... 
 #Array[File] by_exon_report = make_tsv_reports.hist_by_exon

 ###reportes calidad fastp
 Array[File] fastp_rep = fastp_qual.fastp_stats


 #Create a file with a list of the generated histo glob_stats for merge in excel report
  call CreateFoFN {
    input:
      array_of_files = make_tsv_reports.hist_global,#bams_stat_depth_global_coverage_stats,
      fofn_name = experiment_name
     
  }

 #Create a file with a list of the generated output_global_report
  call CreateFoFN as CreateFoFN_samtools{
    input:
      array_of_files = stat_alineamiento,
      fofn_name = experiment_name
     
  }

 call CreateFoFN as CreateFoFN_fastp{
    input:
      array_of_files = fastp_rep,
      fofn_name = experiment_name
     
  }


 ####### esto mergea archivos de distintas muestras
  call merge_reports {

    input:  
    files_to_merge = CreateFoFN.fofn_list,
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
    tabla3 = merge_reports.merged_report,
    pestana3 = "Profundidad-en-libreria",
    ngs_toolpath = ngs_toolpath

  }

  Array[File] reportes_salidas = ["${make_excel.reporte_excel}"]
  Array[Pair[String,File]] samples_x_files = cross (path_save, reportes_salidas)
  scatter (pairs in samples_x_files) {
    call symlink_important_files {
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
    }
  }  

 output {
        
    Array[File] depth_global_cov_stats = histo_cob.histo_global ###estadistica del alinamiento...
    Array[File] bams_stat_depth_global_coverage_stats = make_tsv_reports.hist_global
    ##plot distribution
    Array[File] plot_distribution = make_tsv_reports.distributions_plot
#tsv_results
    Array[File] tsv_exon = make_tsv_reports.hist_by_exon


    Array[File] by_exon_depth = histo_cob.histo_exon
    File excel_qual_report = make_excel.reporte_excel
    ####for pdf_report
    Array[File] bams_sex_prediction = histo_cob.sex_prediction

    Array[File] fastp_rep_out = fastp_qual.fastp_stats
 }

}
