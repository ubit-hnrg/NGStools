#import '../qualityControl.wdl' as qc 


###########################TASKS
task mkdir {
    String path_softlink

    command{
        mkdir -p ${path_softlink}
    }
}


task mkdir_samplename {
    String path_softlink
    String samplename

    command{
        mkdir -p ${path_softlink}${samplename}
    }

    output {
        String path_out_softlink = "${path_softlink}" + "${samplename}"+"/"
}
}


task copy2data {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}

task coord_generator {

  File intervalo_captura
  File chromosome_length
  File generic_exon_coords
  Int padding
  Int merge_tolerance
  String toolpath
  String gatk_jar
  #File ref_dict
  String path_save
  String library_name = basename(intervalo_captura, ".bed" )

  
   
  command <<<
    #!/bin/bash
    set -e
    set -o pipefail
    
    ${toolpath}bedtools2/bin/slopBed -i ${intervalo_captura} -g ${chromosome_length} -b ${padding} | sort -k1,1 -k2,2n -V > ${library_name}_padded_${padding}.bed 

       
    ####Exon_restricted interval for quality_control  ${library_name}_padded_${padding}.bed | sort -k1,1 -k2,2n -V 
    ${toolpath}bedtools2/bin/intersectBed -a ${generic_exon_coords} -b ${intervalo_captura} > exon_restricted2_${library_name}.bed
  

    cp -L ${intervalo_captura} ${path_save}
    cp -L ${library_name}_padded_${padding}.bed ${path_save}
    #cp -L ${library_name}_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list ${path_save}
    #cp -L ${library_name}_padded_${padding}.interval_list ${path_save}
    cp -L exon_restricted2_${library_name}.bed ${path_save}
  >>>

  output {

    File padded_coord = "${library_name}_padded_${padding}.bed"
    File exon_restricted = "exon_restricted2_${library_name}.bed" ##for quality_control
    #File interval_list = "${library_name}_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list"
    #File eval_interval_list = "${library_name}_padded_${padding}.interval_list"

  }

}

task read_file_of_tabulated_inputs {
  File tabulatedSampleFilePaths
  Array[String] my_lines = read_lines(tabulatedSampleFilePaths)

  command <<<
      cut -f1 -d'|' ${tabulatedSampleFilePaths} >  sample_id.list
      cut -f2 -d'|' ${tabulatedSampleFilePaths} >  sampleNAME.list
      cut -f3 -d'|' ${tabulatedSampleFilePaths} >  R1_fastq.list
      cut -f4 -d'|' ${tabulatedSampleFilePaths} >  R2_fastq.list
      cut -f1 -d'|' ${tabulatedSampleFilePaths}| sort| uniq > unique_sample_id.list

    >>>
        #cut -f1 -d'|' ${tabulatedSampleFilePaths}| sort| uniq >  unique_samples.list
    output {
      File samplenames = 'sample_id.list'
      Array[String] info_name = read_lines('sampleNAME.list')
      Array[String] array_of_samples = read_lines('sample_id.list')
    	Array[File] array_of_R1_files = read_lines('R1_fastq.list')
      Array[File] array_of_R2_files = read_lines('R2_fastq.list')
      File unique_samples = 'unique_sample_id.list'
    }  
}

## cleaning fastq files
task fastp {

  String sample_name
  File R1_fastq_gz
  File R2_fastq_gz
  String R1_stripped_basename = basename(R1_fastq_gz, ".fastq.gz")
  String R2_stripped_basename = basename(R2_fastq_gz, ".fastq.gz")
  String report_name = basename(R2_fastq_gz,"_R2_001.fastq.gz")
  String toolpath
  Int trim_front
  Int trim_tail


  command {
    ${toolpath}fastp -i ${R1_fastq_gz} -I ${R2_fastq_gz} -o ${R1_stripped_basename}_cleaned.fastq.gz -O ${R2_stripped_basename}_cleaned.fastq.gz -h ${report_name}_fastp.html -j ${report_name}_fastp.json --trim_front1=${trim_front} --trim_tail1=${trim_tail} -A
  rm ${R1_stripped_basename}_cleaned.fastq.gz ${R2_stripped_basename}_cleaned.fastq.gz
  
  }

  output {
    #File fastq_cleaned_R1 = "${R1_stripped_basename}_cleaned.fastq.gz"
    #File fastq_cleaned_R2 = "${R2_stripped_basename}_cleaned.fastq.gz"
    File fastp_json_report = "${report_name}_fastp.json"
    File fastp_html_report = "${report_name}_fastp.html"  
  }

}




task CreateFoFN2 {
  # Command parameters
  Array[File] array_of_files
  String fofn_name
  #Map [String, String] lista = {"array_of_files":"fofn_name"} 
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list \
   
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }
}

task symlink_important_files {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}

#####quality control
###fastp

task fastp_qual {
  File inputs_json_report
  String report_name = basename(inputs_json_report, ".txt")

  #${sep=' -I ' input_bqsr_reports}
  command <<<
  /home/hnrg/NGStools/pipeline_wdl/qualityControl/estadistica_fastp_V2.py -i ${inputs_json_report} -o ${report_name}_fastp_report.tsv -bb ${report_name}_N_bases_before.txt -ba ${report_name}_N_bases_after.txt -ra ${report_name}_N_reads_after.txt
  >>>

  output {
    File fastp_stats = "${report_name}_fastp_report.tsv"
    File bases_after = "${report_name}_N_bases_after.txt"
    File bases_before = "${report_name}_N_bases_before.txt"
    File reads_after = "${report_name}_N_reads_after.txt"

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
        String path_save 
    
    
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

        ####samtools stat
        ${toolpath}samtools stats ${input_bam} -t ${intervalo_captura} > ${sample_name}_samtools.stats
         /usr/local/bin/plot-bamstats ${sample_name}_samtools.stats -p ${path_save}samtools_plots/${sample_name}


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
        File samtools_stat_experiment_bam = "${sample_name}_samtools.stats"

    }

}

task samtools_reports_file {

  String sampleID
  String N_total_reads
  String N_bases_before
  String N_bases_after##from fastp_report
  File samtools_library_report
  String ngs_toolpath

  command {
  ${ngs_toolpath}/pipeline_wdl/qualityControl/samtools_stats_report_V2.py -N=${N_total_reads}  -l=${samtools_library_report} -ba ${N_bases_after} -bb ${N_bases_before} -o=${sampleID}_samtools_report.tsv

  }

  output {
 
  File output_global_report = "${sampleID}_samtools_report.tsv" 

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

task Create_inputs_for_preprocesing {
  File bams_sample_names
  File ubams_paths 


  command <<<  
  python <<CODE 

  with open("${bams_sample_names}", "r") as sf:
      samples = sf.readlines()
      samples =[i.strip('\n') for i in samples]
      if samples[-1]=='':
          samples = samples[:-1]
        

  with open("${ubams_paths}", "r") as ubf:
      ubams = ubf.readlines()
      ubams =[i.strip('\n') for i in ubams]
      if ubams[-1]=='':
          ubams = ubams[:-1]
      
  open_files = []
  for i in range(len(samples)):
      sample = samples[i]
      ubam = ubams[i]
    
      filename ='%s.txt'%sample
      if sample not in open_files:
          with open(filename,'w') as f:
            f.write("%s\n"%ubam)
            open_files.append(sample)
          f.close()
      else:
          with open(filename,'a') as f:
            f.write("%s\n"%ubam)
          f.close()

  CODE
  >>>

  output {
    Array[File] ubam_samples = glob("*.txt")
  }
}




task CreateFoFN {
  # Command parameters
  Array[File]+ array_of_files
  String fofn_name
  
  command {
    for f in $array_of_files; do print $f;done
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


workflow upgrade_statistics {

File tabulatedSampleFilePaths ##samples
String path_softlink

File list_bams
Array[File] bams = read_lines(list_bams)
File list_bam_index
Array[File] bams_index = read_lines(list_bam_index)


 ##### parametros trimmeado fastp
  Int trim_front_fastp = "5" 
  Int trim_tail_fastp = "5"


###GATK
  String gatk_jar = "gatk-package-4.0.8.1-local.jar"
  String toolpath = "/home/hnrg/HNRG-pipeline-V0.1/tools/"
  String ngs_toolpath = "/home/hnrg/NGStools"

  ###coordenadas exonicas (usamos ENSEMBL)
  File generic_exon_coords = "/home/hnrg/HNRG-pipeline-V0.1/libraries/intervalos/ensembl_canonicos_GRCh37_0based.tsv"
###################### inputs para crear intervalo
  File intervalo_captura
  File chromosome_length = "/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/chromosome_lengths_hg19.txt"
  Int padding = "100"
  Int merge_tolerance = "200"
  String experiment_name 


###################calls 
    call mkdir {
        input: 
        path_softlink = path_softlink
    }



call coord_generator {
        input:
        intervalo_captura = intervalo_captura,
        chromosome_length = chromosome_length,
        padding = padding,
        merge_tolerance = merge_tolerance,
        toolpath = toolpath,
        gatk_jar = gatk_jar,
        #ref_dict = ref_dict,
        #padded_interval = coord_generator.padded_coord,
        generic_exon_coords = generic_exon_coords,
        path_save = path_softlink
    }

call read_file_of_tabulated_inputs {
      input:
      tabulatedSampleFilePaths = tabulatedSampleFilePaths
  }

  ###Convert multiple pairs of input fastqs in parallel
  scatter (i in range(length(read_file_of_tabulated_inputs.array_of_samples ))) {

    call fastp {
      input: 
      trim_front = trim_front_fastp ,
      trim_tail = trim_tail_fastp,
      sample_name = read_file_of_tabulated_inputs.array_of_samples[i],
      R1_fastq_gz = read_file_of_tabulated_inputs.array_of_R1_files[i],
      R2_fastq_gz = read_file_of_tabulated_inputs.array_of_R2_files[i],
      toolpath = toolpath
    }

} 

   #Create a file with a list of the generated fastp_json file
  call CreateFoFN2 as FoFN_fastp_json {
    input:
    array_of_files = fastp.fastp_json_report,
    fofn_name = "fastp_reports_json" 
  }
 
  

call Create_inputs_for_preprocesing as fastp_report_files {
    input:
    ubams_paths = FoFN_fastp_json.fofn_list,
    bams_sample_names = read_file_of_tabulated_inputs.samplenames
  }

#Array[File] muestras  =  Create_inputs_for_preprocesing.ubam_samples

 Array[File] fastp_json_reports  =  fastp_report_files.ubam_samples
 #Array[File] fastp_html = fastp.fastp_html_report
 
scatter (samples in fastp_json_reports){
call fastp_qual {
      input:
      inputs_json_report = samples#fastp_report_files.ubam_samples
    }
}
## fin scatter fastp
 ####qual control

 Array[File] N_total_reads_bam = fastp_qual.reads_after ###ahora es sobre N_bases
  Array[File] N_bases_after_filtering = fastp_qual.bases_after
 Array[File]  N_bases_before_filtering = fastp_qual.bases_before
 
  

#####################scatter por los bams reducidos
   scatter (idx in range(length(bams))){

      call mkdir_samplename {
    input: 
     path_softlink = path_softlink,
     samplename = basename(bams[idx], '.bam')#sample_name
    }
    call histo_cob {
      input: 
      input_bam = bams[idx],#bams_ready,
      input_bam_index = bams_index,
      #pipeline_version = pipeline_v,
      intervalo_captura = intervalo_captura,
      ensembl2intervalo_captura = coord_generator.exon_restricted, #exon_coords,
      toolpath = toolpath,
      ngs_toolpath = ngs_toolpath,
      path_save = mkdir_samplename.path_out_softlink
     
    }

  call samtools_reports_file {

  input: 
  sampleID = basename(bams[idx], '.bam'),#base_file_name,
  N_total_reads = read_string(N_total_reads_bam[idx]), ###ahora es sobre N_bases
  N_bases_after = read_string(N_bases_after_filtering[idx]),
  N_bases_before = read_string(N_bases_before_filtering[idx]), 

  #samtools_global_report = samtools_stat.samtools_stat_original_bam,
  samtools_library_report = histo_cob.samtools_stat_experiment_bam,
  ngs_toolpath = ngs_toolpath

  }

    ###crear tsv 
    call make_tsv_reports {
        input:
        by_exon_cov =   histo_cob.histo_exon,
        global_cov = histo_cob.histo_global,
        ngs_toolpath = ngs_toolpath,
        sample_name = basename(bams[idx], '.bam') #basename(fastp_qual.fastp_stats[idx],'_fastp_report.tsv')#basename(analysis_readybam[idx], '.bam')
    }
 
  
 }

Array[String] path_save = mkdir_samplename.path_out_softlink

#Create a file with a list of the generated histo glob_stats for merge in excel report
  call CreateFoFN {
    input:
      array_of_files = make_tsv_reports.hist_global,#bams_stat_depth_global_coverage_stats,
      fofn_name = experiment_name,# basename(bams[idx], '.bam')#experiment_name
     
  }

 Array[File] samtools_global = samtools_reports_file.output_global_report
 #Create a file with a list of the generated output_global_report
  call CreateFoFN as CreateFoFN_samtools{
    input:
      array_of_files = samtools_global,#samtools_reports_file.output_global_report,#stat_alineamiento,
      fofn_name = experiment_name #basename(bams[idx], '.bam')#experiment_name
     
  }

 call CreateFoFN as CreateFoFN_fastp{
    input:
      array_of_files = fastp_qual.fastp_stats,#fastp_rep,
      fofn_name = experiment_name # basename(bams[idx], '.bam')#experiment_name
     
  }


 ###### esto mergea archivos de distintas muestras
  call merge_reports {

    input:  
    files_to_merge = CreateFoFN.fofn_list,
    experiment_name = experiment_name,#basename(bams[idx], '.bam'),#experiment_name,
    ngs_toolpath = ngs_toolpath

  } 

  call merge_reports as merge_samtools_reports{

    input:  
      files_to_merge = CreateFoFN_samtools.fofn_list,
      ngs_toolpath = ngs_toolpath,
      experiment_name = experiment_name#basename(bams[idx], '.bam'),#experiment_name
  } 


  call merge_reports as merge_fastp_reports{

    input:  
      files_to_merge = CreateFoFN_fastp.fofn_list,
      experiment_name = experiment_name,#basename(bams[idx], '.bam'),#experiment_name,
      ngs_toolpath = ngs_toolpath
    } 

#  call make_excel { 
#     input:
#     experiment_name = experiment_name,#basename(bams[idx], '.bam'),#experiment_name,
#     tabla1 = merge_fastp_reports.merged_report,
#     pestana1 = "Filtrado",
#     tabla2 = merge_samtools_reports.merged_report, 
#     pestana2 = "Alineamiento",
#     tabla3 = merge_reports.merged_report,
#     pestana3 = "Profundidad-en-libreria",
#     ngs_toolpath = ngs_toolpath

#   }


#  ###samtools_stat
#   #File excel_report = make_excel.reporte_excel

#    Array[File] reportes_salidas = ["${make_excel.reporte_excel}"]
#   Array[Pair[String,File]] samples_x_files = cross (path_save, reportes_salidas)
#   scatter (pairs in samples_x_files) {
#     call symlink_important_files {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }  
#     Array[File]+ hist_glob = make_tsv_reports.hist_global
#    Array[Pair[String,File]] global_out = zip (path_save, hist_glob)
#   scatter (pairs in global_out) {
#     call symlink_important_files as save_global {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }
#       Array[File]+ distrib_plot = make_tsv_reports.distributions_plot
#    Array[Pair[String,File]] distri_out = zip (path_save, distrib_plot)
#   scatter (pairs in distri_out) {
#     call symlink_important_files as save_distri_plot{
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }

#         Array[File]+ hist_exon = make_tsv_reports.hist_by_exon
#    Array[Pair[String,File]] hist_exon_out = zip (path_save, hist_exon)
#   scatter (pairs in hist_exon_out) {
#     call symlink_important_files as save_hist_exon{
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }
  #  Array[Pair[String,File]] fastp_html_out = zip (path_save, fastp_html)
  # scatter (pairs in fastp_html_out) {
  #   call symlink_important_files as save_html_fastp {
  #       input:
  #       output_to_save = pairs.right,
  #       path_save = pairs.left
  #   }
  # }
 




####end
}