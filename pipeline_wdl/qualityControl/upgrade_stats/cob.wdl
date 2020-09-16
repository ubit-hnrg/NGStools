#####quality control workflow



task cobertura {
    
        File intervalo_captura
        File input_bam
        File input_bam_index
        String sample_name = basename( input_bam,'.bam')
        String toolpath
        String path_save
        String ngs_toolpath
        File ensembl2intervalo_captura
 
    
    command {   
        #!/bin/bash
        set -e
        set -o pipefail

        sort -k1,1V -k2,2n ${intervalo_captura} > intervalo_sorted.bed
        # esto reporta la cobertura en cada intervalo de captura y hace un histograma global también con el keyword "all"
        ${toolpath}/bedtools2/bin/sort ${input_bam} -m 1G | ${toolpath}/bedtools2/bin/coverageBed -a intervalo_sorted.bed -b - -sorted -hist > ${sample_name}.hist.aux
        echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        cat header.txt ${sample_name}.hist.aux > ${sample_name}.hist 
        rm ${sample_name}.hist.aux header.txt

        #histograma global del bam restringido a toda la librería
        grep '^all' ${sample_name}.hist > global.hist
        echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global.header.txt
        cat global.header.txt global.hist > ${sample_name}.global.hist
        rm global.header.txt global.hist
        
        ###sex prediction
        #${ngs_toolpath}/python_scripts/bam_sex_xy.py -b ${input_bam} > ${sample_name}_sex.txt
        
         ####samtools stat
        ${toolpath}samtools stats ${input_bam} -t ${intervalo_captura} > ${sample_name}_samtools.stats
         /usr/local/bin/plot-bamstats ${sample_name}_samtools.stats -p ${path_save}samtools_plots/${sample_name}
        
        #### COBERTURA  ##################################
        #### EXONES     ##################################
        #histograma restringido a cada exon de ensembl que está en la librería de captura
        
        ${toolpath}/bedtools2/bin/sort ${input_bam} -m 1G | ${toolpath}/bedtools2/bin/coverageBed -a ${ensembl2intervalo_captura} -b - -sorted -hist > ${sample_name}.ENS.hist.aux1
        echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        grep -v '^all' ${sample_name}.ENS.hist.aux1 > ${sample_name}.ENS.hist.aux2
        cat header.txt ${sample_name}.ENS.hist.aux2 > ${sample_name}.ENS.hist
        rm ${sample_name}.ENS.hist.aux1 ${sample_name}.ENS.hist.aux2 header.txt
        ##
        


        cp -L ${sample_name}.global.hist ${path_save}
        cp -L ${sample_name}_samtools.stats ${path_save}
        cp -L ${sample_name}.ENS.hist ${path_save}
        #cp -L ${sample_name}_sex.txt ${path_save} 

        }
    output {
        File histo_global ="${sample_name}.global.hist"
        File samtools_stat_experiment_bam = "${sample_name}_samtools.stats"
         File hist_exon = "${sample_name}.ENS.hist"
        #File sex_prediction = "${sample_name}_sex.txt"

    }
 
 } ###fin


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
#       ${ngs_toolpath}/python_scripts/bam_sex_xy.py -b ${input_bam} > ${sample_name}_sex.txt
#       cp -L ${sample_name}_sex.txt ${path_save} 
      
#     >>>

#   output {
#     File sex_prediction = "${sample_name}_sex.txt"
#     }
#     }
    
#     task samtools_experiment_stat {
#       File intervalo_captura
#         File input_bam
#         File input_bam_index
#         String sample_name = basename( input_bam,'.bam')
#         String toolpath
#         String path_save 
    
#     command <<<
#         #!/bin/bash
#         set -e
#         set -o pipefail

#         ####samtools stat
#         ${toolpath}samtools stats ${input_bam} -t ${intervalo_captura} > ${sample_name}_samtools.stats
#          /usr/local/bin/plot-bamstats ${sample_name}_samtools.stats -p ${path_save}samtools_plots/${sample_name}
#         cp -L ${sample_name}_samtools.stats ${path_save}
#         >>>

#         output {
        
#         File samtools_stat_experiment_bam = "${sample_name}_samtools.stats"
#         }
# }
# task cob_exones {

#       File input_bam
#       File input_bam_index
#       File ensembl2intervalo_captura
#       String sample_name = basename( input_bam,'.bam')
#       String toolpath
#       String path_save 

#     command <<<
#         #!/bin/bash
#         set -e
#         set -o pipefail

#         #### COBERTURA  ##################################
#         #### EXONES     ##################################
#         #histograma restringido a cada exon de ensembl que está en la librería de captura
        
#         ${toolpath}/bedtools2/bin/coverageBed -a ${ensembl2intervalo_captura} -b ${input_bam}  -hist > ${sample_name}.ENS.hist.aux1
#         echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
#         grep -v '^all' ${sample_name}.ENS.hist.aux1 > ${sample_name}.ENS.hist.aux2
#         cat header.txt ${sample_name}.ENS.hist.aux2 > ${sample_name}.ENS.hist
#         rm ${sample_name}.ENS.hist.aux1 ${sample_name}.ENS.hist.aux2 header.txt
#         ##
#         cp -L ${sample_name}.ENS.hist ${path_save}
#       >>>

#     runtime {
#     memory: "12GB"
#     }

#     output { 
#     File hist_exon = "${sample_name}.ENS.hist"
#     }

# }

task samtools_reports_file {

  String sampleID
  String N_total_reads
  String N_bases_before
  String N_bases_after##from fastp_report
  File samtools_library_report
  String ngs_toolpath
  String path_save

  command {
  ${ngs_toolpath}/pipeline_wdl/qualityControl/samtools_stats_report_V2.py -N=${N_total_reads}  -l=${samtools_library_report} -ba ${N_bases_after} -bb ${N_bases_before} -o=${sampleID}_samtools_report.tsv
   cp -L ${sampleID}_samtools_report.tsv ${path_save}
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
        String path_save
    

    command {

        #!/bin/bash
        set -e

        # make global tsv report
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/global_coverage_report_inLibrary.py -i=${global_cov} -o ${sample_name}_experiment_global_report.tsv -op ${sample_name}.distributions.eps -s ${sample_name}

        # make tsv coverage report by exon
        python ${ngs_toolpath}/pipeline_wdl/qualityControl/local_coverage_report_ENS_intersect_Library.py -i=${by_exon_cov} -o ${sample_name}_ENS_local_report.tsv -s=${sample_name}
        
        cp -L ${sample_name}.distributions.eps ${sample_name}_experiment_global_report.tsv  ${sample_name}_ENS_local_report.tsv ${path_save}

    }

    output {
        File hist_by_exon = "${sample_name}_ENS_local_report.tsv" 
        File hist_global = "${sample_name}_experiment_global_report.tsv"
        File distributions_plot = "${sample_name}.distributions.eps"

    }

}




workflow coverage_qual {

###inputs
File bam_in
File bam_index_in
###paths
String gatk_jar
String toolpath 
String ngs_toolpath 
String path_save

###intervalos
File intervalo_captura
File ensembl2intervalo_captura #coord_generator.exon_restricted


###from fastp_qual
String N_total_reads
String N_bases_before
String N_bases_after


    call cobertura {
      input: 
      input_bam = bam_in,#[idx],#bams_ready,
      input_bam_index = bam_index_in,#[idx],
      #pipeline_version = pipeline_v,
      intervalo_captura = intervalo_captura,
      ensembl2intervalo_captura = ensembl2intervalo_captura, #exon_coords,
      toolpath = toolpath,
      ngs_toolpath = ngs_toolpath,
      path_save = path_save #mkdir_samplename.path_out_softlink
     
    }
    # call samtools_experiment_stat {
    #   input:
       
    #   intervalo_captura = intervalo_captura,
    #   input_bam = bams_in,#[idx],#bams_ready,
    #   input_bam_index = bams_index_in,#[idx],
    #   toolpath = toolpath,
    #   path_save = path_save #mkdir_samplename.path_out_softlink
    # }
    # call sex_pred{
    #   input:
    #   input_bam = bams[idx],#bams_ready,
    #   input_bam_index = bams_index[idx],
    #   ngs_toolpath = ngs_toolpath,
    #   path_save = mkdir_samplename.path_out_softlink
    # }

    # call cob_exones {
    #   input:
    #   input_bam = bams_in,#[idx],#bams_ready,
    #   input_bam_index = bams_index_in,#[idx],
    #   ensembl2intervalo_captura = ensembl2intervalo_captura, #coord_generator.exon_restricted, #exon_coords,
    #   toolpath = toolpath,
    #   path_save = path_save #mkdir_samplename.path_out_softlink
    # }

  call samtools_reports_file {

  input: 
  sampleID = basename(bam_in, '.bam'),#base_file_name,
  N_total_reads = read_string(N_total_reads), ###ahora es sobre N_bases
  N_bases_after = read_string(N_bases_after),
  N_bases_before = read_string(N_bases_before), 

  #samtools_global_report = samtools_stat.samtools_stat_original_bam,
  samtools_library_report = cobertura.samtools_stat_experiment_bam,
  ngs_toolpath = ngs_toolpath,
  path_save = path_save

  }

    ###crear tsv 
    call make_tsv_reports {
    input:
        by_exon_cov = cobertura.hist_exon,
        global_cov = cobertura.histo_global,
        ngs_toolpath = ngs_toolpath,
        sample_name = basename(bam_in, '.bam'),
        path_save = path_save #basename(fastp_qual.fastp_stats[idx],'_fastp_report.tsv')#basename(analysis_readybam[idx], '.bam')
    }
 
 output {
 File por_exon_tsv = make_tsv_reports.hist_by_exon
 File global_tsv = make_tsv_reports.hist_global
 File distribution_plot = make_tsv_reports.distributions_plot
 File samtools_global = samtools_reports_file.output_global_report
 }


}