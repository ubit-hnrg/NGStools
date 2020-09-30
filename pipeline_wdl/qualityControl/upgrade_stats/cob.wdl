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
        ${ngs_toolpath}/python_scripts/bam_sex_xy.py -b ${input_bam} > ${sample_name}_sex.txt

        ###septiembre,20: se agrega eliminar duplicados.
        ${toolpath}samtools view -F1024 -u ${input_bam} > bam_nodups.sam
         ${toolpath}samtools stats bam_nodups.sam -t intervalo_sorted.bed > ${sample_name}_nodups.stats

        ####global_hist for no dups_bams
        ${toolpath}bedtools2/bin/coverageBed -a intervalo_sorted.bed -b bam_nodups.sam -hist ${sorted} > ${sample_name}_nodup.hist.aux
        ##${toolpath}bedtools2/bin/coverageBed -a ${intervalo_captura} -b ${input_bam} -sorted -hist > ${sample_name}.hist.aux
        echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header_nodup.txt
        cat header_nodup.txt ${sample_name}_nodup.hist.aux > ${sample_name}_nodup.hist 
        rm ${sample_name}_nodup.hist.aux header_nodup.txt bam_nodups.sam

        #histograma global del bam nodup restringido a toda la librería
        grep '^all' ${sample_name}_nodup.hist > global_nodup.hist
        echo -e 'chr\tDP\tBPs\tIntervalLength\tfrequency' > global_nodup.header.txt
        cat global_nodup.header.txt global_nodup.hist > ${sample_name}_nodup.global.hist
        rm global_nodup.header.txt global_nodup.hist

        
         ####samtools stat
        ${toolpath}samtools stats ${input_bam} -t ${intervalo_captura} > ${sample_name}_samtools.stats
         /usr/local/bin/plot-bamstats ${sample_name}_samtools.stats -p ${path_save}samtools_plots/${sample_name}

         #### COBERTURA  ##################################
        #### EXONES     ##################################
        #histograma restringido a cada exon de ensembl que está en la librería de captura
        
        ##${toolpath}bedtools2/bin/sort ${input_bam} -m 1G | 
        ${toolpath}bedtools2/bin/coverageBed -a ${ensembl2intervalo_captura} -b ${input_bam} -hist ${sorted} > ${sample_name}.ENS.hist.aux1
        echo -e 'chr\tstart\tend\ttranscriptID\tgene\texonNumber\tstrand\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
        grep -v '^all' ${sample_name}.ENS.hist.aux1 > ${sample_name}.ENS.hist.aux2
        cat header.txt ${sample_name}.ENS.hist.aux2 > ${sample_name}.ENS.hist
        rm ${sample_name}.ENS.hist.aux1 ${sample_name}.ENS.hist.aux2 header.txt
        ##
        ##regiones no cubiertas en el intervalo de captura. -bga reporta la profunidad in bedgraph format. reporta las regiones con 0 cobertura. 
        ## por lo que dps se puede filtrar lo no cubierto.-
        bedtools genomecov -ibam ${input_bam} -bga | awk '$4==0'| bedtools intersect -a intervalo_sorted.bed -b - > ${sample_name}.no_cubierto_intervalo.tsv

        
        #cp -L ${sample_name}.global.hist ${path_save}
        cp -L ${sample_name}_samtools.stats ${path_save}
        cp -L ${sample_name}.ENS.hist ${path_save}
        cp -L ${sample_name}_sex.txt ${path_save} 
        cp -L ${sample_name}_nodups.stats ${path_save}
        cp -L ${sample_name}.no_cubierto_intervalo.tsv ${path_save}
       

           
    }

    output {
        #File histo_global ="${sample_name}.global.hist"
        File histo_global_nodup = "${sample_name}_nodup.global.hist"
        File samtools_stat_experiment_bam = "${sample_name}_samtools.stats"
        File histo_exon = "${sample_name}.ENS.hist"
        File sex_prediction = "${sample_name}_sex.txt"
        File nodups = "${sample_name}_nodups.stats"
        File no_cubierto_intervalo = "${sample_name}.no_cubierto_intervalo.tsv"


    }
 
 } ###fin



task samtools_reports_file {

  String sampleID
  String N_total_reads
  String N_bases_before
  String N_bases_after ##from fastp_report
  File samtools_library_report
  File samtools_dup
  String path_save
  String ngs_toolpath

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
   

  call samtools_reports_file {

  input: 
  sampleID = basename(bam_in, '.bam'),#base_file_name,
  N_total_reads = read_string(N_total_reads), ###ahora es sobre N_bases
  N_bases_after = read_string(N_bases_after),
  N_bases_before = read_string(N_bases_before),
  samtools_dup = cobertura.nodups,

  #samtools_global_report = samtools_stat.samtools_stat_original_bam,
  samtools_library_report = cobertura.samtools_stat_experiment_bam,
  ngs_toolpath = ngs_toolpath,
  path_save = path_save

  }

    ###crear tsv 
    call make_tsv_reports {
    input:
        by_exon_cov = cobertura.hist_exon,
        #global_cov = cobertura.histo_global,
        global_cov_nodups = cobertura.histo_global_nodup,
        ngs_toolpath = ngs_toolpath,
        sample_name = basename(bam_in, '.bam'),
        path_save = path_save #basename(fastp_qual.fastp_stats[idx],'_fastp_report.tsv')#basename(analysis_readybam[idx], '.bam')
    }
 
 output {
 File por_exon_tsv = make_tsv_reports.hist_by_exon
 File global_tsv = make_tsv_reports.hist_global_nodups
 File distribution_plot = make_tsv_reports.distributions_plot
 File samtools_global = samtools_reports_file.output_global_report
 File nocubierto_intervalo = cobertura.no_cubierto_intervalo
 }


}