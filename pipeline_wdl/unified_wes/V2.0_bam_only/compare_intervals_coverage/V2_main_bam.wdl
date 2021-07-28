#####uit julio 2021
## pipeline para generar bam sin restricciones y comparar cobertura intervalos.
## 


import './fastq2ubam_V2.wdl' as fastq2ubam #1
import './only_bam.wdl' as onlybam #3
import './ubam2bwa_V2.wdl' as ubam2bwa #2
import './qc_bam.wdl' as qual_control #4




###########################TASKS
task mkdir {
    String path_softlink

    command{
        mkdir -p ${path_softlink}
    }
}

task borrado_fastp {
    File path1
    File path2

    command <<<
    while read -r line; do
    rm "$line"
    done < "${path1}"
    while read -r line2; do
    rm "$line2"
    done < "${path2}"
    >>>

}

task borrado_bam_intermedio {
    File path_bams_intermedios
    

    command <<<
    while read -r line; do
    rm "$line"
    done < "${path_bams_intermedios}"
    >>>

}

task borrado {
String archivo_borrar

command <<<
readlink -f ${archivo_borrar} | xargs rm
>>>
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


task groupingBams_bysample_glob {
  String pattern 
  Array[File] array_of_files

  command{
  }

  output {
    Array[File] subArray_input_ubam2gvcf = glob('../inputs/*/*${pattern}*')
  }
}

task symlink_important_files {
    File output_to_save
    String path_save
    #ln -s ${output_to_save} ${path_save}
    command{
       cp -L ${output_to_save} ${path_save}
    }
}

task coord_generator {

  File intervalo_captura
  ###File intervalo_captura2
  ###File intervalo_captura3
  File chromosome_length
  File generic_exon_coords  ##ensembl
  Int padding
  Int merge_tolerance
  String toolpath
  String gatk_jar
  File ref_dict
  String path_save
  String library_name = basename(intervalo_captura, ".bed" )

  
   
  command <<<
    #!/bin/bash
    set -e
    set -o pipefail
    
    ${toolpath}bedtools2/bin/slopBed -i ${intervalo_captura} -g ${chromosome_length} -b ${padding} | sort -k1,1 -k2,2n -V > ${library_name}_padded_${padding}.bed 

     ####Exon_restricted interval for quality_control  ${library_name}_padded_${padding}.bed | sort -k1,1 -k2,2n -V 
    ${toolpath}bedtools2/bin/intersectBed -a ${generic_exon_coords} -b ${intervalo_captura}|sort -k1,1V -k2,2n - > exon_restricted2_${library_name}.bed

    
    ###merged
     
    ${toolpath}bedtools2/bin/mergeBed -i ${library_name}_padded_${padding}.bed -d ${merge_tolerance} > ${library_name}_padded_${padding}_merged_${merge_tolerance}.bed

    java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=${library_name}_padded_${padding}_merged_${merge_tolerance}.bed -O=${library_name}_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list -SD=${ref_dict}

    java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=${library_name}_padded_${padding}.bed -O=${library_name}_padded_${padding}.interval_list -SD=${ref_dict}
     
  
    cp -L ${intervalo_captura} ${path_save}
    cp -L ${library_name}_padded_${padding}.bed ${path_save}

    cp -L ${library_name}_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list ${path_save}
    cp -L ${library_name}_padded_${padding}.interval_list ${path_save}
    cp -L exon_restricted2_${library_name}.bed ${path_save}
  >>>

  output {
    File padded_coord = "${library_name}_padded_${padding}.bed"
    File exon_restricted = "exon_restricted2_${library_name}.bed" ##for quality_control

    #File interval_restricted = "exon_restricted2interval.bed" ##for quality_control
    #File merged_padded_coord = "intervalo_b37_padded_merged_${merge_tolerance}.bed"
    File interval_list = "${library_name}_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list"
    File eval_interval_list = "${library_name}_padded_${padding}.interval_list"
  }

}



task build_excell_report{
    #File annovar_tsv
    File exon_coverage_report
    #File plof
    String samplename2
    String ngs_toolpath
    #File no_cubierto
    String pipeline_version

    #String original_sample
  
     #/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${sample}.output_xlsx
     #este anda        ${ngs_toolpath}/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${plof}:GnomAD_PLOF ${no_cubierto}:no_cubierto  ${samplename2}_variants_${pipeline_version}.xlsx

    command{

       ${ngs_toolpath}/pipeline_wdl/qualityControl/make_excel_report.py ${exon_coverage_report}:ExonCoverage ${samplename2}_variants_${pipeline_version}.xlsx
   
   }    

    output{
        File excell_report = '${samplename2}_variants_${pipeline_version}.xlsx'
    }
}    


##################################END TASK

########MAIN

workflow main_workflow {

  ###inputs for fastq2ubam workflows
  
  File tabulatedSampleFilePaths ##samples

  String? pipeline_version = "V3" ###no me lo toma el input como default


  ####metadata
  String run_date                   
  String library_name 
  String platform_name = "Illumina"
  String platformmodel = "NextSeq-500"
  String sequencing_center =  "UIT-HNRG" 
  String readlenght 
  String ubam_list_name = "ubamfiles"
  String ref_name = ".hg38"

  ###GATK
  String gatk_jar = "gatk-package-4.0.8.1-local.jar"
  String toolpath = "/home/hnrg/HNRG-pipeline-V0.1/tools/"
  String ngs_toolpath = "/home/hnrg/NGStools"


  ###coordenadas exonicas (usamos ENSEMBL)
  File generic_exon_coords = " '/home/hnrg/HNRG-pipeline-V0.1/libraries/intervalos/hg38/ensembl_transcriptos_canoninos_nobiotype_104_38.tsv'"
  
  ###save location
  String path_softlink

 
  ########################  command line para el bwa ##################
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 4 -Y"
  
  ########## referencia
  File ref_fasta = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
  File ref_fasta_index = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai"
  File ref_dict = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dict"
  File ref_amb = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.amb"
  File ref_ann = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.ann"
  File ref_bwt = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.bwt"
  File ref_pac = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.pac"
  File ref_sa = "/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.64.sa"
    
  
  File dbSNP_vcf = "/data/new_dbs/grch38/dbSNP/GCF_000001405.39.gz"
  File dbSNP_vcf_index = "/data/new_dbs/grch38/dbSNP/GCF_000001405.39.gz.tbi"
  ###bam2gvcf input
  Array[File] known_indels_sites_VCFs = ["/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz","/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"]
  Array[File] known_indels_sites_indices = ["/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi","/data/new_dbs/grch38/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"]

  ##### parametros trimmeado fastp
  Int trim_front_fastp = "5" 
  Int trim_tail_fastp = "5"
    
  Int break_bands_at_multiples_of = "1000000"
  Int haplotype_scatter_count = "2"

  ##################################
  Int compression_level = "1"
  String java_heap_memory_initial = "128m"

  ########optimization flags
  String gatk_gkl_pairhmm_implementation = "LOGLESS_CACHING"
  Int gatk_gkl_pairhmm_threads = "1"
    
  ####optimizacion haplotypecaller
  String smith_waterman_implementation = "AVX_ENABLED"
  Float? contamination = "0"
  String newqual = "true"


  ##################anotacion funcional
  String reference_version = "GRCh38"


  ##annovar
  #String db_annovar = "/home/hnrg/HNRG-pipeline-V0.1/dbs/hg19_annovar/" #path annovar
    
    
    ###esto creo q vuela
  #  File annovar_table_pl #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
  #  File joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py


  ###################### inputs para crear intervalo
  File intervalo_captura
  File chromosome_length = "/data/new_dbs/grch38/hg38/chromosome_lengths_GHRCH38.txt" #"/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/chromosome_lengths_hg19.txt"
  Int padding = "100"
  Int merge_tolerance = "200" ##usamos 5000


    ###################calls 
    call mkdir {
        input: 
        path_softlink = path_softlink
    }

    #### se generan los intervalos para los distintos pasos con padding y merge_tolerance

    call coord_generator {
        input:
        intervalo_captura = intervalo_captura,
        chromosome_length = chromosome_length,
        padding = padding,
        merge_tolerance = merge_tolerance,
        toolpath = toolpath,
        gatk_jar = gatk_jar,
        ref_dict = ref_dict,
        #padded_interval = coord_generator.padded_coord,
        generic_exon_coords = generic_exon_coords,
        path_save = path_softlink
    }

    # ## se limita el intervalo a la libreria del experimento que se va a usar para el control de calidad.
   
    # call restrict_to_TSO {
    #     input:
    #     padded_interval = coord_generator.padded_coord,
    #     generic_exon_coords = generic_exon_coords,
    #     toolpath = toolpath
    # }


    ###raw_data in fastq to uBAM
    call fastq2ubam.ConvertPairedFastQsToUnmappedBamWf {  
        input: 
        trim_front = trim_front_fastp, 
        trim_tail = trim_tail_fastp,
        tabulatedSampleFilePaths = tabulatedSampleFilePaths,
        run_date = run_date,                   
        library_name = library_name,   
        platform_name = platform_name,
        gatk_jar = gatk_jar,
        toolpath = toolpath,
        ubam_list_name = ubam_list_name,
        sequencing_center = sequencing_center,
        read_lenght = readlenght,
        path_softlink = path_softlink,
        platform_model = platformmodel
    }
 
    ##alineamiento y mapeo a la referencia.
    call ubam2bwa.ubamtobwa {
        input:
        array_unmapped_bams = ConvertPairedFastQsToUnmappedBamWf.output_ubams,
        ref_name = ref_name,
        compression_level = compression_level,
        java_heap_memory_initial = java_heap_memory_initial,
        bwa_commandline = bwa_commandline,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        gatk_jar = gatk_jar,
        toolpath = toolpath
    } 

  Array[File] array_of_samples_txt = ConvertPairedFastQsToUnmappedBamWf.muestras

  ##borrado de archivos de fastp para liberar espacio
    call borrado_fastp {
        input:
        path1 = ConvertPairedFastQsToUnmappedBamWf.p_borrar1,
        path2 = ConvertPairedFastQsToUnmappedBamWf.p_borrar2
    }

    call borrado_bam_intermedio {
    input:
     path_bams_intermedios = ubamtobwa.output_path_bam_intermedio
    }
 #inputs_bams is an array of files. Each element is a file containing all the aligned and merged bams of a sample.
 scatter (sample_txt in array_of_samples_txt)  {

   #Array[File] flowcell_mapped_bams = read_lines(flowcell_mapped_bams_listfile)

   String sample_name = basename(sample_txt, ".txt")

   #####subworkflow de fastq2bwa

   call mkdir_samplename {
    input: 
     path_softlink = path_softlink,
     samplename = sample_name
    }

   call groupingBams_bysample_glob {
     input: 
       pattern = sample_name,
       array_of_files = ubamtobwa.output_mergedbam_files

    }

   call onlybam.only_bam_wflw {
     input:
      
      base_file_name =  sample_name,
      lib_restricted = coord_generator.padded_coord,
      path_save = mkdir_samplename.path_out_softlink,
      bams_entrada = groupingBams_bysample_glob.subArray_input_ubam2gvcf,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      wes_calling_interval_list = coord_generator.eval_interval_list,#coord_generator.interval_list,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      haplotype_scatter_count = haplotype_scatter_count,
      compression_level = compression_level,
      gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation,
      gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
      wgs_calling_interval_list = coord_generator.eval_interval_list,#coord_generator.interval_list, 
      wgs_evaluation_interval_list = coord_generator.eval_interval_list,# coord_generator.interval_list,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      ngs_toolpath = ngs_toolpath,

      smith_waterman_implementation = smith_waterman_implementation,
      contamination = contamination,
      newqual = newqual,
      java_heap_memory_initial = java_heap_memory_initial,
      intervalo_captura = intervalo_captura    
      } 
    }
  

 Array[String] uniquesample_name =read_lines(ConvertPairedFastQsToUnmappedBamWf.samplesnames)

Array[File] salidas_json = ConvertPairedFastQsToUnmappedBamWf.fastp_json_reports
 Array[String] array_path_save_json = mkdir_samplename.path_out_softlink
 Array[Pair[String,File]] samples_x_files_json = zip (array_path_save_json, salidas_json)
 scatter (pairs in samples_x_files_json) {
     call symlink_important_files {
       input:
        output_to_save = pairs.right,
        path_save = pairs.left
      }
  }


#call qual_control.qual_control {
#   input: 
#   #stat_alineamiento = only_bam_wflw.reporte_final,
#   #fastp_json_files = ConvertPairedFastQsToUnmappedBamWf.fastp_json_reports,
#   path_save = mkdir_samplename.path_out_softlink,
#   #bams_N_reads = N_reads_bams,#only_bam_wflw.bams_N_reads,
  #  analysis_readybam = only_bam_wflw.analysis_ready_bam,
  #  analysis_readybam_index = only_bam_wflw.analysis_ready_bam_index,
  #  toolpath = toolpath,
  #  ngs_toolpath = ngs_toolpath,
  #  intervalo_captura = intervalo_captura,
  #  pipeline_v= pipeline_version,
  #  experiment_name = basename(tabulatedSampleFilePaths, ".txt"),
  #  exon_coords = generic_exon_coords,#####ensembl crudo para ver que se pierde    ###coord_generator.exon_restricted, #### ensembl vs intervalo_captura
  #  experiment_path = path_softlink
  # }


# Array[File] exon_tsv = qual_control.tsv_exon
#  Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink
#  Array[Pair[String,File]] samples_by_exon = zip (array_path_save_byexon, exon_tsv)
#   scatter (pairs in samples_by_exon) {
#     call symlink_important_files as byexon{
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }


 
 # Array[File] prof_by_exon = quality_control_V2.by_exon_depth##","${coord_generator.padded_coord}"] #"${name}_coverage_statistics_by_exon.tsv"
#Array[File] html_reports_from_fastq = ConvertPairedFastQsToUnmappedBamWf.fastp_html_reports 
#Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink

# Array[Pair[String,File]] html_fastp_by_samplre = zip (array_path_save_byexon, html_reports_from_fastq)
#   scatter (pairs in html_fastp_by_samplre) {
#      call symlink_important_files as html_fastp{
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }


# ####for pdf purpose 
#   Array[File] alineamiento_rep = qual_control.reporte_final_alineamiento#only_bam_wflw.reporte_final ### archivo para mergear... estadistica en la libreria del experimento
#   Array[File] global = qual_control.depth_global_cov_stats 
#   Array[Pair[String,File]] global_report = zip (array_path_save_byexon, global)
#   scatter (pairs in global_report) {
#     call symlink_important_files as global_hist {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }


  
  # ####meter en pdf
  # Array[File] plot_dist = qual_control.plot_distribution
  # Array[Pair[String,File]] plot_hist_save = zip (array_path_save_byexon, plot_dist)
  # scatter (pairs in plot_hist_save) {
  #   call symlink_important_files as save_plot_distrib {
  #       input:
  #       output_to_save = pairs.right,
  #       path_save = pairs.left
  #   }
  # }

#  ###samtools_stat
#   Array[File] samtools_stat_report_from_reduced_bam = qual_control.samtools_stat_report_from_reduced_bam
#    Array[Pair[String,File]] samtools_stat_out = zip (array_path_save_byexon, samtools_stat_report_from_reduced_bam)
#   scatter (pairs in samtools_stat_out) {
#     call symlink_important_files as save_samtools {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }

#  Array[File] tsv_global = qual_control.bams_stat_depth_global_coverage_stats
#   Array[Pair[String,File]] qual_tsv = zip (array_path_save_byexon, tsv_global)
#   scatter (pairs in qual_tsv) {
#     call symlink_important_files as qual_tsv_save {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }



#   #Array[File] fastp_qual = quality_control_V2.fastp_rep_out
#   Array[File] sex_pred= qual_control.bams_sex_prediction
#   Array[File] gene_list = singleGenotypeGVCFs.annovar_gene_list 
#   Array[File] plof = singleGenotypeGVCFs.gene_plof_file 
#   Array[File] exon_distances = FuncionalAnnotationSingle.vcf_exon_distance
#   Array[File] no_cubierto = qual_control.nocubierto



 ####excel_report
#     scatter (idx in range(length(exon_tsv))){
#          String samplename2 = basename(exon_tsv[idx],"_ENS_local_report.tsv")
 
#          call build_excell_report {
#             input:
#             #annovar_tsv = join_annovar_exon_dist.anno_dist,
#             #plof = plof[idx],
#             samplename2 = samplename2,
#             exon_coverage_report = exon_tsv[idx],
#             ngs_toolpath = ngs_toolpath,
#             #no_cubierto = no_cubierto[idx],
#             pipeline_version = pipeline_version
            
#            }
# } #fin scatter


# Array[File?] reporte_variantes = build_excell_report.excell_report
# #Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink
#  Array[Pair[String,File?]] samples_by_variant = zip (array_path_save_byexon, reporte_variantes)
#   scatter (pairs in samples_by_variant) {
#     call symlink_important_files as build_excell_reportbyvariants {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }


#     Array[File] Tsv_annovar = singleGenotypeGVCFs.annovar_tsv_out
#     scatter (idx in range(length(Tsv_annovar))){
       
#        String sample = basename(Tsv_annovar[idx],"multianno_restrict.tsv")
#        String samplename2 = basename(exon_tsv[idx],"_ENS_local_report.tsv")
       
#        #if(sample==samplename2){
#          #mergear tsv_annovar con distancias_exones
        
#       call join_annovar_exon_dist {
#           input:
#             name = samplename2,
#               annovar_variants = Tsv_annovar[idx],
#               exon_dist = exon_distances[idx],
#               ngs_toolpath = ngs_toolpath,
#               pipeline_version = pipeline_version
#          }

#        call build_excell_report {
#             input:
#             annovar_tsv = join_annovar_exon_dist.anno_dist,
#             plof = plof[idx],
#             samplename2 = samplename2,
#             exon_coverage_report = exon_tsv[idx],
#             ngs_toolpath = ngs_toolpath,
#             no_cubierto = no_cubierto[idx],
#             pipeline_version = pipeline_version
            
#            }
#         call pdf_report {
#             input:
#             alineamiento = alineamiento_rep[idx], ##ok
#             name = samplename2,
#             glob_rep =tsv_global[idx], ##ok
#             sex = sex_pred[idx],
#             #fastp_rep = fastp_qual[idx],
#             fastp_rep = qual_control.fastp_rep_out[idx],
#             tso = basename(tabulatedSampleFilePaths, ".txt"),
#             date = run_date,
#             path = array_path_save_json[idx],
#             ngs_toolpath = ngs_toolpath,
#             pipeline_version = pipeline_version
            
#            }  
#         #}

#     }

# Array[File?] reporte_variantes = build_excell_report.excell_report
# #Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink
#  Array[Pair[String,File?]] samples_by_variant = zip (array_path_save_byexon, reporte_variantes)
#   scatter (pairs in samples_by_variant) {
#     call symlink_important_files as build_excell_reportbyvariants {
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }





   # Outputs that will be retained when execution is complete
#   output {
#     Array[File] output_ubams = ConvertPairedFastQsToUnmappedBamWf.output_ubams
#     Array[String] output_ubams_sample_names =  ConvertPairedFastQsToUnmappedBamWf.output_ubams_sample_names
#     File unmapped_ubam_list = ConvertPairedFastQsToUnmappedBamWf.unmapped_ubam_list
#     File samplesnames = ConvertPairedFastQsToUnmappedBamWf.samplesnames
#     Array[File] muestras  =  ConvertPairedFastQsToUnmappedBamWf.muestras
#   # Outputs del workflow only_bam_wflw (only_bam_wflw)  
 
#    Array[File] duplication_metrics = only_bam_wflw.duplication_metrics
#    Array[File] bqsr_report = only_bam_wflw.bqsr_report 
#    Array[File] analysis_ready_bam = only_bam_wflw.analysis_ready_bam
#    Array[File] analysis_ready_bam_index = only_bam_wflw.analysis_ready_bam_index
#    Array[File] analysis_ready_bam_md5 = only_bam_wflw.analysis_ready_bam_md5 
#    Array[File] gvcf_summary_metrics = only_bam_wflw.gvcf_summary_metrics 
#    Array[File] gvcf_detail_metrics = only_bam_wflw.gvcf_detail_metrics 
#    Array[File] output_vcf = only_bam_wflw.output_gvcf
#    Array[File] output_vcf_index = only_bam_wflw.output_gvcf_index
#   ##### output jointgenotype
  
#    ###output_joint_single_vcf
 
#    Array[File?] outputvcf = singleGenotypeGVCFs.outputvcf
#    Array[File?] outputvcfindex =  singleGenotypeGVCFs.outputvcfindex
#    Array[File?] detail_metrics_file =  singleGenotypeGVCFs.metrica1
#    Array[File?] summary_metrics_file = singleGenotypeGVCFs.metrica2
#    #File? intervalo = JointGenotyping.inter
#    Array[File] reporte_final = qual_control.reporte_final ### archivo para mergear... estadistica en la libreria del experimento
#  }

}




  