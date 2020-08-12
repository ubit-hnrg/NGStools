##### V2 pipeline UIT Agosto 2020
## genotipado individual
## reporte en pdf
## 


import './fastq2ubam_V2.wdl' as fastq2ubam #1
import './bam2gvcf_V2.wdl' as bamtogvcf #2
import './ubam2bwa_V2.wdl' as ubam2bwa #3
import './jointgenotype_single_V2.wdl' as single_genotypeGVCF #4
import './quality_control_V2.wdl' as qual_control #5
import './anotaciones_hnrg_single.wdl' as anotacionesSingle #5




###########################TASKS
##creo q vuela este task por el de abajo
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

  File experiment_lib
  File chromosome_length
  File generic_exon_coords
  Int padding
  Int merge_tolerance
  String toolpath
  String gatk_jar
  File ref_dict
  String path_save

  
   
  command <<<
    #!/bin/bash
    set -e
    set -o pipefail
    
    ${toolpath}bedtools2/bin/slopBed -i ${experiment_lib} -g ${chromosome_length} -b ${padding} | sort -k1,1 -k2,2n -V > intervalo_b37_padded_${padding}.bed 

    ###merged
     
    ${toolpath}bedtools2/bin/mergeBed -i intervalo_b37_padded_${padding}.bed -d ${merge_tolerance} > intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed

    java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed -O=intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list -SD=${ref_dict}  

    java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}.bed -O=intervalo_b37_padded_${padding}.interval_list -SD=${ref_dict}
     
    ####TSO_restricted for quality_control
    ${toolpath}bedtools2/bin/intersectBed -wa -a ${generic_exon_coords} -b intervalo_b37_padded_${padding}.bed | sort -k1,1 -k2,2n -V | uniq > exon_restricted2interval.bed
 

    cp -L intervalo_b37_padded_${padding}.bed ${path_save}
    cp -L intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list ${path_save}
    cp -L intervalo_b37_padded_${padding}.interval_list ${path_save}
    cp -L exon_restricted2interval.bed ${path_save}
  >>>

  output {
    File padded_coord = "intervalo_b37_padded_${padding}.bed"
    File interval_restricted = "exon_restricted2interval.bed" ##for quality_control
    #File merged_padded_coord = "intervalo_b37_padded_merged_${merge_tolerance}.bed"
    File interval_list = "intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list"
    File eval_interval_list = "intervalo_b37_padded_${padding}.interval_list"
  }

}


# task restrict_to_TSO {
#   File padded_interval
#   File generic_exon_coords
#   String toolpath

#   command <<<
#    #!/bin/bash
#     set -e
#     set -o pipefail

#   ${toolpath}bedtools2/bin/intersectBed -wa -a ${generic_exon_coords} -b ${padded_interval} | sort -k1,1 -k2,2n -V | uniq > exon_restricted2interval.bed
#   >>>

#   output {

#     File interval_restricted = "exon_restricted2interval.bed"


#   }

# }

task build_excell_report{
    File annovar_tsv
    File exon_coverage_report
    #String sample
    String samplename2
    #String original_sample
  
     #/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${sample}.output_xlsx

    command{

       /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${samplename2}_variants.xlsx
   
   }    

    output{
        File excell_report = '${samplename2}_variants.xlsx'
    }
}    

task pdf_report {
  
  File alineamiento 
  String name 
  File glob_rep 
  File sex 
  File fastp_rep 
  String tso 
  String date
  String path

  command {
    /home/hnrg/NGStools/python_scripts/pdf_report_per_sample.py -fq ${fastp_rep} -aq ${alineamiento} -dp ${glob_rep} -s ${sex} -n ${name} -d ${date} -t ${tso} -o ${path}/${name}_qual_report.pdf
  }


}




##################################END TASK

########MAIN

workflow main_workflow {

  ###inputs for fastq2ubam workflows
  
  File tabulatedSampleFilePaths ##samples

  String pipeline_version = "V2.0"


  ####metadata
  String run_date                   
  String library_name 
  String platform_name = "Illumina"
  String platformmodel = "NextSeq-500"
  String sequencing_center =  "UIT-HNRG" 
  String readlenght 
  String ubam_list_name = "ubamfiles"
  String ref_name = ".b37"

  ###GATK
  String gatk_jar = "gatk-package-4.0.8.1-local.jar"
  String toolpath = "/home/hnrg/HNRG-pipeline-V0.1/tools/"

  ###coordenadas exonicas (usamos ENSEMBL)
  File generic_exon_coords = "/home/hnrg/HNRG-pipeline-V0.1/libraries/intervalos/intersect_ensembl_Tso.bed"
  
  ###save location
  String path_softlink

 
  ########################  command line para el bwa ##################
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 4 -Y"
  
  ########## referencia
  File ref_fasta = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.ann"
  File ref_fasta_index = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.fai"
  File ref_dict = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.dict"
  File ref_amb = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.amb"
  File ref_ann = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.ann"
  File ref_bwt = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.bwt"
  File ref_pac = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.pac"
  File ref_sa = "/home/hnrg/HNRG-pipeline-V0.1/references/hs37d5/hs37d5.fa.sa"
    
  
  File dbSNP_vcf = "/home/hnrg/HNRG-pipeline-V0.1/dbs/preprocessing_dbs/All_20180423.vcf.gz"
  File dbSNP_vcf_index = "/home/hnrg/HNRG-pipeline-V0.1/dbs/preprocessing_dbs/All_20180423.vcf.gz.tbi"
  ###bam2gvcf input
  Array[File] known_indels_sites_VCFs = ["/home/hnrg/HNRG-pipeline-V0.1/dbs/preprocessing_dbs/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"]
  Array[File] known_indels_sites_indices = ["/home/hnrg/HNRG-pipeline-V0.1/dbs/preprocessing_dbs/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi"]

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
  String reference_version = "GRCh37.75"


  ##annovar
  String db_annovar = "/home/hnrg/HNRG-pipeline-V0.1/dbs/hg19_annovar/" #path annovar
    
    
    ###esto creo q vuela
    File annovar_table_pl #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
    File joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py


  ###################### inputs para crear intervalo
  File experiment_lib
  File chromosome_length = "/home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/chromosome_lengths_hg19.txt"
  Int padding = "100"
  Int merge_tolerance = "200"


    ###################calls 
    call mkdir {
        input: 
        path_softlink = path_softlink
    }

    #### se generan los intervalos para los distintos pasos con padding y merge_tolerance

    call coord_generator {
        input:
        experiment_lib = experiment_lib,
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

   call bamtogvcf.bam2gvcf {
     input:
      base_file_name =  sample_name,
      lib_resctricted = coord_generator.padded_coord,
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
      wes_calling_interval_list = coord_generator.interval_list,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      haplotype_scatter_count = haplotype_scatter_count,
      compression_level = compression_level,
      gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation,
      gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
      wgs_calling_interval_list = coord_generator.interval_list, 
      wgs_evaluation_interval_list = coord_generator.interval_list,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      smith_waterman_implementation = smith_waterman_implementation,
      contamination = contamination,
      newqual = newqual,
      java_heap_memory_initial = java_heap_memory_initial,
      tso_bed = coord_generator.padded_coord
    } 

    ######single_genotype
    call single_genotypeGVCF.singleGenotypeGVCFs {
        input:
        #num_gvcfs= cantidad_gvcf,
        eval_interval_list   = coord_generator.eval_interval_list,
        array_path_save = mkdir_samplename.path_out_softlink,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        callset_name = basename(tabulatedSampleFilePaths, ".txt"),
        ref_fasta = ref_fasta,
        ref_fasta_index =ref_fasta_index,
        ref_dict = ref_dict,
        gatk_jar = gatk_jar,
        toolpath = toolpath,
        sample_name = sample_name,
        #input_gvcfs = gvcf.left,
        #input_gvcfs_indices = gvcf.right,
        region_padded_bed = coord_generator.padded_coord,

        ####input del anterior jointgenotyping
        input_gvcf = bam2gvcf.output_gvcf,
        input_gvcf_index = bam2gvcf.output_gvcf_index,
        pipeline_v = pipeline_version,

        ####annovar
        db_annovar = db_annovar,#path annovar
        annovar_table_pl = annovar_table_pl, #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
        joinPY = joinPY
        
    }



    call anotacionesSingle.FuncionalAnnotationSingle {
        input:
        input_vcf = singleGenotypeGVCFs.restricted_vcf,
        path_save = mkdir_samplename.path_out_softlink,
        toolpath = toolpath,
        samplename1 = sample_name,
        java_heap_memory_initial = "12g",
        pipeline_v = pipeline_version,
        exon_coordinates = coord_generator.interval_restricted,
        reference_version = reference_version
        

      }

   } ###fin scatter gvcf



  ### intervalo para quality control -> coord_generator.interval_restricted

 Array[File] archivos_a_borrar3 = bam2gvcf.borrar_SortandFix #,"${}"]

  scatter (archivos in archivos_a_borrar3){
    call borrado as borrado_Sort_and_Fix {
    input:
      archivo_borrar = archivos
    }
  } 

 #Array[String] uniquesample_name =read_lines(ConvertPairedFastQsToUnmappedBamWf.samplesnames)

#Array[File] salidas_json = ConvertPairedFastQsToUnmappedBamWf.fastp_json_reports
 Array[String] array_path_save_json = mkdir_samplename.path_out_softlink
 #Array[Pair[String,File]] samples_x_files_json = zip (array_path_save_json, salidas_json)
#  scatter (pairs in samples_x_files_json) {
  #    call symlink_important_files {
  #      input:
  #       output_to_save = pairs.right,
  #       path_save = pairs.left
  #     }
  # }


call qual_control.quality_control_V2 {
   input: 
   stat_alineamiento = bam2gvcf.reporte_final,
   fastp_json_files = ConvertPairedFastQsToUnmappedBamWf.fastp_json_reports,
   path_save = mkdir_samplename.path_out_softlink,
   analysis_readybam = bam2gvcf.analysis_ready_bam,
   analysis_readybam_index = bam2gvcf.analysis_ready_bam_index,
   toolpath = toolpath,
   Tso_name = basename(tabulatedSampleFilePaths, ".txt"),
   exon_coords = coord_generator.interval_restricted,
   pipeline_v = pipeline_version
   #tso_bed = tso_bed
  }





#  Array[Pair[String,File]] test_save1 = zip (array_path_save_byexon, test1)
#   scatter (pairs in test_save1) {
#     call symlink_important_files as test{
#         input:
#         output_to_save = pairs.right,
#         path_save = pairs.left
#     }
#   }



 
 ##descomentar abajo y borrar la de arriba
 Array[File] prof_by_exon = quality_control_V2.by_exon_depth##","${coord_generator.padded_coord}"] #"${name}_coverage_statistics_by_exon.tsv"
 Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink

 Array[Pair[String,File]] samples_by_exon = zip (array_path_save_byexon, prof_by_exon)
  scatter (pairs in samples_by_exon) {
    call symlink_important_files as byexon{
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
    }
  }

### poner el exon_distance con el vcf anotado_hnrg freq
## anotar con plof, gnomad /home/hnrg/HNRG-pipeline-V0.1/libraries/GRCh37/gnomad_plof_HNRG.tsv
####for pdf purpose 
  Array[File] alineamiento_rep = bam2gvcf.reporte_final ### archivo para mergear... estadistica en la libreria del experimento
  Array[File] global = quality_control_V2.depth_global_cov_stats
  #Array[File] fastp_qual = quality_control_V2.fastp_rep_out
  Array[File] sex_pred= quality_control_V2.bams_sex_prediction
 
#  scatter (indice in range(length(alineamiento_rep))){
       
#        String S1 = basename(alineamiento_rep[indice],"_samtools_report.tsv")
#        String S2 = basename(prof_by_exon[indice],"_global_coverage_statistics_V2.0.tsv")
#        String S3 = basename(fastp[indice],"_fastp_report.tsv")
#        String S4 = basename(sex_pred[indice],"_sex.txt")

 ####excel_report
    Array[File] Tsv_annovar = singleGenotypeGVCFs.annovar_tsv_out
    scatter (idx in range(length(Tsv_annovar))){
       
       String sample = basename(Tsv_annovar[idx],"multianno_restrict.tsv")
       String samplename2 = basename(prof_by_exon[idx],"_coverage_statistics_by_exon_V2.0.tsv")
       
       if(sample==samplename2){
       call build_excell_report {
            input:
            annovar_tsv = Tsv_annovar[idx],
            samplename2 = samplename2,
            exon_coverage_report = prof_by_exon[idx]
            
           }
        call pdf_report {
            input:
            alineamiento = Tsv_annovar[idx],
            name = samplename2,
            glob_rep = global[idx],
            sex = sex_pred[idx],
            #fastp_rep = fastp_qual[idx],
            fastp_rep = quality_control_V2.fastp_rep_out[idx],
            tso = basename(tabulatedSampleFilePaths, ".txt"),
            date = run_date,
            path = array_path_save_json[idx]
            
           }  
          }

    }

Array[File?] reporte_variantes = build_excell_report.excell_report
#Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink
 Array[Pair[String,File?]] samples_by_variant = zip (array_path_save_byexon, reporte_variantes)
 # scatter (pairs in samples_by_variant) {
 #   call symlink_important_files as build_excell_reportbyvariants{
 #       input:
 #       output_to_save = pairs.right,
 #       path_save = pairs.left
 #   }
 # }





   # Outputs that will be retained when execution is complete
  output {
    Array[File] output_ubams = ConvertPairedFastQsToUnmappedBamWf.output_ubams
    Array[String] output_ubams_sample_names =  ConvertPairedFastQsToUnmappedBamWf.output_ubams_sample_names
    File unmapped_ubam_list = ConvertPairedFastQsToUnmappedBamWf.unmapped_ubam_list
    File samplesnames = ConvertPairedFastQsToUnmappedBamWf.samplesnames
    Array[File] muestras  =  ConvertPairedFastQsToUnmappedBamWf.muestras
  # Outputs del workflow bam2gvcf (bam2gvcf)  
 
   Array[File] duplication_metrics = bam2gvcf.duplication_metrics
   Array[File] bqsr_report = bam2gvcf.bqsr_report 
   Array[File] analysis_ready_bam = bam2gvcf.analysis_ready_bam
   Array[File] analysis_ready_bam_index = bam2gvcf.analysis_ready_bam_index
   Array[File] analysis_ready_bam_md5 = bam2gvcf.analysis_ready_bam_md5 
   Array[File] gvcf_summary_metrics = bam2gvcf.gvcf_summary_metrics 
   Array[File] gvcf_detail_metrics = bam2gvcf.gvcf_detail_metrics 
   Array[File] output_vcf = bam2gvcf.output_gvcf
   Array[File] output_vcf_index = bam2gvcf.output_gvcf_index
  ##### output jointgenotype
   
   ###output_joint_single_vcf
 
   Array[File?] outputvcf = singleGenotypeGVCFs.outputvcf
   Array[File?] outputvcfindex =  singleGenotypeGVCFs.outputvcfindex
   Array[File?] detail_metrics_file =  singleGenotypeGVCFs.metrica1
   Array[File?] summary_metrics_file = singleGenotypeGVCFs.metrica2
   #File? intervalo = JointGenotyping.inter

   #Array[File] Samt_bam_stat = bam2gvcf.Samt_bam_stat 
   Array[File] Samt_TSO_stat = bam2gvcf.Samt_TSO_stat
   Array[File] reporte_final = bam2gvcf.reporte_final ### archivo para mergear... estadistica en la libreria del experimento
 }

}




  