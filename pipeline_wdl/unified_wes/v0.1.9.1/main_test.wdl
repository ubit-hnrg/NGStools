#### main para unir workflows

import './hnrg-fastq2ubam_test.wdl' as fastq2ubam 
import './bam2gvcf.wdl' as bamtogvcf
import './ubam2bwa.wdl' as ubam2bwa
import './jointgenotype_sinrecalibracion.wdl' as joint_genotype
import './quality_control.wdl' as qual_control 
import './processMultisampleVCF.wdl' as splitVCF

task borrado_fastp {
File path1
File path2

#command <<<
#mv ${write_lines(path)}  ${t}.txt;
command <<<
while read -r line; do
rm "$line"
done < "${path1}"
while read -r line2; do
rm "$line2"
done < "${path2}"
>>>

}

task build_excell_report{
    File annovar_tsv
    File exon_coverage_report
    String sample
    String samplename2
    #String original_sample
  
    
    command<<<
    if [ "${sample}" == "${samplename2}"]; then 

        /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${sample}.variants.xlsx
    fi
    >>>     

    output{
        File excell_report = '${sample}.variants.xlsx'
    }
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
    command{
       ln -s ${output_to_save} ${path_save}
    }
}




workflow main_workflow {

  ###inputs for fastq2ubam workflows
  File tabulatedSampleFilePaths
  String run_date                   
  String library_name 
  #String platform_unit 
  String platform_name 
  String gatk_jar
  String toolpath
  String ubam_list_name
  String sequencing_center =  "UIT-HNRG" 
  String readlenght
  String platformmodel
 
  File exon_coords
  File tso_bed

  #inputs for preprocesing ( ubam2gvcf)

  ####### MUESTRA
  #File sample_name
  String path_softlink
  String ref_name

  #File flowcell_unmapped_bams_list
  #String unmapped_bam_suffix

  #File unmapped_bam
    
  #String base_file_name 

  #Array[File] flowcell_unmapped_bams
    
 
  ########################  command line para el bwa ##################
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 4 -Y"
  
  ########## referencia
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
    
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  ##### parametros trimmeado fastp
  Int trim_front_fastp = "5" 
  Int trim_tail_fastp = "5"
  ####
  ###Scatter interval
  #archivo de listas de cobertura de intervalos 
  #File wgs_calling_interval_list
  #File wgs_evaluation_interval_list
  #File wes_calling_interval_list
    
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


    String db_annovar #path annovar
    File annovar_table_pl #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
    File joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py


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

  call borrado_fastp {
    input:
     path1 = ConvertPairedFastQsToUnmappedBamWf.p_borrar1,
     path2 = ConvertPairedFastQsToUnmappedBamWf.p_borrar2
  }
 
  call ubam2bwa.ubamtobwa {
   input:
    array_unmapped_bams = ConvertPairedFastQsToUnmappedBamWf.output_ubams,
    #flowcell_unmapped_bams_list = sample,
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

 #inputs_bams is an array of files. Each element is a file containing all the aligned and merged bams of a sample.
 scatter (sample_txt in array_of_samples_txt)  {
   #ubamtobwa.output_mergedbam_files
  
   #File flowcell_mapped_bams_listfile = sample_txt
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
      #base_file_name =  ubamtobwa.nombre_base,
      base_file_name =  sample_name,
      lib_resctricted = tso_bed,
                                          #### ARTEFACTO PARA PROBAR EL WORKFLOW ANTERIOR, ESTO HAY QUE TRABAJRLO Y ARTICULARLO BIEN. 
      #inputs_ubams = ConvertPairedFastQsToUnmappedBamWf.muestras,
      #uniquesample_name = ConvertPairedFastQsToUnmappedBamWf.samplesnames,
      path_save = mkdir_samplename.path_out_softlink,
      #sample_name = sample_name,
      #flowcell_unmapped_bams_list = sample,
      bams_entrada = groupingBams_bysample_glob.subArray_input_ubam2gvcf,
      #bams_entrada = ubamtobwa.output_mergedbam_files[1],                   #### ARTEFACTO PARA PROBAR EL WORKFLOW ANTERIOR, ESTO HAY QUE TRABAJRLO Y ARTICULARLO BIEN. 
      ##bams_entrada = flowcell_mapped_bams,  #array of basm corresponding to ONE sample. 
      #ref_name = ref_name,
      #unmapped_bam_suffix = ".bam",
      #bwa_commandline = bwa_commandline,
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
      #wes_calling_interval_list = wes_calling_interval_list,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      haplotype_scatter_count = haplotype_scatter_count,
      compression_level = compression_level,
      gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation,
      gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
      #wgs_calling_interval_list = wgs_calling_interval_list,
      #wgs_evaluation_interval_list = wgs_evaluation_interval_list,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      smith_waterman_implementation = smith_waterman_implementation,
      contamination = contamination,
      newqual = newqual,
      java_heap_memory_initial = java_heap_memory_initial,
      tso_bed = tso_bed
    } 

    ####control de calidad... y reducir bams

    ##bam2gvcf.analysis_ready_bam ####analysis_ready_bam



  }

 #Array[File] archivos_a_borrar1 = ubamtobwa.output_mergedbam_files#,"${bam2gvcf.borrar_SortandFix}"]

 #scatter (archivos in archivos_a_borrar1){
 #call borrado as borrado_merged_bam {
 #input:
 #archivo_borrar = archivos
 #}
 #}



 #Array[File] archivos_a_borrar2 = groupingBams_bysample_glob.subArray_input_ubam2gvcf#,"${bam2gvcf.borrar_SortandFix}"]

 #scatter (archivos in archivos_a_borrar2){
 #call borrado as borrado_subset_glob {
 #input:
 #archivo_borrar = archivos
 #}
 #}


 
 #Array[File] paths2symb = ["${pepe.dummy_outfile}","${juan.dummy_outfile}"]
 #scatter (paths in paths2symb) {
 #    call symlink_a_lo_guapo {
 #        input:
 #         output_to_save = paths 
 #    }
 #}


 #############borrado de archivos intermedios##################
 ##### borro archivo intermedio mergebam
 #scatter (paths2 in ubam2gvcf.path_mergebam){
 #  call borrado_single {
 #   input: 
 #    path1 = paths2
 #}
 #}

 ##### borro archivos intermedios  mark_dup####sort_and_fix#####gatherbams#########

 #Array[File] paths_borrar = ["${ubam2gvcf.path_markdup}","${ubam2gvcf.path_sortandfix}"]
 #scatter (paths in paths_borrar) {
 #    call borrado_single as borrado_single1 {
 #    input: 
 #    path1 = paths
 #} 
 #}

 #### borro archivo intermedio applybqsr
 #scatter (paths3 in ubam2gvcf.path_applybqsr){
 #  call borrado_single as borrado_single2 {
 #   input: 
 #    path1 = paths3
 #}
 #}


 #Array[File] archivos_a_borrar = bam2gvcf.borrar_Markdup#,"${bam2gvcf.borrar_SortandFix}"]

 #scatter (archivos in archivos_a_borrar){
 #call borrado as borrado_markdup {
 #input:
 #archivo_borrar = archivos
 #}
 #}

 Array[File] archivos_a_borrar3 = bam2gvcf.borrar_SortandFix#,"${}"]

 scatter (archivos in archivos_a_borrar3){
    call borrado as borrado_Sort_and_Fix {
    input:
      archivo_borrar = archivos
    }
  } 

 Array[String] uniquesample_name =read_lines(ConvertPairedFastQsToUnmappedBamWf.samplesnames)



 ##########################################llamada workflow Joint Genotyping
  call joint_genotype.JointGenotyping {
   input:
     array_path_save = mkdir_samplename.path_out_softlink,
     dbSNP_vcf = dbSNP_vcf,
     dbSNP_vcf_index = dbSNP_vcf_index,
     callset_name = basename(tabulatedSampleFilePaths, ".txt"),
     ref_fasta = ref_fasta,
     ref_fasta_index =ref_fasta_index,
     ref_dict = ref_dict,
     gatk_jar = gatk_jar,
     toolpath = toolpath,
     sample_names = uniquesample_name,
     input_gvcfs = bam2gvcf.output_vcf,
     input_gvcfs_indices = bam2gvcf.output_vcf_index
  }


   call splitVCF.processJointVCF {
     input:
     multisampleVCF = JointGenotyping.outputvcf,
     array_path_save = mkdir_samplename.path_out_softlink,

     toolpath = toolpath, 
     region_padded_bed = tso_bed,##Tso_bed
     path_softlink = path_softlink,

    # for annovar prouposes
    db_annovar = db_annovar,#path annovar
    annovar_table_pl = annovar_table_pl, #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
    joinPY = joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py

   }

 #Array[File] salidas = ["${GatherBamFiles.output_bam}","${GatherBamFiles.output_bam_index}","${MergeVCFs.output_vcf}","${MergeVCFs.output_vcf_index}","${CollectGvcfCallingMetrics.summary_metrics}","${CollectGvcfCallingMetrics.detail_metrics}"]

 #scatter (paths in salidas) {
 #    call symlink_important_files {
 #        input:
 #        output_to_save = paths,
 #        path_save = path_save
 #    }
 #}    

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

 #Array[File] salidas_html = ConvertPairedFastQsToUnmappedBamWf.fastp_html_reports
 #Array[String] array_path_save_html = mkdir_samplename.path_out_softlink
 #Array[Pair[String,File]] samples_x_files_html = zip (array_path_save_html, salidas_html)
 #scatter (pairs in samples_x_files_html) {
 #    call symlink_important_files as save_fastp_html{
 #        input:
 #        output_to_save = pairs.right,
 #        path_save = pairs.left
 #    }
 #}
  
#Map[String,String] bams_N_reads = {"bam2gvcf.sampl_name_bam" : "${am2gvcf.N_reads_bam}"}




  call qual_control.quality_control {
   input: 
   stat_alineamiento = bam2gvcf.reporte_final,
   #bams_N_reads = bams_N_reads,
   fastp_json_files = ConvertPairedFastQsToUnmappedBamWf.fastp_json_reports,
   path_save = mkdir_samplename.path_out_softlink,
   analysis_readybam = bam2gvcf.analysis_ready_bam,
   toolpath = toolpath,
   Tso_name = basename(tabulatedSampleFilePaths, ".txt"),
   exon_coords = exon_coords
   #tso_bed = tso_bed
  }

   

 Array[File] prof_by_exon = quality_control.by_exon_depth
 Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink
 Array[Pair[String,File]] samples_by_exon = zip (array_path_save_byexon, prof_by_exon)
  scatter (pairs in samples_by_exon) {
    call symlink_important_files as byexon{
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
    }
  }

#String samplename2 
Array[File] Tsv_annovar = processJointVCF.annovar_tsv_out
    scatter (idx in range(length(Tsv_annovar))){
       call build_excell_report {
            input:
            annovar_tsv = Tsv_annovar[idx],
            exon_coverage_report = prof_by_exon[idx],
            sample = basename(Tsv_annovar[idx],".multianno_multisample.tsv"),
            samplename2 = basename(prof_by_exon[idx],"_coverage_statistics_by_exon.tsv")
           }
      
      

    }

Array[File] reporte_variantes = build_excell_report.excell_report
#Array[String] array_path_save_byexon = mkdir_samplename.path_out_softlink
 Array[Pair[String,File]] samples_by_variant = zip (array_path_save_byexon, reporte_variantes)
  scatter (pairs in samples_by_variant) {
    call symlink_important_files as byvariants{
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
    }
  }

 #Array[File] archivos_Apply = bam2gvcf.borrar_Applybqsr
 #scatter (archivos_app in archivos_Apply ){
 #call borrado as borrado_Apply{
 #input:
 #archivo_borrar = archivos_app
 #} 
 #}

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
   Array[File] output_vcf = bam2gvcf.output_vcf
   Array[File] output_vcf_index = bam2gvcf.output_vcf_index
  ##### output jointgenotype
   
   File? outputvcf = JointGenotyping.outputvcf
   File? outputvcfindex =  JointGenotyping.outputvcfindex
   File? detail_metrics_file =  JointGenotyping.metrica1
   File? summary_metrics_file = JointGenotyping.metrica2
   File? intervalo = JointGenotyping.inter

   #Array[File] Samt_bam_stat = bam2gvcf.Samt_bam_stat 
   Array[File] Samt_TSO_stat = bam2gvcf.Samt_TSO_stat
   Array[File] reporte_final = bam2gvcf.reporte_final ### archivo para mergear... estadistica en la libreria del experimento

  
  }



  
}

