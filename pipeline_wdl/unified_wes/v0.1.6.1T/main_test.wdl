#### main para unir workflows


import './hnrg-fastq2ubam_test.wdl' as fastq2ubam 
import './bam2gvcf.wdl' as bamtogvcf
import './ubam2bwa.wdl' as ubam2bwa
import './jointgenotype_sinrecalibracion.wdl' as joint_genotype

task borrado {
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

task borrado_single {
Array[File] path1 
#File path2

#command <<<
#mv ${write_lines(path)}  ${t}.txt;
command <<<
while read -r line; do
rm "$line"
done < "${path1}"
>>>

}

task mkdir_samplename {
    String path_softlink
    String samplename

    command{
        mkdir  ${path_softlink}${samplename}
    }

    output {
        String path_out_softlink = "${path_softlink}" + "${samplename}"+"/"
}
}


task subset_array_glob {
  String pattern 
  Array[File] array_of_files

  command{
  }

  output {
    Array[File] subArray_input_ubam2gvcf = glob('../inputs/*/*${pattern}*')
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


  call fastq2ubam.ConvertPairedFastQsToUnmappedBamWf {  
      input: 
      tabulatedSampleFilePaths = tabulatedSampleFilePaths,
      run_date = run_date,                   
      library_name = library_name,   
      platform_name = platform_name,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      ubam_list_name = ubam_list_name,
      sequencing_center = sequencing_center,
      read_lenght = readlenght,
      platform_model = platformmodel

      }

  call borrado {
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

   call subset_array_glob {
     input: 
       pattern = sample_name,
       array_of_files = ubamtobwa.output_mergedbam_files

   }

   call bamtogvcf.bam2gvcf {
      input:
    #base_file_name =  ubamtobwa.nombre_base,
    base_file_name =  sample_name,                                          #### ARTEFACTO PARA PROBAR EL WORKFLOW ANTERIOR, ESTO HAY QUE TRABAJRLO Y ARTICULARLO BIEN. 
    #inputs_ubams = ConvertPairedFastQsToUnmappedBamWf.muestras,
    #uniquesample_name = ConvertPairedFastQsToUnmappedBamWf.samplesnames,
    path_save = mkdir_samplename.path_out_softlink,
    #sample_name = sample_name,
    #flowcell_unmapped_bams_list = sample,
    bams_entrada = subset_array_glob.subArray_input_ubam2gvcf,
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
    java_heap_memory_initial = java_heap_memory_initial
    } 

  }

#Array[File] paths2symb = ["${pepe.dummy_outfile}","${juan.dummy_outfile}"]
#scatter (paths in paths2symb) {
#    call symlink_a_lo_guapo {
#        input:
#        output_to_save = paths 
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
  
  }



  
}

