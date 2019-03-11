
#### main para unir workflows


import "./hnrg-fastq2ubam_test.wdl" as fastq2ubam 
import "./ubam_to_gvcf_array.wdl" as preprocesing

workflow main_workflow {

###inputs for fastq2ubam workflows
  File tabulatedSampleFilePaths
  String run_date                   
  String library_name 
  String platform_unit 
  String platform_name 
  String gatk_jar
  String toolpath
  String ubam_list_name
  String sequencing_center =  "UIT-HNRG" 

#inputs for preprocesing ( ubam2gvcf)

    ####### MUESTRA
    #File sample_name
    
    String ref_name

    #File flowcell_unmapped_bams_list
    #String unmapped_bam_suffix

    #File unmapped_bam
    
    #String base_file_name 

    #Array[File] flowcell_unmapped_bams
    
 
    ########################command line para el bwa 
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
    File wgs_calling_interval_list
    File wgs_evaluation_interval_list
    File wes_calling_interval_list
    
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
 
  

  ####genotipado
    File unpadded_intervals_file
    
    #Array[File] inputs_ubams ###ubams del fastq2ubam


  
  call fastq2ubam.ConvertPairedFastQsToUnmappedBamWf {  
      input: 
      tabulatedSampleFilePaths = tabulatedSampleFilePaths,
      run_date = run_date,                   
      library_name = library_name,  
      platform_unit = platform_unit, 
      platform_name = platform_name,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      ubam_list_name = ubam_list_name,
      sequencing_center = sequencing_center

       }

Array[File] inputs_ubams = ConvertPairedFastQsToUnmappedBamWf.muestras
File uniquesample_name = ConvertPairedFastQsToUnmappedBamWf.samplesnames

scatter (sample in inputs_ubams)  {

    String sample_name = basename(sample, ".txt")

   call preprocesing.pipeV0 {
      input:
    #inputs_ubams = ConvertPairedFastQsToUnmappedBamWf.muestras,
    #uniquesample_name = ConvertPairedFastQsToUnmappedBamWf.samplesnames,
    
    sample_name = sample_name,
    flowcell_unmapped_bams_list = sample,
    ref_name = ref_name,
    unmapped_bam_suffix = ".bam",
    bwa_commandline = bwa_commandline,
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
    wes_calling_interval_list = wes_calling_interval_list,
    break_bands_at_multiples_of = break_bands_at_multiples_of,
    haplotype_scatter_count = haplotype_scatter_count,
    compression_level = compression_level,
    gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation,
    gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
    wgs_calling_interval_list = wgs_calling_interval_list,
    wgs_evaluation_interval_list = wgs_evaluation_interval_list,
    unpadded_intervals_file = unpadded_intervals_file,
    gatk_jar = gatk_jar,
    toolpath = toolpath,
    smith_waterman_implementation = smith_waterman_implementation,
    contamination = contamination,
    newqual = newqual,
    java_heap_memory_initial = java_heap_memory_initial
    } 

  }
    
    # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bams = ConvertPairedFastQsToUnmappedBamWf.output_bams
    Array[String] output_bams_sample_names =  ConvertPairedFastQsToUnmappedBamWf.output_bams_sample_names
    File unmapped_bam_list = ConvertPairedFastQsToUnmappedBamWf.unmapped_bam_list
   File samplesnames = ConvertPairedFastQsToUnmappedBamWf.samplesnames
    Array[File] muestras  =  ConvertPairedFastQsToUnmappedBamWf.muestras
  # Outputs del workflow pipeV0 (ubam2gvcf)  
 
   Array[File] duplication_metrics = pipeV0.duplication_metrics
   Array[File] bqsr_report = pipeV0.bqsr_report 
   Array[File] analysis_ready_bam = pipeV0.analysis_ready_bam
   Array[File] analysis_ready_bam_index = pipeV0.analysis_ready_bam_index
   Array[File] analysis_ready_bam_md5 = pipeV0.analysis_ready_bam_md5 
   Array[File] gvcf_summary_metrics = pipeV0.gvcf_summary_metrics 
   Array[File] gvcf_detail_metrics = pipeV0.gvcf_detail_metrics 
   Array[File]output_vcf = pipeV0.output_vcf
   Array[File] output_vcf_index = pipeV0.output_vcf_index
  }



  
}

