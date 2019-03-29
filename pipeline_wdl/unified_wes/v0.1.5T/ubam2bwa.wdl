######subworkflow del task fastqtobwa
###
###

workflow ubamtobwa {
    
    #String path_save
 
    ### PATH local de las herramientas sacadas de docker
    String gatk_jar
    String toolpath


    ####### Reference name, .b37 , .hg19, etc. 
    String ref_name
    #String base_file_name
   
    
    
 
    ########################command line para el bwa 
    #String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 4 -Y"
    String bwa_commandline



    ########## referencia
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
   

    ##################################
    Int compression_level
    String java_heap_memory_initial

  


  String sample_name_ubam
  
  File flowcell_unmapped_bams_list ### archivo txt con ubams de 1 sample (muchos ubam)
  
  
  String base_file_name = sample_name_ubam  + ref_name
  
  Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)
  #String unmapped_bam_suffix
  #String base_file_name
  #File unmapped_bam


    call GetBwaVersion {
      input:
      toolpath = toolpath
    }


    
    

    # Map reads to reference
    call SamToFastqAndBwaMem2 {
      input:
        array_input_bam = array_unmapped_bams,
        bwa_commandline = bwa_commandline,
        compression_level = compression_level,
        java_heap_memory_initial = java_heap_memory_initial,
        #output_bam_basename = bam_basename + ".unmerged",
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


  ## Align flowcell-level unmapped input bams in parallel
  #scatter (unmapped_bam in flowcell_unmapped_bams) {

  ## Get the basename, i.e. strip the filepath and the extension
  #String bam_basename = basename(unmapped_bam, '.bam')



  # Merge original uBAM and BWA-aligned BAM 
  #  call MergeBamAlignment {
  #    input:
  #      array_unmapped_bam = unmapped_bam,
  #      #output_bam_basename = bam_basename + ".aligned.unsorted",

#        bwa_commandline = bwa_commandline,
#        bwa_version = GetBwaVersion.version,
#        aligned_bam = SamToFastqAndBwaMem.output_bam,
#        ref_fasta = ref_fasta,
#        ref_fasta_index = ref_fasta_index,
#        ref_dict = ref_dict,
#        ref_amb = ref_amb,
#        ref_ann = ref_ann,
#        ref_bwt = ref_bwt,
#        ref_pac = ref_pac,
#        ref_sa = ref_sa,
#        compression_level = compression_level,
#        gatk_jar = gatk_jar,
#        toolpath = toolpath     
#    }
#output {
#  Array[File] bam_alineado = MergeBamAlignment.output_bam
#  String nombre_base = base_file_name  

#}
#}

}


task GetBwaVersion {
String toolpath


#/tools/bwa-0.7.17/./bwa 2>&1 | \
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
   ${toolpath}./bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
   
  output {
    String version = read_string(stdout())
  }
}




# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  Array[File] array_input_bams
  
  String gatk_jar
  String toolpath
  String bwa_commandline

  Int compression_level
  String java_heap_memory_initial

  #String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa  
        

  command <<<
    set -o pipefail
    set -e

  for ubamfile in ${sep=' ' array_input_bams}  ; do
    output_filename=$(basename $ubamfile)


   java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${toolpath}${gatk_jar} \
        SamToFastq \
        --INPUT=$ubamfile \
        -F=/dev/stdout \
        --INTERLEAVE=true \
        --NON_PF=true | \
        ${toolpath}${bwa_commandline} ${ref_fasta} /dev/stdin - 2> >(tee $output_filename.unmerged.bwa.stderr.log >&2) | \
        ${toolpath}samtools view -1 - > $output_filename.unmerged.bam
  
  done        
      
  >>>
  output {
    Array[File] output_bam_files = glob(*.bam)
    #File output_bam = "${output_bam_basename}.bam"
    Array[File] bwa_stderr_log = glob(*log)
  }
}


#All in one
task SamToFastqAndBwaMem2 {
  Array[File] array_input_bams
  
  String gatk_jar
  String toolpath
  String bwa_commandline
  String bwa_version
  Int compression_level
  String java_heap_memory_initial
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa  
        

  command <<<
    set -o pipefail
    set -e

  for ubamfile in ${sep=' ' array_input_bams}  ; do
    output_bwa_filename=$(basename $ubamfile)
    #output_filename=$(basename $ubamfile)

   java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${toolpath}${gatk_jar} \
        SamToFastq \
        --INPUT=$ubamfile \
        -F=/dev/stdout \
        --INTERLEAVE=true \
        --NON_PF=true | \
        ${toolpath}${bwa_commandline} ${ref_fasta} /dev/stdin - 2> >(tee $output_filename.unmerged.bwa.stderr.log >&2) | \
        ${toolpath}samtools view -1 - > $output_bwa_filename.unmerged.bam
  
    java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar \
    ${toolpath}${gatk_jar} \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY=SILENT \
        --EXPECTED_ORIENTATIONS=FR \
        --ATTRIBUTES_TO_RETAIN=X0 \
        --ALIGNED_BAM=$output_bwa_filename.unmerged.bam \
        --UNMAPPED_BAM=$ubamfile \
        -O=$output_bwa_filename.unmerged.aligned.unsorted.bam \
        -R=${ref_fasta} \
        --SORT_ORDER="unsorted" \
        --IS_BISULFITE_SEQUENCE=false \
        --ALIGNED_READS_ONLY=false \
        --CLIP_ADAPTERS=false \
        --MAX_RECORDS_IN_RAM=2000000 \
        --ADD_MATE_CIGAR=true \
        --MAX_INSERTIONS_OR_DELETIONS=-1 \
        --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS=true \
        -PG="bwamem" \
        --PROGRAM_GROUP_NAME="bwamem" \
        --PROGRAM_GROUP_VERSION="${bwa_version}" \
        --PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline} ${ref_fasta}" \
        --UNMAP_CONTAMINANT_READS=true \
        --UNMAPPED_READ_STRATEGY=COPY_TO_TAG\
        --ADD_PG_TAG_TO_READS=false
  done        
  >>>
  
  output {
    Array[File] output_bwa_files = glob(*unmerged.bam)
    Array[File] output_mergedbam_files = glob(*aligned.unsorted.bam)
    #File output_bam = "${output_bam_basename}.bam"
    Array[File] bwa_stderr_log = glob(*log)
  }
}



task MergeBamAlignment {
  File unmapped_bam
  String bwa_commandline
  String bwa_version

  File aligned_bam
  String output_bam_basename

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  String gatk_jar
  String toolpath
  Int compression_level

  command <<<
   java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar \
   ${toolpath}${gatk_jar} \
       MergeBamAlignment \
      --VALIDATION_STRINGENCY=SILENT \
      --EXPECTED_ORIENTATIONS=FR \
      --ATTRIBUTES_TO_RETAIN=X0 \
      --ALIGNED_BAM=${aligned_bam} \
      --UNMAPPED_BAM=${unmapped_bam} \
      -O=${output_bam_basename}.bam \
      -R=${ref_fasta} \
      --SORT_ORDER="unsorted" \
      --IS_BISULFITE_SEQUENCE=false \
      --ALIGNED_READS_ONLY=false \
      --CLIP_ADAPTERS=false \
      --MAX_RECORDS_IN_RAM=2000000 \
      --ADD_MATE_CIGAR=true \
      --MAX_INSERTIONS_OR_DELETIONS=-1 \
      --PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      --ALIGNER_PROPER_PAIR_FLAGS=true \
      -PG="bwamem" \
      --PROGRAM_GROUP_NAME="bwamem" \
      --PROGRAM_GROUP_VERSION="${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline} ${ref_fasta}" \
      --UNMAP_CONTAMINANT_READS=true \
      --UNMAPPED_READ_STRATEGY=COPY_TO_TAG\
      --ADD_PG_TAG_TO_READS=false
    >>>
 
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}
