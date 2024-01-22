######subworkflow del task fastqtobwa
###
###

workflow ubamtobwa {
    
    #Main input, array of ubams of (all samples) x (all lanes). i.e. containing all ubams of readgroups.  
    Array[File] array_unmapped_bams

    ### PATH local de las herramientas sacadas de docker
    String gatk_jar
    String toolpath
    ####### Reference name, .b37 , .hg19, etc. 
    String ref_name
    ########################command line para el bwa 
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
    Int compression_level
    String java_heap_memory_initial
  


    call GetBwaVersion {
      input:
      toolpath = toolpath
    }
    

    # Map reads to reference
    call Serial_SamToFastq_BwaMem_MergeBamAlignment {
      input:
      array_input_ubams = array_unmapped_bams,
      bwa_commandline = bwa_commandline,
      compression_level = compression_level,
      bwa_version = GetBwaVersion.version,
      java_heap_memory_initial = java_heap_memory_initial,
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

  output {
    Array[File] output_mergedbam_files = Serial_SamToFastq_BwaMem_MergeBamAlignment.output_mergedbam_files

  }

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






#All in one task
task Serial_SamToFastq_BwaMem_MergeBamAlignment {
  
  Array[File] array_input_ubams  
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

  for ubamfile in ${sep=' ' array_input_ubams}  ; do
        output_bwa_prefix=$(basename $ubamfile .unmapped.bam)  

   java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${toolpath}${gatk_jar} \
        SamToFastq \
        --INPUT $ubamfile \
        -F /dev/stdout \
        --INTERLEAVE true \
        --NON_PF true | \
        ${toolpath}${bwa_commandline} ${ref_fasta} /dev/stdin - 2> >(tee $output_filename.unmerged.bwa.stderr.log >&2) | \
        ${toolpath}samtools view -1 - > $output_bwa_prefix.aligned.unmerged.bam
  
    java -Dsamjdk.compression_level=${compression_level} -Xms3000m -jar \
    ${toolpath}${gatk_jar} \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ATTRIBUTES_TO_REMOVE NM \
        --ATTRIBUTES_TO_REMOVE MD \
        --ALIGNED_BAM $output_bwa_prefix.aligned.unmerged.bam \
        --UNMAPPED_BAM $ubamfile \
        -O $output_bwa_prefix.merged.unsorted.bam \
        -R ${ref_fasta} \
        --PAIRED_RUN true \
        --SORT_ORDER "unsorted" \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        -PG "bwamem2" \
        --PROGRAM_GROUP_NAME "bwamem2" \
        --PROGRAM_GROUP_VERSION "${bwa_version}" \
        --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline} ${ref_fasta}" \
        --UNMAP_CONTAMINANT_READS true \
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG\
        --ADD_PG_TAG_TO_READS false
  done        
  >>>
  
  
  #####se agrega de gatk2017
# ATTRIBUTES_TO_REMOVE=NM \
# ATTRIBUTES_TO_REMOVE=MD \
# PAIRED_RUN=true \
#
#
#
#######
  output {
    Array[File] output_mergedbam_files = glob("*merged.unsorted.bam")
  }
}

