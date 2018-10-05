############################################################################################
###      Este pipeline NO esta optimizado (Octubre 2018)                                 ###
###   Implementa pre-procesamiento de datos de acuerdo a las best practices de GATK 2016 ###
###                                                                                      ###
###                                                                                      ###
### REQUERIMIENTOS:                                                                      ###
### datos de secuencias Pair-end en formato uBAM                                         ###
###                                                                                      ###
### SALIDA:                                                                              ###
###  archivo BAM limpio con su indice, adecuado para descubrimiento de variantes         ###
############################################################################################








task GetBwaVersion {
  
  
  String path_herramientas

  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    ${path_herramientas}/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
   
  output {
    String version = read_string(stdout())
  }
}

task validarbam {

String sample_name
File input_bam
String path_herramientas



command {
    
      java -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar ValidateSamFile -I=${input_bam} -M SUMMARY > >(tee ${sample_name}.Validatesamfile.stderr.log) | tail -n1 
 }

output {    
    String valor = read_string(stdout())
    File stderr_log = "${sample_name}.Validatesamfile.stderr.log"
}

}

task SamToFastq {
  File input_bam
  String sample_name
  
  
  String path_herramientas

  Int compression_level
  String java_heap_memory_initial
   
  command <<<
    set -o pipefail
    set -e

      java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
        SamToFastq \
        --INPUT=${input_bam} \
        -F=${sample_name}.f1.fastq \
        --INTERLEAVE=true \
        --NON_PF=true       
      
  >>>
 
  output {
    File FASTQ1 = "${sample_name}.f1.fastq"
  }
}

task BwaMem {
  
  String sample_name

  File fastaentrada 
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa


  String path_herramientas


    command <<<
    set -o pipefail
    set -e

    ${path_herramientas}/./bwa mem -K 100000000 -p -v 3 -t 4 ${ref_fasta} ${fastaentrada} > ${output_bam_basename}.bam
    >>>
    

    output {
    File output_bam = "${output_bam_basename}.bam"
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
  String path_herramientas
  String java_opt
  Int compression_level
  String sample_name
  

  command {
   java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar \
   ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
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
      --PROGRAM_RECORD_ID="bwamem" \
      --PROGRAM_GROUP_VERSION="${bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE="./bwa mem -K 100000000 -p -v 3 -t 4 ${ref_fasta}" \
      --PROGRAM_GROUP_NAME="bwamem" \
      --ALIGNER_PROPER_PAIR_FLAGS=true \
      
    }  
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {

  String sample_name  
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  Int compression_level
  String java_heap_memory_initial
  
  String path_herramientas
  



 # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files. 

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      MarkDuplicates \
      --INPUT=${sep=' INPUT=' input_bams} \
      --OUTPUT=${output_bam_basename}.bam \
      --METRICS_FILE=${metrics_filename} \
      --VALIDATION_STRINGENCY=SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      --ASSUME_SORT_ORDER="queryname" \
      --CLEAR_DT="false" \
      --ADD_PG_TAG_TO_READS=false \
    }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File input_bam
  String sample_name
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String path_herramientas
  
  Int compression_level

  String java_opt_sort
  String java_opt_fix

  command {
    set -o pipefail

    java -Dsamjdk.compression_level=${compression_level} ${java_opt_sort} -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    java -Dsamjdk.compression_level=${compression_level} ${java_opt_fix} -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ${output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ${ref_fasta}
  }

  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict  
  
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}


# Generate Base Quality Score Recalibration (BQSR) model

######dbSNP varia si es hg38 o b37
task BaseRecalibrator {
  String sample_name
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  
  String path_herramientas
  String java_opt

  command { 
    java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }

  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
task GatherBqsrReports {

  String sample_name
  Array[File] input_bqsr_reports
  String output_report_filename

  String path_herramientas
  String java_opt

  command {
    java ${java_opt} -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
 
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}


# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  String sample_name
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String path_herramientas
  String java_opt

  command {  
    java ${java_opt} -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      ApplyBQSR \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_bam_basename}.bam \
      -L ${sep=" -L " sequence_group_interval} \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }

  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename
  Int compression_level
  ##### son 3g, cambiar dps en el json
  #String java_heap_memory_initial
  String path_herramientas
  String sample_name
  
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xmx3g -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      GatherBamFiles \
      -I=${sep=' -I=' input_bams} \
      -O=${output_bam_basename}.bam \
      --CREATE_INDEX=true \
      --CREATE_MD5_FILE=true
    }
  
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }

}
########################################3 fin de pre processing#######################










##############################################    WORKFLOW
workflow pipeV0 {


    ### PATH local de las herramientas sacadas de docker
    String path_herramientas

    ####### MUESTRA
    String sample_name

    File unmapped_bam
    ########################command line para el bwa 
    String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 4 -Y"
  
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

    ##################################
    Int compression_level
    String java_heap_memory_initial
    
    call GetBwaVersion {
        input: path_herramientas = path_herramientas
    }

    call validarbam {
    input:
    input_bam = unmapped_bam,
    path_herramientas = path_herramientas,
    sample_name = sample_name
    }

######Si valor es = 0, el archivo ubam no tiene errores y prosigue. 
#### hay que meter un booleano en vez del string Valor1 para poder seguir el else y llamar a una funcion que muestre en pantalla "el archivo bam de entrada esta malformado"
##### la comparacion del if(valor1=0) debe ser booleana

String valor1 = validarbam.valor
if (valor1=="0") {
    call SamToFastq {
      input:
        input_bam = unmapped_bam,
        sample_name = sample_name,
        path_herramientas = path_herramientas,
        compression_level = compression_level,
        java_heap_memory_initial = java_heap_memory_initial
  }

  

  call BwaMem {
      input: 
        fastaentrada = SamToFastq.FASTQ1,
        bwa_commandline = bwa_commandline,
        sample_name = sample_name,
        output_bam_basename = sample_name + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        path_herramientas = path_herramientas        
  }

  # Merge original uBAM and BWA-aligned BAM 
    call MergeBamAlignment {
      input:
        unmapped_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        bwa_version = GetBwaVersion.version,
        aligned_bam = BwaMem.output_bam,
        sample_name = sample_name,
        output_bam_basename = sample_name + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        
        path_herramientas = path_herramientas,  
        compression_level = compression_level     
    }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
      sample_name = sample_name,
      path_herramientas = path_herramientas, 
      input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = sample_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_name + ".duplicate_metrics",
      # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
      # and the merged output.
      #disk_size = (md_disk_multiplier * SumFloats.total_size) + small_additional_disk,
      compression_level = compression_level,
	  
      java_heap_memory_initial = java_heap_memory_initial
  }


# Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags {
    input:
      input_bam = MarkDuplicates.output_bam,
      sample_name = sample_name,
      output_bam_basename = sample_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      compression_level = compression_level,
      path_herramientas = path_herramientas
      
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
    
  }


  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        sample_name = sample_name,
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        recalibration_report_filename = sample_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        path_herramientas = path_herramientas
        
    }  
}
###fin scatter1

# Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      sample_name = sample_name,
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = sample_name + ".recal_data.csv",
      path_herramientas = path_herramientas
      
}

scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = sample_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        path_herramientas = path_herramientas,
        sample_name = sample_name
        
    }
} 

# Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = sample_name,
      # Multiply the input bam size by two to account for the input and output
      #disk_size = (2 * agg_bam_size) + small_additional_disk,
      compression_level = compression_level,
	    path_herramientas = path_herramientas,
      sample_name = sample_name
}


# Outputs that will be retained when execution is complete  
  output {
   File? duplication_metrics = MarkDuplicates.duplicate_metrics
   File? bqsr_report = GatherBqsrReports.output_bqsr_report
   File? analysis_ready_bam = GatherBamFiles.output_bam
   File? analysis_ready_bam_index = GatherBamFiles.output_bam_index
   File? analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
} 

}

############################ fin data preprocessing ##############################






########fin if



}