####bam2gvcf workflow

#####task del summary_metrics de samtools
task samtools_stat{

  String toolpath
  File intervalo_captura #./TruSight_One_v1_padded_100_GRCh37.bed
  File input_bam_reducido
  String name
  String path_save


  command <<<
    set -e
    set -o pipefail
 
    ${toolpath}samtools stats ${input_bam_reducido}  -t ${intervalo_captura} > ${name}_experiment_lib_samtools.stats
    /usr/local/bin/plot-bamstats ${name}_experiment_lib_samtools.stats -p ${path_save}samtools_plots/${name}
  

  >>>
  output {

  File samtools_stat_experiment_bam = "${name}_experiment_lib_samtools.stats" 

  }

}




task reduce_bam {
   File input_bam
   String toolpath
   String output_bam_basename
   File lib_restricted 
 


   command <<<

   echo |${toolpath}samtools view -c ${input_bam} > bams_reads_markdup.txt

   ${toolpath}bedtools2/bin/intersectBed -a ${input_bam} -b ${lib_restricted} -wa > ${output_bam_basename}_lib_restricted.bam 

   >>>
    output {
    File output_reduced_bam = "${output_bam_basename}_lib_restricted.bam"
    String N_reads = read_string("bams_reads_markdup.txt")
    
   }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {

  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  Int compression_level
  String java_heap_memory_initial
  
  String gatk_jar
  String toolpath
  
   # Aggregate aligned+merged flowcell BAM files and mark duplicates
   # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
   # to avoid having to spend time just merging BAM files. 

   # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
   # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
   # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
   #-Xmx${java_heap_memory_initial}
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar ${toolpath}${gatk_jar} \
      MarkDuplicates \
      --INPUT=${sep=' --INPUT=' input_bams} \
      --OUTPUT=${output_bam_basename}.bam \
      --METRICS_FILE=${metrics_filename} \
      --VALIDATION_STRINGENCY=SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      --ASSUME_SORT_ORDER="queryname" \
      --CREATE_MD5_FILE=true \
      --CLEAR_DT="false" \
      --ADD_PG_TAG_TO_READS=false
    }

    ####se agrega
    ##       CLEAR_DT="false" \
    ##  ADD_PG_TAG_TO_READS=false
    ##
    ##
    #######
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
   File input_bam
   String output_bam_basename
   File ref_dict
   File ref_fasta
   File ref_fasta_index
   String gatk_jar
   String toolpath
  
   Int compression_level


    command {
     set -o pipefail

     java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar ${toolpath}${gatk_jar} \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
      | \
     java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar ${toolpath}${gatk_jar} \
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
  
  String gatk_jar
  String toolpath

  command { 
    java -Xms4000m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar ${toolpath}${gatk_jar} \
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

  Array[File] input_bqsr_reports
  String output_report_filename

  String gatk_jar
  String toolpath


  command {
    java -Xms3000m -jar ${toolpath}${gatk_jar} \
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
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_jar
  String toolpath
 

  command {  
    java -Xms3000m -jar ${toolpath}${gatk_jar} \
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
  String gatk_jar
  String toolpath
  
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xmx3g -jar ${toolpath}${gatk_jar} \
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

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  File interval_list

  #Scatter_count = 50
  Int scatter_count
  ## break_bands_at_multiples_of= 1000000
  Int break_bands_at_multiples_of
  Int compression_level 
  
  String gatk_jar
  String toolpath
  

  
  command <<<
    set -e
    mkdir out
    java -Dsamjdk.compression_level=${compression_level} -Xmx1g -jar ${toolpath}${gatk_jar} \
      IntervalListTools \
      -I=${interval_list} \
      -O=out \
      --SCATTER_COUNT=${scatter_count} \
      --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      --UNIQUE=true \
      --SORT=true \
      --BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} 
      

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals),file=open("error.txt", "a"))  
    CODE
  >>>

  
  output {
    Array[File] out = glob("out/*/*.interval_list")
    #Int interval_count = read_int(stdout())
    Int interval_count = read_int("error.txt")
    #File error = "error.txt"
  }

}

#########new
# PRIVATE #
# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
#   contamination estimate for use in variant calling

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  File input_bam
  File input_bam_index
  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File ref_fasta
  File ref_fasta_index
  String output_prefix
  Float contamination_underestimation_factor
  String toolpath

  command <<<
    set -e

    # creates a ${output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    ${toolpath}VerifyBamID/bin/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ${output_prefix} \
    --BamFile ${input_bam} \
    --Reference ${ref_fasta} \
    --UDPath ${contamination_sites_ud} \
    --MeanPath ${contamination_sites_mu} \
    --BedPath ${contamination_sites_bed} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/${contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
 
  output {
    File selfSM = "${output_prefix}.selfSM"
    File depthSM = "${output_prefix}.depthSM"
    File log = "${output_prefix}.log"
    Float contamination = read_float(stdout())
  }
}


########




# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String gatk_gkl_pairhmm_implementation
  Int gatk_gkl_pairhmm_threads
  #Int compression_level
  String gatk_jar
  String toolpath


  String smith_waterman_implementation
  Float? contamination
  String newqual

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them. (no se usa, usa 100)
  #
  command <<<
      
      java -Xmx2g -jar ${toolpath}${gatk_jar} \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${gvcf_basename}.g.vcf.gz \
      -L ${interval_list} \
      -ip 100 \
      --contamination ${default=0 contamination} \
      --max-alternate-alleles 3 \
      -ERC GVCF \
      --pair-hmm-implementation ${gatk_gkl_pairhmm_implementation} \
      --native-pair-hmm-threads ${gatk_gkl_pairhmm_threads} \
      --smith-waterman ${smith_waterman_implementation} \
      --use-new-qual-calculator ${newqual} \
      --bam-output= ${gvcf_basename}_haplotype.bam ###realigned reads \
      ## -new-qual \
      #-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \

       
  >>>

  output {
    File output_gvcf = "${gvcf_basename}.g.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.g.vcf.gz.tbi"
    File bam_haplotypecaller = "${gvcf_basename}_haplotype.bam"

  }
}


task GatherBamFilesHaplotype {
  Array[File] input_bams
  String output_bam_basename
  Int compression_level
  
  ##### son 3g, cambiar dps en el json
  #String java_heap_memory_initial
String gatk_jar
  String toolpath
  
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xmx3g -jar ${toolpath}${gatk_jar} \
      GatherBamFiles \
      -I=${sep=' -I=' input_bams} \
      -O=${output_bam_basename}_haplotype.bam \
      --CREATE_INDEX=true \
      --CREATE_MD5_FILE=true
    }
  
  output {
    File output_bam = "${output_bam_basename}_haplotype.bam"
    File output_bam_index = "${output_bam_basename}_haplotype.bai"
    File output_bam_md5 = "${output_bam_basename}_haplotype.bam.md5"
  }

}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  String gatk_jar
  String toolpath


  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar ${toolpath}${gatk_jar} \
      MergeVcfs \
      -I=${sep=' -I=' input_vcfs} \
      -O=${output_vcf_name} 
  }

  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

############################# GATK 2018, ver  si funciona. 

task HardfilterVCF {
   
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File interval_list
    Int preemptible_tries
    String gatk_jar
    String toolpath
  

  String output_vcf_name = vcf_basename + ".filtered.vcf.gz"
  command {
    java -Xms3000m -jar ${toolpath}${gatk_jar} \
      VariantFiltration \
      -V ${input_vcf} \
      -L ${interval_list} \
      --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
      --filter-name "HardFiltered" \
      -O ${output_vcf_name}
  }
  output {
      File output_vcf = "${output_vcf_name}"
      File output_vcf_index = "${output_vcf_name}.tbi"
      }

}

#####CNNscorevariants es experimental. ojo.
task CNNScoreVariants {

  
    File? bamout
    File? bamout_index
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    #Int preemptible_tries
    String gatk_jar
    String toolpath
  

  String base_vcf = basename(input_vcf)
  #Boolean is_compressed = basename(base_vcf, "gz") != base_vcf
  String vcf_suffix = ".vcf.gz"
  String vcf_index_suffix =  ".tbi" 
  String output_vcf = base_vcf + ".scored" + vcf_suffix
  String output_vcf_index = output_vcf + vcf_index_suffix

  #String bamout_param = if defined(bamout) then "-I ~{bamout}" else ""
  #String tensor_type = if defined(bamout) then "read-tensor" else "reference" ver diferencias
#${tensor_type}
  command {
     java -Xmx10g -jar ${toolpath}${gatk_jar} \
      CNNScoreVariants \
       -V ${input_vcf} \
       -R ${ref_fasta} \
       -O ${output_vcf} \
       -I ${bamout} \
       -tensor-type "read-tensor" 
  }

  output {
  
    File scored_vcf = "${output_vcf}"
    File scored_vcf_index = "${output_vcf_index}"
  }
}

# task FilterVariantTranches {

  
#    File input_vcf
#    File input_vcf_index
#    String vcf_basename
#    Array[String] snp_tranches
#    Array[String] indel_tranches
#    File hapmap_resource_vcf
#    File hapmap_resource_vcf_index
#    File omni_resource_vcf
#    File omni_resource_vcf_index
#    File one_thousand_genomes_resource_vcf
#    File one_thousand_genomes_resource_vcf_index
#    File dbsnp_resource_vcf
#    File dbsnp_resource_vcf_index
#    String info_key
#    Int preemptible_tries
  


#  command {

#    gatk --java-options -Xmx6g FilterVariantTranches \
#      -V ${input_vcf} \
#      -O ${vcf_basename}.filtered.vcf.gz \
#      ${sep=" " prefix("--snp-tranche ", snp_tranches)} \
#      ${sep=" " prefix("--indel-tranche ", indel_tranches)} \
#      --resource ${hapmap_resource_vcf} \
#      --resource ${omni_resource_vcf} \
#      --resource ${one_thousand_genomes_resource_vcf} \
#      --resource ${dbsnp_resource_vcf} \
#      --info-key ${info_key} \
#      --create-output-variant-index true
#  }

#  output {
#    File filtered_vcf = "${vcf_basename}.filtered.vcf.gz"
#    File filtered_vcf_index = "${vcf_basename}.filtered.vcf.gz.tbi"
#  }
#}

########################### fin gatk nuevo 

# Validate a GVCF with -gvcf specific validation
task ValidateGVCF {
  File input_vcf
  File input_vcf_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbSNP_vcf
  File dbSNP_vcf_index
  File wgs_calling_interval_list
  String gatk_jar
  String toolpath



  command {
    java -Xms3000m -jar ${toolpath}${gatk_jar} \
      ValidateVariants \
      -V ${input_vcf} \
      --reference ${ref_fasta} \
      -L ${wgs_calling_interval_list} \
      -gvcf \
      --validation-type-to-exclude null \
      --dbsnp ${dbSNP_vcf}
  }


}

# Collect variant calling metrics from GVCF output
task CollectGvcfCallingMetrics {
  File input_vcf
  File input_vcf_index
  String metrics_basename
  File dbSNP_vcf
  File dbSNP_vcf_index
  File ref_dict
  File wgs_evaluation_interval_list
  String gatk_jar
  String toolpath


  command {
    java -Xms2000m -jar ${toolpath}${gatk_jar} \
      CollectVariantCallingMetrics \
      -I=${input_vcf} \
      -O=${metrics_basename} \
      --DBSNP=${dbSNP_vcf} \
      --SEQUENCE_DICTIONARY=${ref_dict} \
      --TARGET_INTERVALS=${wgs_evaluation_interval_list} \
      --GVCF_INPUT=true
  }
 
  output {
    File summary_metrics = "${metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "${metrics_basename}.variant_calling_detail_metrics"
  }
}


#############task nuevos
#
# Check that the fingerprints of separate readgroups all match

task CrossCheckFingerprints {
  Array[File] input_bams
  Array[File] input_bam_indexes
  File? haplotype_database_file
  String metrics_filename
  String java_heap_memory_initial
  String toolpath
  String gatk_jar
  
  command <<<
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx${java_heap_memory_initial} \
      -jar ${toolpath}${gatk_jar} \
      CrosscheckReadGroupFingerprints \
      OUTPUT=${metrics_filename} \
      HAPLOTYPE_MAP=${haplotype_database_file} \
      EXPECT_ALL_READ_GROUPS_TO_MATCH=true \
      INPUT=${sep=' INPUT=' input_bams} \
      LOD_THRESHOLD=-20.0
  >>>
  output {
    File metrics = "${metrics_filename}" 
     }
}

# Check that the fingerprint of the sample BAM matches the sample array
task CheckFingerprint {
  File input_bam
  File input_bam_index
  String output_basename
  File? haplotype_database_file
  File? genotypes
  String sample
  String java_heap_memory_initial
  String gatk_jar
  String toolpath
  
  command <<<
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx${java_heap_memory_initial}  \
      -jar ${toolpath}${gatk_jar} \
      CheckFingerprint \
      INPUT=${input_bam} \
      OUTPUT=${output_basename} \
      GENOTYPES=${genotypes} \
      HAPLOTYPE_MAP=${haplotype_database_file} \
      SAMPLE_ALIAS="${sample}" \
      IGNORE_READ_GROUPS=true
  >>>
 
  output {
    File summary_metrics = "${output_basename}.fingerprinting_summary_metrics"
    File detail_metrics = "${output_basename}.fingerprinting_detail_metrics" 
   }
}


task symlink_important_files {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}


#Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  File input_bam
  String output_bam_prefix = basename(input_bam, ".merged.unsorted.bam")  
  String gatk_jar
  String toolpath


  command {
    java -Xms5000m -jar ${toolpath}${gatk_jar} \
      CollectMultipleMetrics \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectBaseDistributionByCycle" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="MeanQualityByCycle" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="ALL_READS"

    touch ${output_bam_prefix}.insert_size_metrics
    touch ${output_bam_prefix}.insert_size_histogram.pdf
  }
 
  output {
    File base_distribution_by_cycle_pdf = "${output_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "${output_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "${output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "${output_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "${output_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "${output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${output_bam_prefix}.quality_distribution_metrics"
  }
}

##############################################    WORKFLOW
workflow bam2gvcf {

  String path_save
 
  ### PATH local de las herramientas sacadas de docker
  String gatk_jar
  String toolpath

  ########## referencia
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? ref_amb
  File? ref_ann
  File? ref_bwt
  File? ref_pac
  File? ref_sa
    
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  File wes_calling_interval_list
  Int break_bands_at_multiples_of
  Int haplotype_scatter_count

  ##################################
  Int compression_level
  String java_heap_memory_initial

  ########optimization flags
  String gatk_gkl_pairhmm_implementation
  Int gatk_gkl_pairhmm_threads

  File wgs_calling_interval_list

  File wgs_evaluation_interval_list

  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu

  #####opt de haplotypecaller
  String smith_waterman_implementation
  Float? contamination
  String newqual
  File intervalo_captura
    
  
  String base_file_name
  #Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)
  File lib_restricted
  #./TruSight_One_v1_padded_100_GRCh37.bed 

 #########agrego 2023 hg38
 #esto es con verify bam id. 
  #File? contamination_sites_vcf
  #File? contamination_sites_vcf_index
  #File? fingerprint_genotypes_file # if this file is empty (0-length) the workflow should not do fingerprint comparison (as there are no fingerprints for the sample)
  #File haplotype_database_file  
 
  Array[File] bams_entrada

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  
  #File bam_markdup


  
  call MarkDuplicates {
    input:
      
      input_bams = bams_entrada,
      #input_bams = reduced_bams,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      compression_level = compression_level,
      java_heap_memory_initial = java_heap_memory_initial,
      gatk_jar = gatk_jar,
        toolpath = toolpath
  }




##########esto lo tengo que checkear
  call reduce_bam {
    input:
    input_bam = MarkDuplicates.output_bam, 
    #input_bam = bam_markdup,
    toolpath = toolpath,
    output_bam_basename = base_file_name, 
    lib_restricted = lib_restricted #intervalo_captura
    #./TruSight_One_v1_padded_100_GRCh37.bed 
  }
  



  # Sort aggregated+deduped BAM file and fix tags
  ############### hay una version de wdl en la web que usa SamtoolsSort as SortSampleBam
  call SortAndFixTags {
    input:
    #input_bam = MarkDuplicates.output_bam,
    input_bam = reduce_bam.output_reduced_bam,
    output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
    ref_dict = ref_dict,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    compression_level = compression_level,
    gatk_jar = gatk_jar,
    toolpath = toolpath      
  }

########esto es nuevo
# if (defined(haplotype_database_file)) {
#   # Check identity of fingerprints across readgroups
#    call CrossCheckFingerprints {
#      input:
#        input_bams = SortAndFixTags.output_bam,
#        input_bam_indexes = SortAndFixTags.output_bam_index,
#        haplotype_database_file = haplotype_database_file,
#        metrics_filename = sample_name + ".crosscheck",
#		    tool_path = tool_path
#    }
#  }
#
######### terminar de definir checkfingerprint


  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
    ref_dict = ref_dict
    
  }

  # Estimate level of cross-sample contamination
  #call CheckContamination {
  #  input:
  #    input_bam = SortAndFixTags.output_bam,
  #    input_bam_index = SortAndFixTags.output_bam_index,
  #    contamination_sites_ud = contamination_sites_ud,
  #    contamination_sites_bed = contamination_sites_bed,
  #    contamination_sites_mu = contamination_sites_mu,
  #    ref_fasta = ref_fasta,
  #    ref_fasta_index = ref_fasta_index,
  #    output_prefix = base_file_name + ".preBqsr",
  #    contamination_underestimation_factor = 0.75,
  #    toolpath = toolpath
  #}



  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
      input_bam = SortAndFixTags.output_bam,
      input_bam_index = SortAndFixTags.output_bam_index,
      recalibration_report_filename = base_file_name + ".recal_data.csv",
      sequence_group_interval = subgroup,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gatk_jar = gatk_jar,
      toolpath = toolpath        
    }  
  }
  ###fin scatter1

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      
    input_bqsr_reports = BaseRecalibrator.recalibration_report,
    output_report_filename = base_file_name + ".recal_data.csv",
    gatk_jar = gatk_jar,
    toolpath = toolpath      
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        gatk_jar = gatk_jar,
        toolpath = toolpath    
        
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
    input_bams = ApplyBQSR.recalibrated_bam,
    output_bam_basename = base_file_name,
    compression_level = compression_level,
    gatk_jar = gatk_jar,
    toolpath = toolpath     
      
  }



##########agregamos 2023
#
 #if (defined(haplotype_database_file) && defined(fingerprint_genotypes_file)) {
    # Check the sample BAM fingerprint against the sample array
    #call CheckFingerprint {
    #  input:
    #    input_bam = GatherBamFiles.output_bam,
    #    input_bam_index = GatherBamFiles.output_bam_index,
    #    haplotype_database_file = haplotype_database_file,
    #    genotypes = fingerprint_genotypes_file,
    #    output_basename = base_file_name,
    #    sample = base_file_name,
		#    tool_path = tool_path
    #}
  #}
####


######
  ############################ fin data preprocessing ##############################  
  ## Output :
  ## - A clean BAM file and its index, suitable for variant discovery analyses.
  ##################################################################################
  # call samtools_stat {
  #   input:
  #   toolpath = toolpath,
  #   name = base_file_name, 
  #   TSO_bed = tso_bed, #./TruSight_One_v1_padded_100_GRCh37.bed
  #   input_bam_reducido = GatherBamFiles.output_bam,
  #   #input_bam_reducido = reduce_bam.output_reduced_bam

  # }
   
  # call samtools_reports_file {

  # input: 
  # sampleID = base_file_name,
  # N_total_reads_bam = bams_reads.N_reads,
  # #samtools_global_report = samtools_stat.samtools_stat_original_bam,
  # samtools_library_report = samtools_stat.samtools_stat_TSO_bam,
  # toolpath = toolpath

  # }

  #BQSR bins the qualities which makes a significantly smaller bam
  #Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB")
  call ScatterIntervalList {
    input:
    interval_list = wes_calling_interval_list,
    scatter_count = haplotype_scatter_count,
    break_bands_at_multiples_of = break_bands_at_multiples_of,
    compression_level = compression_level,
    gatk_jar = gatk_jar,
    toolpath = toolpath     
  }

  #Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  #Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  #Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    ###Generate GVCF by interval
    call HaplotypeCaller {
      input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      interval_list = ScatterIntervalList.out[index],
      gvcf_basename = base_file_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation, 
      gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      smith_waterman_implementation = smith_waterman_implementation,
      contamination =contamination, # CheckContamination.contamination, #contamination,
      ##agrego esto
      #contamination = PreBqsrCheckContamination.contamination,
      
      newqual = newqual
     }
  }

  call GatherBamFilesHaplotype {
    input:
      input_bams = HaplotypeCaller.bam_haplotypecaller,
      output_bam_basename = base_file_name,
      compression_level = compression_level,
      gatk_jar = gatk_jar,
      toolpath = toolpath     
      
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs {
  input:
  input_vcfs = HaplotypeCaller.output_gvcf,
  input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
  output_vcf_name = base_file_name + ".g.vcf.gz",
  gatk_jar = gatk_jar,
  toolpath = toolpath
  }

  # Validate the GVCF output of HaplotypeCaller
  call ValidateGVCF {
  input:
  input_vcf = MergeVCFs.output_vcf,
  input_vcf_index = MergeVCFs.output_vcf_index,
  dbSNP_vcf = dbSNP_vcf,
  dbSNP_vcf_index = dbSNP_vcf_index,
  ref_fasta = ref_fasta,
  ref_fasta_index = ref_fasta_index,
  ref_dict = ref_dict,
  wgs_calling_interval_list = wes_calling_interval_list, #wgs_calling_interval_list,
  gatk_jar = gatk_jar,
  toolpath = toolpath
  }

  # QC the GVCF
  call CollectGvcfCallingMetrics {
  input:
      
  input_vcf = MergeVCFs.output_vcf,
  input_vcf_index = MergeVCFs.output_vcf_index,
  metrics_basename = base_file_name, #+ ".g.vcf.gz",
  dbSNP_vcf = dbSNP_vcf,
  dbSNP_vcf_index = dbSNP_vcf_index,
  ref_dict = ref_dict,
  wgs_evaluation_interval_list = wgs_evaluation_interval_list,
  gatk_jar = gatk_jar,
  toolpath = toolpath
  }
  
  #####

  




  Array[File] salidas = ["${GatherBamFilesHaplotype.output_bam}","${GatherBamFilesHaplotype.output_bam_index}","${GatherBamFiles.output_bam}","${GatherBamFiles.output_bam_index}","${MergeVCFs.output_vcf}","${MergeVCFs.output_vcf_index}","${CollectGvcfCallingMetrics.summary_metrics}","${CollectGvcfCallingMetrics.detail_metrics}"]# ,"${samtools_stat.samtools_stat_TSO_bam}","${samtools_reports_file.output_global_report}"]

  scatter (paths in salidas) {
    call symlink_important_files {
    input:
    output_to_save = paths,
    path_save = path_save
    }
  }


 #   Outputs that will be retained when execution is complete  
  output {
   #####outputs workflow ubam2gvcf

   File duplication_metrics = MarkDuplicates.duplicate_metrics
   File bqsr_report = GatherBqsrReports.output_bqsr_report
   File analysis_ready_bam = GatherBamFiles.output_bam
   File analysis_ready_bam_index = GatherBamFiles.output_bam_index
   File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
   File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
   File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics
   File output_gvcf = MergeVCFs.output_vcf
   File output_gvcf_index = MergeVCFs.output_vcf_index

   #Array[File] borrar_Applybqsr = ApplyBQSR.recalibrated_bam 
   #File borrar_Markdup = MarkDuplicates.output_bam
   #File borrar_SortandFix = SortAndFixTags.output_bam
   #Int N_reads_bam = reduce_bam.N_reads
   #String sampl_name_bam = bams_reads.sampl 
   #String N_reads_bam = bams_reads.N_reads 
   #File Samt_bam_stat = samtools_stat.samtools_stat_original_bam 
   #File Samt_TSO_stat = samtools_stat.samtools_stat_TSO_bam
   #File reporte_final = samtools_reports_file.output_global_report ### archivo para mergear... estadistica en la libreria del experimento


   #File? cross_check_fingerprints_metrics = CrossCheckFingerprints.metrics
   #File? fingerprint_summary_metrics = CheckFingerprint.summary_metrics
   #File? fingerprint_detail_metrics = CheckFingerprint.detail_metrics




  } 


}
