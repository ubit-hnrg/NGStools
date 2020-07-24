#####task del summary_metrics de samtools
task samtools_stat{

###herramientas
#String gatk_jar

String toolpath
File TSO_bed #./TruSight_One_v1_padded_100_GRCh37.bed
File input_bam_reducido
String name

command <<<
    set -e
    set -o pipefail
 #bgzip ${input_bam_reducido} | tabix -p sam /dev/stdout | bgzip -d | 
#${toolpath}samtools sort ${input_bam_reducido} | 
##input es el bam recortado  y el intervalo de captura
${toolpath}samtools stats ${input_bam_reducido}  -t ${TSO_bed} > ${name}_TSO_samtools.stats

>>>
output {

#File samtools_stats = 
#File samtools_reduced_bam = $name'_samtools_reduced.stats'
#File samtools_stat_original_bam = "${name}_orig_samtools.stats"
File samtools_stat_TSO_bam = "${name}_TSO_samtools.stats"

#File 

}

}


task samtools_reports_file {

String sampleID
Int N_total_reads_bam
#File samtools_global_report ##no va mas, necesita el numero total de reads
File samtools_library_report
String toolpath

#String path_salida -T=${N_total_reads_bam}

command {
/home/hnrg/NGStools/pipeline_wdl/qualityControl/samtools_stats_report_v1.0.py -N=${N_total_reads_bam}  -l=${samtools_library_report} -o=${sampleID}_samtools_report.tsv

}

output {
 
File output_global_report = "${sampleID}_samtools_report.tsv" 

}

}




task bams_reads {
   File bam
   String toolpath
   String sample_name
   #  echo ${sample_name}\t$({toolpath}samtools view -c ${bam})

   command {
  
    echo |${toolpath}samtools view -c ${bam} 
   }

   output {
   String sampl = sample_name
   Int N_reads = read_int(stdout())
   }

  }


  task reduce_bam {
   File input_bam
   String toolpath
   String output_bam_basename
   File lib_resctricted 
 


   command <<<
   ${toolpath}bedtools2/bin/intersectBed -a ${input_bam} -b ${lib_resctricted} -wa > ${output_bam_basename}_lib_resctricted.bam 

   >>>
    output {
    File output_reduced_bam = "${output_bam_basename}_lib_resctricted.bam"
    
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
    }
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


task symlink_important_files {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}


##############################################    WORKFLOW
workflow bam2gvcf {

   String path_save
 
    ### PATH local de las herramientas sacadas de docker
    String gatk_jar
    String toolpath


    ####### MUESTRA
   # String ref_name
   
    #String unmapped_bam_suffix




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
    File wes_calling_interval_list
    Int break_bands_at_multiples_of
    Int haplotype_scatter_count

    ##################################
    Int compression_level
    String java_heap_memory_initial

    ########optimization flags
    String gatk_gkl_pairhmm_implementation
    Int gatk_gkl_pairhmm_threads


  #validate gvcf
  
  

  #archivo de listas de cobertura de intervalos 
  #File wgs_coverage_interval_list
  File wgs_calling_interval_list

  File wgs_evaluation_interval_list


  #####opt de haplotypecaller
  String smith_waterman_implementation
  Float? contamination
  String newqual
  File tso_bed
  


  #String sample_name
  
  #File flowcell_unmapped_bams_list ### archivo txt con ubams de 1 sample (muchos ubam)
  
  
  String base_file_name
  #Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)
  File lib_resctricted
 #./TruSight_One_v1_padded_100_GRCh37.bed 

 
 
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
      # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
      # and the merged output.
      #disk_size = (md_disk_multiplier * SumFloats.total_size) + small_additional_disk,
      compression_level = compression_level,
      java_heap_memory_initial = java_heap_memory_initial,
      gatk_jar = gatk_jar,
        toolpath = toolpath
  }




 #scatter (bam in bams_N_reads){
   call bams_reads {
     input:
    bam = MarkDuplicates.output_bam,
     toolpath = toolpath,
     sample_name = base_file_name
    
    }
    #}

     call reduce_bam {
       input:
       input_bam = MarkDuplicates.output_bam, 
       #input_bam = bam_markdup,
       toolpath = toolpath,
       output_bam_basename = base_file_name, 
       lib_resctricted = lib_resctricted
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


  #call borrado as borrar_Markdup {

   #input:
   #  archivo_borrar = MarkDuplicates.output_bam
 #}


    


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
#
Array[File] borrado_apply = ["${ApplyBQSR.recalibrated_bam}"]


} 

  #call path_borrado as borrar_SortandFix {
#
#   input:
#     path1 = SortAndFixTags.output_bam
# }


# Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      # Multiply the input bam size by two to account for the input and output
      #disk_size = (2 * agg_bam_size) + small_additional_disk,
      compression_level = compression_level,
      gatk_jar = gatk_jar,
      toolpath = toolpath     
      
}

############################ fin data preprocessing ##############################
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##################################################################################
#call samtools_stat {
#      input:
#      toolpath = toolpath,
#      name = base_file_name, 
#      TSO_bed = tso_bed, #./TruSight_One_v1_padded_100_GRCh37.bed
#      input_bam_reducido = GatherBamFiles.output_bam,
      #input_bam_reducido = reduce_bam.output_reduced_bam

#     }
   
#  call samtools_reports_file {

#    input: 
#    sampleID = base_file_name,
#    N_total_reads_bam = bams_reads.N_reads,
    #samtools_global_report = samtools_stat.samtools_stat_original_bam,
#    samtools_library_report = samtools_stat.samtools_stat_TSO_bam,
#    toolpath = toolpath

#  }




#BQSR bins the qualities which makes a significantly smaller bam
#Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB")
#call ScatterIntervalList {
#  input:
#      interval_list = wes_calling_interval_list,
#      scatter_count = haplotype_scatter_count,
#      break_bands_at_multiples_of = break_bands_at_multiples_of,
#      compression_level = compression_level,
#      gatk_jar = gatk_jar,
#        toolpath = toolpath     
#}

#Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
#Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  #Call variants in parallel over WGS calling intervals
#  scatter (index in range(ScatterIntervalList.interval_count)) {
    ###Generate GVCF by interval
#    call HaplotypeCaller {
#      input:
#        input_bam = GatherBamFiles.output_bam,
#        input_bam_index = GatherBamFiles.output_bam_index,
#        interval_list = ScatterIntervalList.out[index],
#        gvcf_basename = base_file_name,
#        ref_dict = ref_dict,
#        ref_fasta = ref_fasta,
#        ref_fasta_index = ref_fasta_index,
        # Divide the total output GVCF size and the input bam size to account for the smaller scattered input and output.
        #disk_size = ((binned_qual_bam_size + GVCF_disk_size) / hc_divisor) + ref_size + small_additional_disk,
#        compression_level = compression_level, 
#        gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation, 
#        gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
#        gatk_jar = gatk_jar,
#        toolpath = toolpath,
#        smith_waterman_implementation = smith_waterman_implementation,
#        contamination = contamination,
#        newqual = newqual
         
        
		
#     }
#}

#  call GatherBamFilesHaplotype {
#    input:
#      input_bams = HaplotypeCaller.bam_haplotypecaller,
#      output_bam_basename = base_file_name,
      # Multiply the input bam size by two to account for the input and output
      #disk_size = (2 * agg_bam_size) + small_additional_disk,
#      compression_level = compression_level,
#      gatk_jar = gatk_jar,
#      toolpath = toolpath     
      
#}

# Combine by-interval GVCFs into a single sample GVCF file
#  call MergeVCFs {
#    input:
      
  
#      input_vcfs = HaplotypeCaller.output_gvcf,
#      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
#      output_vcf_name = base_file_name + ".g.vcf.gz",
#      gatk_jar = gatk_jar,
#      toolpath = toolpath
#      }

# Validate the GVCF output of HaplotypeCaller
#  call ValidateGVCF {
#    input:
#      input_vcf = MergeVCFs.output_vcf,
#      input_vcf_index = MergeVCFs.output_vcf_index,
#      dbSNP_vcf = dbSNP_vcf,
#      dbSNP_vcf_index = dbSNP_vcf_index,
#      ref_fasta = ref_fasta,
#      ref_fasta_index = ref_fasta_index,
#      ref_dict = ref_dict,
#      wgs_calling_interval_list = wgs_calling_interval_list,
#      gatk_jar = gatk_jar,
#        toolpath = toolpath
#    }

# QC the GVCF
#  call CollectGvcfCallingMetrics {
#    input:
      
#      input_vcf = MergeVCFs.output_vcf,
#      input_vcf_index = MergeVCFs.output_vcf_index,
#      metrics_basename = base_file_name + ".g.vcf.gz",
#      dbSNP_vcf = dbSNP_vcf,
#      dbSNP_vcf_index = dbSNP_vcf_index,
#      ref_dict = ref_dict,
#      wgs_evaluation_interval_list = wgs_evaluation_interval_list,
#      gatk_jar = gatk_jar,
#        toolpath = toolpath
#  }
  


#Array[File] salidas = ["${GatherBamFilesHaplotype.output_bam}","${GatherBamFilesHaplotype.output_bam_index}","${GatherBamFiles.output_bam}","${GatherBamFiles.output_bam_index}","${MergeVCFs.output_vcf}","${MergeVCFs.output_vcf_index}","${CollectGvcfCallingMetrics.summary_metrics}","${CollectGvcfCallingMetrics.detail_metrics}","${samtools_stat.samtools_stat_TSO_bam}","${samtools_reports_file.output_global_report}"]

Array[File] salidas = ["${GatherBamFiles.output_bam}","${GatherBamFiles.output_bam_index}"]
scatter (paths in salidas) {
    call symlink_important_files {
        input:
        output_to_save = paths,
        path_save = path_save
    }
}


# Outputs that will be retained when execution is complete  
  output {



  #####outputs workflow ubam2gvcf

#   File duplication_metrics = MarkDuplicates.duplicate_metrics
#   File bqsr_report = GatherBqsrReports.output_bqsr_report
   File analysis_ready_bam = GatherBamFiles.output_bam
   File analysis_ready_bam_index = GatherBamFiles.output_bam_index
   File analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
#   File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
#   File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics
#   File output_vcf = MergeVCFs.output_vcf
#   File output_vcf_index = MergeVCFs.output_vcf_index

#   Array[File] borrar_Applybqsr = ApplyBQSR.recalibrated_bam 
   #File borrar_Markdup = MarkDuplicates.output_bam
#   File borrar_SortandFix = SortAndFixTags.output_bam
   #String sampl_name_bam = bams_reads.sampl 
   #String N_reads_bam = bams_reads.N_reads 
   #File Samt_bam_stat = samtools_stat.samtools_stat_original_bam 
#   File Samt_TSO_stat = samtools_stat.samtools_stat_TSO_bam
#   File reporte_final = samtools_reports_file.output_global_report ### archivo para mergear... estadistica en la libreria del experimento


   #"samtools_stat.samtools_stat_TSO_bam","samtools_reports_file.output_global_report"

} 


}
