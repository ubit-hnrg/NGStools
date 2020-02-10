
## Requirements/expectations :
## - One or more human whole-genome per-sample GVCF files
##


workflow singleGenotypeGVCFs {
  File eval_interval_list
  
  String array_path_save
  String callset_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String gatk_jar
  String toolpath

  String sample_names
  File input_gvcfs 
  Array[File] input_gvcfs_indices 

  File dbSNP_vcf
  File dbSNP_vcf_index
  File region_padded_bed

  


  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69



    call GenotypeGVCFs {
      input:
        gvcf = input_gvcfs,
        gvcf_index = input_gvcfs_indices,
        output_vcf_filename = basename(input_gvcfs,"g.vcf.gz") + "single.output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        gatk_jar = gatk_jar,
        toolpath = toolpath
    
    }


    call HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GenotypeGVCFs.output_vcf,
        vcf_index = GenotypeGVCFs.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + ".single.variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + ".single.sites_only.variant_filtered.vcf.gz",
        gatk_jar = gatk_jar,
        toolpath = toolpath
   
    }
  


    call GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs_fofn  = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
        input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
        output_vcf_name = basename(input_gvcfs,"g.vcf.gz") + "single" + ".vcf.gz",
        gatk_jar = gatk_jar,
        toolpath = toolpath

    }

    call CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        gatk_jar = gatk_jar,
        toolpath = toolpath

    }
  

    call restrict_multisample_vcf{
        input:
        multisampleVCF  = FinalGatherVcf.output_vcf,
        region_padded_bed = region_padded_bed,
        toolpath = toolpath
    }

    call symlink_important_files {
       input:
         final_gath = FinalGatherVcf.output_vcf,
         final_gath_idx = FinalGatherVcf.output_vcf_index,
         metrica1 = CollectMetricsOnFullVcf.detail_metrics_file,
         metrica2 = CollectMetricsOnFullVcf.summary_metrics_file,
        path_save = array_path_save
   }
#}

  output {
    # outputs from the small callset path through the wdl
   File? output_singlevcf = FinalGatherVcf.output_vcf
   File? output_singlevcfindex =  FinalGatherVcf.output_vcf_index
   File? metrica1 =  CollectMetricsOnFullVcf.detail_metrics_file
   File? metrica2 = CollectMetricsOnFullVcf.summary_metrics_file
   File restricted_vcf = restrict_multisample_vcf.multisampleVCF_restricted
  }
}

task symlink_important_files {
    File final_gath
    File final_gath_idx
    File metrica1
    File metrica2
    String path_save
    command{
       cp -L ${final_gath} ${path_save}
       cp -L ${final_gath_idx} ${path_save}
       cp -L ${metrica1} ${path_save}
       cp -L ${metrica2} ${path_save}
    }
}

task GetNumberOfSamples {
  File sample_name_map

  command <<<
    wc -l ${sample_name_map} | awk '{print $1}'
  >>>

  output {
    Int sample_count = read_int(stdout())
  }
}


task GenotypeGVCFs {
  #String interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index    
 
  String gatk_jar
  String toolpath
  File gvcf
  Array[File] gvcf_index

  command <<<
    set -e


    java -Xmx1g -Xms1g -jar ${toolpath}${gatk_jar}\
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbSNP_vcf} \
     -G StandardAnnotation \
     --use-new-qual-calculator \
     --variant ${gvcf}
     
  >>>
#-L ${interval}
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename

  String gatk_jar
  String toolpath

  command {
    set -e

    java -Xmx3g -Xms3g -jar ${toolpath}${gatk_jar}\
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    java -Xmx3g -Xms3g -jar ${toolpath}${gatk_jar} \
      MakeSitesOnlyVcf \
     -I=${variant_filtered_vcf_filename} \
     -O=${sites_only_vcf_filename}

  }

  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  }
}



task GatherVcfs {
  File input_vcfs_fofn
  File input_vcf_indexes_fofn
  String output_vcf_name
  
  String gatk_jar
  String toolpath


    command <<<
    set -e
    set -o pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    java -Xmx6g -Xms6g -jar ${toolpath}${gatk_jar} \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    -I ${sep=" -I " input_vcfs_fofn} \
    --output ${output_vcf_name}

    java -Xmx6g -Xms6g -jar ${toolpath}${gatk_jar}\
    IndexFeatureFile \
   --feature-file ${output_vcf_name}
  >>>
#--input 
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {
  File input_vcf
  File input_vcf_index

  String metrics_filename_prefix
  File dbSNP_vcf
  File dbSNP_vcf_index
  File interval_list
  File ref_dict

  String gatk_jar
  String toolpath



  command {
    ##### --SEQUENCE_DICTIONARY ${ref_dict}
    java -Xmx4g -Xms4g -jar ${toolpath}${gatk_jar} \
      CollectVariantCallingMetrics \
      --INPUT ${input_vcf} \
      --DBSNP ${dbSNP_vcf} \
      --OUTPUT ${metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${interval_list}
  }
  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }

}

task GatherMetrics {
  File input_details_fofn
  File input_summaries_fofn

  String output_prefix

  String gatk_jar
  String toolpath
  

  command <<<
    set -e
    set -o pipefail

    
    java -Xmx2g -Xms2g -jar ${toolpath}${gatk_jar} \
    AccumulateVariantCallingMetrics \
    --INPUT ${sep=" --INPUT " input_details_fofn} \
    --OUTPUT ${output_prefix}
  >>>

  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}


task restrict_multisample_vcf{
    File multisampleVCF
    File region_padded_bed
    String toolpath
    String base = basename(multisampleVCF,'.vcf.gz')
    
    command{

        zless ${multisampleVCF} | java -jar ${toolpath}/SnpSift.jar intervals ${region_padded_bed} > ${base}_restricted.vcf
    }

    output {
        File multisampleVCF_restricted = "${base}_restricted.vcf"
    }

}