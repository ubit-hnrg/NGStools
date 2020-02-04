
## Requirements/expectations :
## - One or more human whole-genome per-sample GVCF files
##


workflow JointGenotyping {
  File unpadded_intervals_file ##lista intervalos
  File eval_interval_list
  
  Array[String] array_path_save
  String callset_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String gatk_jar
  String toolpath

  Array[String] sample_names
  Array[File] input_gvcfs 
  Array[File] input_gvcfs_indices 

  File dbSNP_vcf
  File dbSNP_vcf_index
  
  


  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold


  Int num_of_original_intervals = length(read_lines(unpadded_intervals_file))
  Int num_gvcfs = length(input_gvcfs)

  # Make a 2.5:1 interval number to samples in callset ratio interval list
  Int possible_merge_count = floor(num_of_original_intervals / num_gvcfs / 2.5)
  Int merge_count = if possible_merge_count > 1 then possible_merge_count else 1

  
  call DynamicallyCombineIntervals {
    input:
      intervals = unpadded_intervals_file,
      merge_count = merge_count
  }

  Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)

  scatter (idx in range(length(unpadded_intervals))) {
    # the batch_size value was carefully chosen here as it
    # is the optimal value for the amount of memory allocated
    # within the task; please do not change it without consulting
    # the Hellbender (GATK engine) team!
    call ImportGVCFs {
      input:
        sample_names = sample_names,
        interval = unpadded_intervals[idx],
        workspace_dir_name = "genomicsdb",
        batch_size = 50,
        input_gvcfs = input_gvcfs,
        input_gvcfs_indices = input_gvcfs_indices,
        gatk_jar = gatk_jar,
        toolpath = toolpath
    
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
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
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        gatk_jar = gatk_jar,
        toolpath = toolpath
   
    }
  }


  # For small callsets (fewer than 1000 samples) we can gather the VCF shards and collect metrics directly.
  # For anything larger, we need to keep the VCF sharded and gather metrics collected from them.
  Boolean is_small_callset = num_gvcfs <= 1000

  # for small callsets we can gather the VCF shards and then collect metrics on it
  if (is_small_callset) {
    call GatherVcfs as FinalGatherVcf {
      input:
        ##input_vcfs_fofn = GenotypeGVCFs.output_vcf,
        input_vcfs_fofn  = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
        #input_vcfs_fofn = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
        ##input_vcf_indexes_fofn = GenotypeGVCFs.output_vcf_index,
        input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
        #input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index,
        ##input_vcfs_fofn = SitesOnlyGatherVcf.output_vcf,
        ##input_vcf_indexes_fofn = SitesOnlyGatherVcf.output_vcf_index,
        output_vcf_name = callset_name + ".vcf.gz",
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
  }

Array[File] salidas = ["${FinalGatherVcf.output_vcf}","${FinalGatherVcf.output_vcf_index}","${CollectMetricsOnFullVcf.detail_metrics_file}","${CollectMetricsOnFullVcf.summary_metrics_file}"]
Array[Pair[String,File]] samples_x_files = cross (array_path_save, salidas)
scatter (pairs in samples_x_files) {
    call symlink_important_files {
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
    }
}

  output {
    # outputs from the small callset path through the wdl
   File? outputvcf = FinalGatherVcf.output_vcf
   File? outputvcfindex =  FinalGatherVcf.output_vcf_index
   File? metrica1 =  CollectMetricsOnFullVcf.detail_metrics_file
   File? metrica2 = CollectMetricsOnFullVcf.summary_metrics_file

    # outputs from the large callset path through the wdl
    # (note that we do not list ApplyRecalibration here because it is run in both paths)
    #GatherMetrics.detail_metrics_file
    #GatherMetrics.summary_metrics_file

    # output the interval list generated/used by this run workflow
    File? inter = DynamicallyCombineIntervals.output_intervals
  }
}


task symlink_important_files {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
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

task ImportGVCFs {
  Array[String] sample_names
  Array[File] input_gvcfs
  Array[File] input_gvcfs_indices


  String interval
  String gatk_jar
  String toolpath

  String workspace_dir_name


  Int batch_size

  command <<<
    set -e
    set -o pipefail
    
    python << CODE
    gvcfs = ['${sep="','" input_gvcfs}']
    sample_names = ['${sep="','" sample_names}']

    if len(gvcfs)!= len(sample_names):
      exit(1)

    with open("inputs.list", "w") as fi:
      for i in range(len(gvcfs)):
       fi.write(sample_names[i] + "\t" + gvcfs[i] + "\n") 
    
    CODE
    
    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    #-ip 500
    java -Xmx1g -Xms1g -jar ${toolpath}${gatk_jar} \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${interval} \
    --sample-name-map inputs.list  \
    --reader-threads 4 \
    

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}
  
  >>>

  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index    
 
  String gatk_jar
  String toolpath

  command <<<
    set -e

    tar -xf ${workspace_tar}
    WORKSPACE=$( basename ${workspace_tar} .tar)

    java -Xmx1g -Xms1g -jar ${toolpath}${gatk_jar}\
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbSNP_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals \
     --use-new-qual-calculator \
     -V gendb://$WORKSPACE \
     -L ${interval}
  >>>

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
  Array[File] input_vcfs_fofn
  Array[File] input_vcf_indexes_fofn
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
  --input ${output_vcf_name}
  >>>
 
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
  Array[File] input_details_fofn
  Array[File] input_summaries_fofn

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

task DynamicallyCombineIntervals {
  File intervals
  Int merge_count


  command {
    python << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ${merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("${intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
  }


  output {
    File output_intervals = "out.intervals"
  }
}
