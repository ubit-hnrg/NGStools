##### pipeline para obtener un GVCF a partir de un uBAM (enero 2019)
##### se unio fastqtosam y bwa
### se agrega genotypado
#### febrero 2019 se agrega array de entrada de ubams 


task GetBwaVersion {
  
  
  String path_herramientas

  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /tools/bwa-0.7.17/./bwa 2>&1 | \
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
    
      java -jar ${path_herramientas} ValidateSamFile -I=${input_bam} -M SUMMARY > >(tee ${sample_name}.Validatesamfile.stderr.log) | tail -n1 
    }

    output {    
        String valor = read_string(stdout())
        File stderr_log = "${sample_name}.Validatesamfile.stderr.log"
    }

}



# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File input_bam
  #String sample_name
  String path_herramientas

  Int compression_level
  String java_heap_memory_initial

  String output_bam_basename
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
    
    # set the bash variable needed for the command-line
    #bash_ref_fasta=${ref_fasta}

    java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${path_herramientas} \
        SamToFastq \
        --INPUT=${input_bam} \
        -F=/dev/stdout \
        --INTERLEAVE=true \
        --NON_PF=true | \
        /tools/bwa-0.7.17/./bwa mem -K 100000000 -p -v 3 -t 4 ${ref_fasta} /dev/stdin - 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
        samtools view -1 - > ${output_bam_basename}.bam
      
  >>>
 
  output {
    
    File output_bam = "${output_bam_basename}.bam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
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
  #String output_bam_basename
  

  command {
   java -Dsamjdk.compression_level=${compression_level} ${java_opt} -jar \
   ${path_herramientas} \
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
      --UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      --UNMAP_CONTAMINANT_READS=true
      
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
    java -Dsamjdk.compression_level=${compression_level} -Xmx${java_heap_memory_initial} -jar ${path_herramientas} \
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

    java -Dsamjdk.compression_level=${compression_level} ${java_opt_sort} -jar ${path_herramientas} \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    java -Dsamjdk.compression_level=${compression_level} ${java_opt_fix} -jar ${path_herramientas} \
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
    java -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -jar ${path_herramientas} \
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

  #String sample_name
  Array[File] input_bqsr_reports
  String output_report_filename

  String path_herramientas
  String java_opt

  command {
    java ${java_opt} -jar ${path_herramientas} \
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
    java ${java_opt} -jar ${path_herramientas} \
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
  #String sample_name
  
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xmx3g -jar ${path_herramientas} \
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
  String java_heap_memory_initial
  String path_herramientas
  #String sample_name
  String memory

  
  command <<<
    set -e
    mkdir out
    java -Dsamjdk.compression_level=${compression_level} -Xmx1g -jar ${path_herramientas} \
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

  #runtime {
   #   memory: memory
  #}
  output {
    Array[File] out = glob("out/*/*.interval_list")
    #Int interval_count = read_int(stdout())
    Int interval_count = read_int("error.txt")
    #File error = "error.txt"
  }

}

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
  Int compression_level
  String printreads_java_heap_memory_initial
  String haplotypecaller_java_heap_memory_initial
  String path_herramientas
  String sample_name
  String smith_waterman_implementation
  Float? contamination
  String newqual

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them. (no se usa, usa 100)
  #
  command <<<
      
      java -Xmx${haplotypecaller_java_heap_memory_initial} -jar ${path_herramientas} \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${gvcf_basename}.vcf.gz \
      -L ${interval_list} \
      -ip 100 \
      -contamination ${default=null contamination} \
      --max-alternate-alleles 3 \
      -ERC GVCF \
      --pair-hmm-implementation ${gatk_gkl_pairhmm_implementation} \
      --native-pair-hmm-threads ${gatk_gkl_pairhmm_threads} \
      --smith-waterman ${smith_waterman_implementation} \
      --use-new-qual-calculator ${newqual} \
      --sample-name ${sample_name}  
>>>

  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}



# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  #String sample_name
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  String path_herramientas

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar ${path_herramientas} \
      MergeVcfs \
      -I=${sep=' -I=' input_vcfs} \
      -O=${output_vcf_name}
  }

  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

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
  String path_herramientas


  command {
    java -Xms3000m -jar ${path_herramientas} \
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
  String sample_name
  File input_vcf
  File input_vcf_index
  String metrics_basename
  File dbSNP_vcf
  File dbSNP_vcf_index
  File ref_dict
  File wgs_evaluation_interval_list
  String path_herramientas


  command {
    java -Xms2000m -jar ${path_herramientas} \
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

#################################3 GENOTIPADO

task ImportGVCFs {
  String sample_names
  Array[File] input_gvcfs
  Array[File] input_gvcfs_indices
  File interval

  String workspace_dir_name

  String path_herramientas
  
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


    #padding (-ip) = 500, yo puse 100
    java -Xmx4g -Xms4g -jar ${path_herramientas} \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${interval} \
    --sample-name-map inputs.list \
    --reader-threads 5 \
    -ip 100

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

    >>>
 
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  File workspace_tar
  File interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index

  String path_herramientas

  command <<<
    set -e

    tar -xvf ${workspace_tar} 
    WORKSPACE=$(basename ${workspace_tar} .tar)

    java -Xmx4g -Xms4g -jar ${path_herramientas} \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbSNP_vcf} \
     -G StandardAnnotation \
     -new-qual true\
     -V gendb://$WORKSPACE \
     -L ${interval}
  >>>

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
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








##############################################    WORKFLOW
workflow pipeV0 {

    ### PATH local de las herramientas sacadas de docker
    String path_herramientas

    ####### MUESTRA
    String sample_name
    
    String ref_name

    File flowcell_unmapped_bams_list
    String unmapped_bam_suffix

    #File unmapped_bam
    
    String base_file_name = sample_name + "." + ref_name

    Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)
 
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


  ####genotipado
  File unpadded_intervals_file



    
    call GetBwaVersion {
        input: path_herramientas = path_herramientas
    }

  
    ## Align flowcell-level unmapped input bams in parallel
     scatter (unmapped_bam in flowcell_unmapped_bams) {

    ## Get the basename, i.e. strip the filepath and the extension
    String bam_basename = basename(unmapped_bam, unmapped_bam_suffix)

    #call validarbam {
    #    input:
    #     input_bam = unmapped_bam,
    #     path_herramientas = path_herramientas,
    #     sample_name = sample_name
    #}

  
######Si valor es = 0, el archivo ubam no tiene errores y prosigue. 
#### hay que meter un booleano en vez del string Valor1 para poder seguir el else y llamar a una funcion que muestre en pantalla "el archivo bam de entrada esta malformado"
##### la comparacion del if(valor1=0) debe ser booleana

    #String valor1 = validarbam.valor
    #if (valor1=="0") {
    # Map reads to reference
      call SamToFastqAndBwaMem {
      input:
        input_bam = unmapped_bam,
        path_herramientas = path_herramientas,
        compression_level = compression_level,
        java_heap_memory_initial = java_heap_memory_initial,
        output_bam_basename = bam_basename + ".unmerged",
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
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = bam_basename + ".aligned.unsorted",
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
  #}
 }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
      sample_name = sample_name,
      path_herramientas = path_herramientas, 
      input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
      # and the merged output.
      #disk_size = (md_disk_multiplier * SumFloats.total_size) + small_additional_disk,
      compression_level = compression_level,
	  
      java_heap_memory_initial = java_heap_memory_initial
  }


# Sort aggregated+deduped BAM file and fix tags
############### hay una version de wdl en la web que usa SamtoolsSort as SortSampleBam
  call SortAndFixTags {
    input:
      input_bam = MarkDuplicates.output_bam,
      sample_name = sample_name,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
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
        recalibration_report_filename = base_file_name + ".recal_data.csv",
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
      
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      path_herramientas = path_herramientas
      
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
        path_herramientas = path_herramientas,
        sample_name = sample_name
        
    }
} 

# Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      # Multiply the input bam size by two to account for the input and output
      #disk_size = (2 * agg_bam_size) + small_additional_disk,
      compression_level = compression_level,
	  path_herramientas = path_herramientas,
      #sample_name = sample_name,
      
      

}



############################ fin data preprocessing ##############################
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##################################################################################



#BQSR bins the qualities which makes a significantly smaller bam
#Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB")
call ScatterIntervalList {
  input:
      interval_list = wes_calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      compression_level = compression_level,
      #sample_name = sample_name,
      path_herramientas = path_herramientas,
     
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
        gvcf_basename = sample_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # Divide the total output GVCF size and the input bam size to account for the smaller scattered input and output.
        #disk_size = ((binned_qual_bam_size + GVCF_disk_size) / hc_divisor) + ref_size + small_additional_disk,
        compression_level = compression_level, 
        gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation, 
        gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
        path_herramientas = path_herramientas,
        sample_name = sample_name
        
		
     }
}

# Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs {
    input:
      #sample_name = sample_name,
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_name = sample_name + ".g.vcf.gz",
      path_herramientas = path_herramientas
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
      wgs_calling_interval_list = wgs_calling_interval_list,
      path_herramientas = path_herramientas
  }

# QC the GVCF
  call CollectGvcfCallingMetrics {
    input:
      sample_name = sample_name,
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = sample_name + ".g.vcf.gz",
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_dict = ref_dict,
      wgs_evaluation_interval_list = wgs_evaluation_interval_list,
      path_herramientas = path_herramientas
  }


  ###############################genotipado
#  Int num_of_original_intervals = length(read_lines(unpadded_intervals_file))
  ##Int num_gvcfs = length(MergeVCFs.output_vcf)
#  Int num_gvcfs = 1

  # Make a 2.5:1 interval number to samples in callset ratio interval list
#  Int possible_merge_count = floor(num_of_original_intervals / num_gvcfs / 2.5)
#  Int merge_count = if possible_merge_count > 1 then possible_merge_count else 1
  
#  call DynamicallyCombineIntervals {
#    input:
#      intervals = unpadded_intervals_file,
#      merge_count = merge_count
 # }

#Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)

#scatter (idx in range(length(unpadded_intervals))) {

#call ImportGVCFs {
#      input:
#        sample_names = sample_name,
#        interval = unpadded_intervals[idx],
#        workspace_dir_name = "mibasededatos",
#        input_gvcfs = MergeVCFs.output_vcf,
#        input_gvcfs_indices = MergeVCFs.output_vcf_index,
#        path_herramientas = path_herramientas,
#        batch_size = 50
#    }

#    call GenotypeGVCFs {
#      input:
#        workspace_tar = ImportGVCFs.output_genomicsdb,
#        interval = unpadded_intervals[idx],
#        output_vcf_filename = "output.vcf.gz",
#        ref_fasta = ref_fasta,
#        ref_fasta_index = ref_fasta_index,
#        ref_dict = ref_dict,
#        dbSNP_vcf = dbSNP_vcf,
#        dbSNP_vcf_index = dbSNP_vcf_index,
#        path_herramientas = path_herramientas
#}
#}










# Outputs that will be retained when execution is complete  
  output {
   File duplication_metrics = MarkDuplicates.duplicate_metrics
   File bqsr_report = GatherBqsrReports.output_bqsr_report
   File analysis_ready_bam = GatherBamFiles.output_bam
   File analysis_ready_bam_index = GatherBamFiles.output_bam_index
   File? analysis_ready_bam_md5 = GatherBamFiles.output_bam_md5
   File? gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
   File? gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics
   File output_vcf = MergeVCFs.output_vcf
   File output_vcf_index = MergeVCFs.output_vcf_index



} 






}
