task GenerateSubsettedContaminationResources {
  input {
    String bait_set_name
    File target_interval_list
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    Int preemptible_tries
  }

  String output_ud = bait_set_name + "." + basename(contamination_sites_ud)
  String output_bed = bait_set_name + "." + basename(contamination_sites_bed)
  String output_mu = bait_set_name + "." + basename(contamination_sites_mu)
  String target_overlap_counts = "target_overlap_counts.txt"

  command <<<
    set -e -o pipefail

    grep -vE "^@" ~{target_interval_list} |
       awk -v OFS='\t' '$2=$2-1' |
       /app/bedtools intersect -c -a ~{contamination_sites_bed} -b - |
       cut -f6 > ~{target_overlap_counts}

    function restrict_to_overlaps() {
        # print lines from whole-genome file from loci with non-zero overlap
        # with target intervals
        WGS_FILE=$1
        EXOME_FILE=$2
        paste ~{target_overlap_counts} $WGS_FILE |
            grep -Ev "^0" |
            cut -f 2- > $EXOME_FILE
        echo "Generated $EXOME_FILE"
    }

    restrict_to_overlaps ~{contamination_sites_ud} ~{output_ud}
    restrict_to_overlaps ~{contamination_sites_bed} ~{output_bed}
    restrict_to_overlaps ~{contamination_sites_mu} ~{output_mu}

  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "3.5 GiB"
    disks: "local-disk 10 HDD"
    docker: "us.gcr.io/broad-gotc-prod/bedtools:2.27.1"
  }
  output {
    File subsetted_contamination_ud = output_ud
    File subsetted_contamination_bed = output_bed
    File subsetted_contamination_mu = output_mu
  }
}

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
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Float contamination_underestimation_factor
    Boolean disable_sanity_check = false
  }

  Int disk_size = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB")) + 30

  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/~{contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "7.5 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:c1cba76e979904eb69c31520a0d7f5be63c72253-1553018888"
    cpu: 2
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}



# Check that the fingerprints of separate readgroups all match
task CrossCheckFingerprints {
  Array[File] input_bams
  Array[File] input_bam_indexes
  File? haplotype_database_file
  String metrics_filename
  String java_heap_memory_initial
  String memory
  Int cpu
  String tool_path
  
  command <<<
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx${java_heap_memory_initial} \
      -jar ${tool_path}/picard.jar \
      CrosscheckReadGroupFingerprints \
      OUTPUT=${metrics_filename} \
      HAPLOTYPE_MAP=${haplotype_database_file} \
      EXPECT_ALL_READ_GROUPS_TO_MATCH=true \
      INPUT=${sep=' INPUT=' input_bams} \
      LOD_THRESHOLD=-20.0 && \
      rm -rf ../inputs
  >>>
  runtime {
    memory: memory
    cpu: cpu
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

workflow contamin {


    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File target_interval_list
    File bait_interval_list
    File bait_set_name

 String cross_check_fingerprints_by = "READGROUP"


  call Processing.GenerateSubsettedContaminationResources {
    input:
        bait_set_name = bait_set_name,
        target_interval_list = target_interval_list,
        contamination_sites_bed = references.contamination_sites_bed,
        contamination_sites_mu = references.contamination_sites_mu,
        contamination_sites_ud = references.contamination_sites_ud,
        preemptible_tries = papi_settings.preemptible_tries
  }

  if (defined(haplotype_database_file) && defined(fingerprint_genotypes_file)) {
    # Check the sample BAM fingerprint against the sample array
    call CheckFingerprint {
      input:
        input_bam = GatherBamFiles.output_bam,
        input_bam_index = GatherBamFiles.output_bam_index,
        haplotype_database_file = haplotype_database_file,
        genotypes = fingerprint_genotypes_file,
        output_basename = base_file_name,
        sample = sample_name,
		tool_path = tool_path
    }
  }



# Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      input_bam = SortSampleBam.output_bam,
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      output_prefix = sample_and_unmapped_bams.base_file_name + ".preBqsr",
      preemptible_tries = papi_settings.agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }


}