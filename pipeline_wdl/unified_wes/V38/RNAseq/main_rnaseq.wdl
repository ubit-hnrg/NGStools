## Copyright Broad Institute, 2018
##
## Workflows for processing RNA data for germline short variant discovery with GATK (v3+v4) and related tools 
##
## Requirements/expectations :
## - BAM 
##
## Output :
## - A BAM file and its index.
## - A VCF file and its index. 
## - A Filtered VCF file and its index. 
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs. 
 
 workflow RNAseq {

	File inputBam
	String sampleName = basename(inputBam,".bam")

	File refFasta
	File refFastaIndex
	File refDict

	
	File wgsCallingIntervalList

	Array[File] knownVcfs
	Array[File] knownVcfsIndices

	File dbSnpVcf
	File dbSnpVcfIndex

	Int? minConfidenceForVariantCalling

	## Inputs for STAR
	Int? readLength
	File? zippedStarReferences
	File annotationsGTF
  
	## Optional user optimizations
	Int? haplotypeScatterCount
	Int scatterCount = select_first([haplotypeScatterCount, 6])
	Boolean? use_gatk4_for_all_tools
	Boolean use_all_gatk4 = select_first([use_gatk4_for_all_tools, false])


	call RevertSam {
		input:
			input_bam = inputBam,
			base_name = sampleName + ".reverted",
			sort_order = "queryname",
			gatk_path = gatk_path
	}

	call SamToFastq {
		input:
			unmapped_bam = RevertSam.output_bam,
			base_name = sampleName,
			gatk_path = gatk_path
	}

	if (!defined(zippedStarReferences)) {

		call StarGenerateReferences { 
			input:
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				annotations_gtf = annotationsGTF,
				read_length = readLength
		}
	}

	File starReferences = select_first([zippedStarReferences,StarGenerateReferences.star_genome_refs_zipped,""])

	call StarAlign { 
		input: 
			star_genome_refs_zipped = starReferences,
			fastq1 = SamToFastq.fastq1,
			fastq2 = SamToFastq.fastq2,
			base_name = sampleName + ".star",
			read_length = readLength,
	}

	call MergeBamAlignment {
		input: 
			unaligned_bam = RevertSam.output_bam,
			star_bam = StarAlign.output_bam,
			base_name = ".merged",
			ref_fasta = refFasta,
			ref_dict = refDict,
			gatk_path = gatk_path
	}

	call MarkDuplicates {
		input:
			input_bam = MergeBamAlignment.output_bam,
			base_name = sampleName + ".dedupped"
	}
	if (!use_all_gatk4){
		call SplitNCigarReads {
			input:
				input_bam = MarkDuplicates.output_bam,
				input_bam_index = MarkDuplicates.output_bam_index,
				base_name = sampleName + ".split",
				ref_fasta = refFasta,
				ref_fasta_index = refFastaIndex,
				ref_dict = refDict,
				interval_list = wgsCallingIntervalList
		}
	}

	if (use_all_gatk4){
               call SplitNCigarReads_GATK4 {
                        input:
                                input_bam = MarkDuplicates.output_bam,
                                input_bam_index = MarkDuplicates.output_bam_index,
                                base_name = sampleName + ".split",
                                ref_fasta = refFasta,
                                ref_fasta_index = refFastaIndex,
                                ref_dict = refDict,
                                interval_list = wgsCallingIntervalList,
                                gatk_path = gatk_path
                }
	}

	call BaseRecalibrator {
		input:
			input_bam = select_first([SplitNCigarReads_GATK4.output_bam, SplitNCigarReads.output_bam]),
			input_bam_index = select_first([SplitNCigarReads_GATK4.output_bam_index, SplitNCigarReads.output_bam_index]),
			recal_output_file = sampleName + ".recal_data.csv",
  			dbSNP_vcf = dbSnpVcf,
  			dbSNP_vcf_index = dbSnpVcfIndex,
  			known_indels_sites_VCFs = knownVcfs,
  			known_indels_sites_indices = knownVcfsIndices,
  			ref_dict = refDict,
  			ref_fasta = refFasta,
  			ref_fasta_index = refFastaIndex,
			gatk_path = gatk_path
	}

	call ApplyBQSR {
		input:
			input_bam =  select_first([SplitNCigarReads_GATK4.output_bam, SplitNCigarReads.output_bam]),
			input_bam_index = select_first([SplitNCigarReads_GATK4.output_bam_index, SplitNCigarReads.output_bam_index]),
			base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			gatk_path = gatk_path
	}
        if (!use_all_gatk4){
		call ScatterIntervalList {
			input:
				interval_list = wgsCallingIntervalList,
				scatter_count = scatterCount
		}
	}

        if (use_all_gatk4){
                call ScatterIntervalList_GATK4 {
                        input:
				interval_list = wgsCallingIntervalList,
				scatter_count = scatterCount,
				gatk_path = gatk_path
                }
        }

	scatter (i in range(select_first([ScatterIntervalList_GATK4.interval_count, ScatterIntervalList.interval_count]))) {
		if (!use_all_gatk4){
			call HaplotypeCaller {
				input:
					input_bam = ApplyBQSR.output_bam,
					input_bam_index = ApplyBQSR.output_bam_index,
					base_name = sampleName + ".hc",
					interval_list = select_first([ScatterIntervalList_GATK4.out,ScatterIntervalList.out])[i],
					ref_fasta = refFasta,
					ref_fasta_index = refFastaIndex,
					ref_dict = refDict,
					dbSNP_vcf = dbSnpVcf,
					dbSNP_vcf_index = dbSnpVcfIndex,
					stand_call_conf = minConfidenceForVariantCalling,
			}
		}
		if (use_all_gatk4){
			call HaplotypeCaller_GATK4 {
        			input:
	                        	input_bam = ApplyBQSR.output_bam,
	                                input_bam_index = ApplyBQSR.output_bam_index,
        	                        base_name = sampleName + ".hc",
                	                interval_list = select_first([ScatterIntervalList_GATK4.out,ScatterIntervalList.out])[i],
                        	        ref_fasta = refFasta,
                                	ref_fasta_index = refFastaIndex,
	                                ref_dict = refDict,
        	                        dbSNP_vcf = dbSnpVcf,
                	                dbSNP_vcf_index = dbSnpVcfIndex,
                        	        stand_call_conf = minConfidenceForVariantCalling,
					gatk_path = gatk_path
			}
		}
		File HaplotypeCallerOutputVcf = select_first([HaplotypeCaller_GATK4.output_vcf, HaplotypeCaller.output_vcf])
		File HaplotypeCallerOutputVcfIndex = select_first([HaplotypeCaller_GATK4.output_vcf_index, HaplotypeCaller.output_vcf_index])
	}

        call MergeVCFs {
                input:
                        input_vcfs = HaplotypeCallerOutputVcf,
                        input_vcfs_indexes =  HaplotypeCallerOutputVcfIndex,
                        output_vcf_name = sampleName + ".g.vcf.gz",
                        preemptible_count = preemptible_count,
                        docker = gatk4_docker,
                        gatk_path = gatk_path
        }
	
	call VariantFiltration {
		input:
			input_vcf = MergeVCFs.output_vcf,
			input_vcf_index = MergeVCFs.output_vcf_index,
			base_name = sampleName + ".variant_filtered.vcf.gz",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			preemptible_count = preemptible_count,
			docker = gatk4_docker,
			gatk_path = gatk_path
	}

	output {
		File recalibrated_bam = ApplyBQSR.output_bam
		File recalibrated_bam_index = ApplyBQSR.output_bam_index
		File merged_vcf = MergeVCFs.output_vcf
		File merged_vcf_index = MergeVCFs.output_vcf_index
		File variant_filtered_vcf = VariantFiltration.output_vcf
		File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
	}
}


#NOTE: assuming aggregated bams & paired end fastqs
task SamToFastq {
	File unmapped_bam
	String base_name

  String gatk_path

  

	command <<<
	 	${gatk_path} \
	 	    SamToFastq \
	 	    --INPUT ${unmapped_bam} \
	 	    --VALIDATION_STRINGENCY SILENT \
	 	    --FASTQ ${base_name}.1.fastq.gz \
	 	    --SECOND_END_FASTQ ${base_name}.2.fastq.gz
	>>>

	output {
		File fastq1 = "${base_name}.1.fastq.gz"
    		File fastq2 = "${base_name}.2.fastq.gz"
	}


}

task StarGenerateReferences {
	File ref_fasta
	File ref_fasta_index
	File annotations_gtf
	Int? read_length  ## Should this be an input, or should this always be determined by reading the first line of a fastq input

	Int? num_threads
	Int threads = select_first([num_threads, 8])
    


	command <<<
		set -e
		mkdir STAR2_5

		STAR \
		--runMode genomeGenerate \
		--genomeDir STAR2_5 \
		--genomeFastaFiles ${ref_fasta} \
		--sjdbGTFfile ${annotations_gtf} \
		${"--sjdbOverhang "+(read_length-1)} \
		--runThreadN ${threads}

		ls STAR2_5

		tar -zcvf star-HUMAN-refs.tar.gz STAR2_5
	>>>

	output {
		Array[File] star_logs = glob("*.out")
		File star_genome_refs_zipped = "star-HUMAN-refs.tar.gz"
	}

}


task StarAlign {
	File star_genome_refs_zipped
	File fastq1
	File fastq2
	String base_name
	Int? read_length

	Int? num_threads
	Int threads = select_first([num_threads, 8])
	Int? star_mem_max_gb
	Int star_mem = select_first([star_mem_max_gb, 45])
	#Is there an appropriate default for this?
	Int? star_limitOutSJcollapsed

	

	command <<<
		set -e

		tar -xvzf ${star_genome_refs_zipped}

		STAR \
		--genomeDir STAR2_5 \
		--runThreadN ${threads} \
		--readFilesIn ${fastq1} ${fastq2} \
		--readFilesCommand "gunzip -c" \
		${"--sjdbOverhang "+(read_length-1)} \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
		--limitBAMsortRAM ${star_mem+"000000000"} \
		--limitOutSJcollapsed ${default=1000000 star_limitOutSJcollapsed} \
		--outFileNamePrefix ${base_name}.
	>>>

	output {
		File output_bam = "${base_name}.Aligned.sortedByCoord.out.bam"
		File output_log_final = "${base_name}.Log.final.out"
		File output_log = "${base_name}.Log.out"
		File output_log_progress = "${base_name}.Log.progress.out"
		File output_SJ = "${base_name}.SJ.out.tab"
	}


}

task MergeBamAlignment {

    File ref_fasta
    File ref_dict

    File unaligned_bam
    File star_bam
    String base_name

    String gatk_path


    #Using default for max_records_in_ram
 
    command <<<
        ${gatk_path} \
            MergeBamAlignment \
            --REFERENCE_SEQUENCE ${ref_fasta} \
            --UNMAPPED_BAM ${unaligned_bam} \
            --ALIGNED_BAM ${star_bam} \
            --OUTPUT ${base_name}.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --PAIRED_RUN False \
            --VALIDATION_STRINGENCY SILENT
    >>>
 
    output {
        File output_bam="${base_name}.bam"
    }


}

task MarkDuplicates {

 	File input_bam
 	String base_name

  String gatk_path



 	command <<<
 	    ${gatk_path} \
 	        MarkDuplicates \
 	        --INPUT ${input_bam} \
 	        --OUTPUT ${base_name}.bam  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE ${base_name}.metrics
 	>>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 		File metrics_file = "${base_name}.metrics"
 	}


}

## Not validated in GATK4 
task SplitNCigarReads {

	File input_bam
	File input_bam_index
	String base_name
	File interval_list

	File ref_fasta
	File ref_fasta_index
	File ref_dict


    command <<<
    	java -jar /usr/gitc/GATK35.jar \
    		-T SplitNCigarReads \
    		-R ${ref_fasta} \
    		-I ${input_bam} \
    		-o ${base_name}.bam \
    		-rf ReassignOneMappingQuality \
    		-RMQF 255 \
    		-RMQT 60 \
    		-U ALLOW_N_CIGAR_READS
    >>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 	}

}

task SplitNCigarReads_GATK4 {

  File input_bam
  File input_bam_index
  String base_name
  File interval_list

  File ref_fasta
  File ref_fasta_index
  File ref_dict

	String gatk_path
	String docker
        Int preemptible_count

    command <<<
        ${gatk_path} \
                SplitNCigarReads \
                -R ${ref_fasta} \
                -I ${input_bam} \
                -O ${base_name}.bam 
    >>>

        output {
                File output_bam = "${base_name}.bam"
                File output_bam_index = "${base_name}.bai"
        }


}

task BaseRecalibrator {

    File input_bam
    File input_bam_index
    String recal_output_file

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path


    command <<<
        ${gatk_path} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${recal_output_file} \
            -known-sites ${dbSNP_vcf} \
            -known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
    >>>

    output {
        File recalibration_report = recal_output_file
    }


}


task ApplyBQSR {

    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path


    command <<<
        ${gatk_path} \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${base_name}.bam \
            --bqsr-recal-file ${recalibration_report}
    >>>

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

}

task HaplotypeCaller {

	File input_bam
	File input_bam_index
	String base_name

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index

	Int? stand_call_conf

	command <<<
		java -jar /usr/gitc/GATK35.jar \
		    -T HaplotypeCaller \
		    -R ${ref_fasta} \
		    -I ${input_bam} \
		    -L ${interval_list} \
		    -dontUseSoftClippedBases \
		    -stand_call_conf ${default=20 stand_call_conf} \
		    --dbsnp ${dbSNP_vcf} \
		    -o ${base_name}.vcf.gz
	>>>

    output {
        File output_vcf = "${base_name}.vcf.gz"
        File output_vcf_index = "${base_name}.vcf.gz.tbi"
    }

}

task HaplotypeCaller_GATK4 {

	File input_bam
	File input_bam_index
	String base_name

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index

	String gatk_path


	Int? stand_call_conf

	command <<<
		${gatk_path} --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R ${ref_fasta} \
		-I ${input_bam} \
		-L ${interval_list} \
		-O ${base_name}.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
		--dbsnp ${dbSNP_vcf}
	>>>

	output {
		File output_vcf = "${base_name}.vcf.gz"
		File output_vcf_index = "${base_name}.vcf.gz.tbi"
	}

}

task VariantFiltration {

	File input_vcf
	File input_vcf_index
	String base_name

 	File ref_dict
 	File ref_fasta
 	File ref_fasta_index


	command <<<
		 ${gatk_path} \
		    VariantFiltration \
			--R ${ref_fasta} \
			--V ${input_vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O ${base_name}
	>>>

	output {
    	File output_vcf = "${base_name}"
    	File output_vcf_index = "${base_name}.tbi"
	}

	
}

task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    String gatk_path



    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command <<<
        ${gatk_path} --java-options "-Xms2000m"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT=' input_vcfs} \
            --OUTPUT ${output_vcf_name}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }


}

task ScatterIntervalList {

    File interval_list
    Int scatter_count

    String docker
    Int preemptible_count

    command <<<
        set -e
        mkdir out
        java -Xms1g -jar /usr/gitc/picard.jar \
            IntervalListTools \
            SCATTER_COUNT=${scatter_count} \
            SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            UNIQUE=true \
            SORT=true \
            INPUT=${interval_list} \
            OUTPUT=out

        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int(stdout())
    }


}

task ScatterIntervalList_GATK4 {

	File interval_list
	Int scatter_count
	String gatk_path
	String docker
	Int preemptible_count

    command <<<
        set -e
        mkdir out
        ${gatk_path} --java-options "-Xms1g" \
            IntervalListTools \
            --SCATTER_COUNT=${scatter_count} \
            --SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            --UNIQUE=true \
            --SORT=true \
            --INPUT=${interval_list} \
            --OUTPUT=out
	
        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        f = open("interval_count.txt", "w+")
        f.write(str(len(intervals)))
        f.close()
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int("interval_count.txt")
    }

}

task RevertSam {
    File input_bam
    String base_name
    String sort_order

    String gatk_path

    String docker
    Int preemptible_count

    command <<<
        ${gatk_path} \
        	RevertSam \
        	--INPUT ${input_bam} \
        	--OUTPUT ${base_name}.bam \
            --VALIDATION_STRINGENCY SILENT \
        	--ATTRIBUTE_TO_CLEAR FT \
        	--ATTRIBUTE_TO_CLEAR CO \
        	--SORT_ORDER ${sort_order}
    >>>

    output {
        File output_bam = "${base_name}.bam"
    }

   
}
