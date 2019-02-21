#############################
#### General description ####
#
# This script perform two main tasks
# 1) Prefiltering reads based on 
#   * Read's Quality
#   * Read's Length
#   * Trimming of Adapters (unsupervised, or optionally providing an adapter sequence)
#   * Quality Report of input reads.
#
# 2) Ubam generation, one per fastq file, (one per Flowcell and Lane)


############################
##### Mandatory Inputs #####
#
# samplename
# array of foward reads (R1.fastq.gz fastqfiles)
# array of reverse reads (R2.fastq.gz fastqfiles)

#for optional inputs see json file.

###########################
#### #expected output  ####
#
## one ubam file per each PE files 
#  i.e. : samplename_flowcell-ID_Lane-N.ubam
#



####################################################################################################################
#the expected file must be 3 columns tab separated: One line per multiplexed lane (and flowcell if would be the case)
#####################################################################################################################

#sample_1   fastq_L1_R1__abspath    fastq_L1_R2_abspath
#sample_1   fastq_L2_R1_abspath    fastq_L2_R2_abspath
#sample_1   fastq_L3_R1_abspath    fastq_L3_R2_abspath
#sample_2   fastq_L1_R1_abspath    fastq_L1_R2__abspath

task read_file_of_tabulated_inputs {
    File ffields
    Array[String] my_lines = read_lines(ffields)
#    String sep = '|'

    command <<<
        cut -f1  ${ffields} >  samples.list
        cut -f2  ${ffields} >  R1_fastq.list
        cut -f3  ${ffields} >  R2_fastq.list

    >>>
    
    output {
    	Array[String] array_of_samples = read_lines('samples.list')
    	Array[File] array_of_R1_files = read_lines('R1_fastq.list')
        Array[File] array_of_R2_files = read_lines('R2_fastq.list')

    }  
}


########################################
## Trimming filtering and quality report
########################################

task fastp {

String sample_name
File R1_fastq_gz
File R2_fastq_gz
String R1_stripped_basename = basename(R1_fastq_gz, ".fastq.gz")
String R2_stripped_basename = basename(R2_fastq_gz, ".fastq.gz")


command {
    fastp -i ${R1_fastq_gz} -I ${R2_fastq_gz} -o ${R1_stripped_basename}_cleaned.fastq.gz -O ${R2_stripped_basename}_cleaned.fastq.gz -h report.html -j report.json
}

output {
    File fastq_cleaned_R1 = "${R1_stripped_basename}_cleaned.fastq.gz"
    File fastq_cleaned_R2 = "${R2_stripped_basename}_cleaned.fastq.gz"
    }

}


########################################
########  fastq.gz to ubam    ##########
########################################

task fastq2ubam {
    String sample_name
    File R1_fastq_cleaned_gz
    File R2_fastq_cleaned_gz
    String toolpath
    String gatk4_jar
    String Flowcell
    Int Lane_Number
    String LIBRARY
    String PLATFORM = "illumina"

    command <<<
        set -e
        rgid=$(zcat ${R1_fastq_cleaned_gz} |head -n1| cut -f3,4 -d':' |sed -e 's/:/.Lane/g')
        java -Xmx8G -jar ${toolpath}${gatk4_jar} FastqToSam \
        --FASTQ=${R1_fastq_cleaned_gz} \
        --FASTQ2=${R2_fastq_cleaned_gz} \
        --OUTPUT=${sample_name}_${Flowcell}_Lane${Lane_Number}_u.bam \
        --READ_GROUP_NAME=$rgid --SAMPLE_NAME=${sample_name} \
        --LIBRARY_NAME=${LIBRARY} \
        --PLATFORM=${PLATFORM}
    >>>

    output {
        File ubam_1 = "${sample_name}_${Flowcell}_Lane${Lane_Number}_u.bam"
        }

    }

########################################
########  Validate ubam    ##########
########################################

task validate_ubam {

	String sample_name
	File input_bam
	String toolpath
    String gatk4_jar



	command {
    
  	java -Xmx8G -jar ${toolpath}${gatk4_jar} ValidateSamFile -I=${input_bam} -M SUMMARY > >(tee ${sample_name}.Validatesamfile.stderr.log) | tail -n1
	}

	output {    
    	String flag = read_string(stdout())
    	File stderr_log = "${sample_name}.Validatesamfile.stderr.log"
	}

}


workflow prefiltering {

    String sample_name
    File R1_fastq_gz
    File R2_fastq_gz
    String gatk4_jar
    String toolpath

    call fastp {
        input: 
        sample_name = sample_name,
        R1_fastq_gz = R1_fastq_gz,
        R2_fastq_gz = R2_fastq_gz
    }

    call fastq2ubam {
        input: 
        sample_name = sample_name,
        R1_fastq_cleaned_gz=fastp.fastq_cleaned_R1,
        R2_fastq_cleaned_gz=fastp.fastq_cleaned_R2,
        toolpath = toolpath,
        gatk4_jar = gatk4_jar

    }

    call validate_ubam{
        input:
        sample_name = sample_name,
        input_bam = fastq2ubam.ubam_1,
        toolpath = toolpath, 
        gatk4_jar = gatk4_jar
    }

}

