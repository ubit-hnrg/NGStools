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
    File fastp_output_1 = "${R1_stripped_basename}_cleaned.fastq.gz"
    File fastp_output_2 = "${R2_stripped_basename}_cleaned.fastq.gz"
    }

}


########################################
########  fastq.gz to ubam    ##########
########################################

task fastq2ubam {
String sample_name
File R1_fastq_cleaned_gz
File R2_fastq_cleaned_gz
String toolpath = "/home/ariel/tools/"
String gatk4_jar
String Flowcell
Int Lane_Number
String LIBRARY
String PLATFORM = "illumina"

command{
   command #<<<
    #set -e 
    #touch ${PLATFORM}.txt    
    #touch ${gatk4_jar}.txt
    java -Xmx8G -jar ${toolpath}${gatk4_jar} FastqToSam --FASTQ=${R1_fastq_cleaned_gz} --FASTQ2=${R2_fastq_cleaned_gz} --OUTPUT=${sample_name}_${Flowcell}_Lane${Lane_Number}_u.bam --READ_GROUP_NAME=${Flowcell}_Lane${Lane_Number} --SAMPLE_NAME=${sample_name} --LIBRARY_NAME=${LIBRARY} --PLATFORM=${PLATFORM}
    #>>>
}

output {
    File fastp_output_1 = "${sample_name}_${Flowcell}_Lane${Lane_Number}_u.bam"
    }


}


workflow prefiltering {

String sample_name
File R1_fastq_gz
File R2_fastq_gz
String gatk4_jar

call fastp {
    input: 
    sample_name = sample_name,
    R1_fastq_gz = R1_fastq_gz,
    R2_fastq_gz = R2_fastq_gz
}

call fastq2ubam {
    input: 
    sample_name = sample_name,
    R1_fastq_cleaned_gz = R1_fastq_gz,
    R2_fastq_cleaned_gz = R2_fastq_gz,
    gatk4_jar = gatk4_jar
}


}

