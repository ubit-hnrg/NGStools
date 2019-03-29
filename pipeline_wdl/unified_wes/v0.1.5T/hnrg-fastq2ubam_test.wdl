##Copyright Broad Institute, 2018
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - The following metada descriptors per sample:
## ```readgroup   fastq_pair1_file_path   fastq_pair2_file_path   sample_name   library_name   platform_unit   run_date   platform_name   sequecing_center``` 
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs
##
## Cromwell version support 
## - Successfully tested on v32
## - Does not work on versions < v23 due to output syntax
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

# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToUnmappedBamWf {

  File tabulatedSampleFilePaths
  #Array[String] sample_name        # This will be extracted from tabulaltedSampleFilePaths
  #Array[String] fastq_1            # This will be extracted from tabulaltedSampleFilePaths
  #Array[String] fastq_2            # This will be extracted from tabulaltedSampleFilePaths
  #Array[String] readgroup_name     # THIS WILL COMPUTED AUTOMATICALLY FROM FASTQFILES AS: "FLOWCELLID_LANE{NUMBER}"
  String run_date                   
  String library_name 
  String platform_unit 
  String platform_name 
  String sequencing_center
  String gatk_jar
  String toolpath

  String ubam_list_name


  #String? gatk_docker_override
  #String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:latest"])
  #String? gatk_path_override
  #String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])
  #Int? preemptible_attempts

  #obtain array variables from tabulated file
  call read_file_of_tabulated_inputs {
      input:
      tabulatedSampleFilePaths = tabulatedSampleFilePaths
  }
    

  # Convert multiple pairs of input fastqs in parallel
  scatter (i in range(length( read_file_of_tabulated_inputs.array_of_samples ))) {

    call fastp {
        input: 
        sample_name = read_file_of_tabulated_inputs.array_of_samples[i],
        R1_fastq_gz = read_file_of_tabulated_inputs.array_of_R1_files[i],
        R2_fastq_gz = read_file_of_tabulated_inputs.array_of_R2_files[i],
        toolpath = toolpath
    }

    call get_read_group_name{
        input:
        fastq_1 = fastp.fastq_cleaned_R1
    }

    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        sample_name = read_file_of_tabulated_inputs.array_of_samples[i],
        fastq_1 = fastp.fastq_cleaned_R1,
        fastq_2 = fastp.fastq_cleaned_R2,

        rgfile = get_read_group_name.rgname,
        library_name = library_name,
        platform_unit = platform_unit,
        run_date = run_date,   ### ojo aca es la fecha de coorrida del seqcuenciador!!! Esto puede que deba ser una columna mas a parsear en el archivo tabulado y convertirse en array
        platform_name = platform_name,
        sequencing_center = sequencing_center,
        gatk_jar = gatk_jar,
        toolpath = toolpath
        #docker = gatk_docker,
        #preemptible_attempts = preemptible_attempts
    }
  
  
#call Creando_inputs {
 #input:
  #ubam_path = PairedFastQsToUnmappedBAM.output_bam,
  #ubam_name = read_file_of_tabulated_inputs.array_of_samples[i],

 #}

} ## fin scatter
 
  #Create a file with a list of the generated ubams
  call CreateFoFN {
    input:
      array_of_files = PairedFastQsToUnmappedBAM.output_ubam,
      fofn_name = ubam_list_name,
      #docker = gatk_docker
  }

  
  call Create_inputs_for_preprocesing {
    input:
      ubams_paths = CreateFoFN.fofn_list,
      bams_sample_names = read_file_of_tabulated_inputs.samplenames,

  }

  call path_borrado {
#
      input:
     path1 = fastp.fastq_cleaned_R1,
     path2 = fastp.fastq_cleaned_R2
 }

  # Outputs that will be retained when execution is complete
  output {
    #String path_borrado = write_lines(fastp.fastq_cleaned_R1)
    File p_borrar1 = path_borrado.path_borrar1 
    File p_borrar2 = path_borrado.path_borrar2
    Array[File] output_ubams = PairedFastQsToUnmappedBAM.output_ubam
    Array[String] output_ubams_sample_names =  read_file_of_tabulated_inputs.array_of_samples
    File unmapped_ubam_list = CreateFoFN.fofn_list 
    File samplesnames = read_file_of_tabulated_inputs.unique_samples ##samples names unicos
    Array[File] muestras  =  Create_inputs_for_preprocesing.ubam_samples ###array de samples: S1,S2,..,Sn
  }


}

# TASK DEFINITIONS

####################################################################################################################
#the expected file must be 3 columns tab separated: One line per multiplexed lane (and flowcell if would be the case)
#####################################################################################################################

#sample_1   fastq_L1_R1__abspath    fastq_L1_R2_abspath
#sample_1   fastq_L2_R1_abspath    fastq_L2_R2_abspath
#sample_1   fastq_L3_R1_abspath    fastq_L3_R2_abspath
#sample_2   fastq_L1_R1_abspath    fastq_L1_R2__abspath

task read_file_of_tabulated_inputs {
    File tabulatedSampleFilePaths
    Array[String] my_lines = read_lines(tabulatedSampleFilePaths)
#    String sep = '|'

    command <<<
        cut -f1 -d'|' ${tabulatedSampleFilePaths} >  samples.list
        cut -f2 -d'|' ${tabulatedSampleFilePaths} >  R1_fastq.list
        cut -f3 -d'|' ${tabulatedSampleFilePaths} >  R2_fastq.list
        cut -f1 -d'|' ${tabulatedSampleFilePaths}| sort| uniq > unique_samples.list

    >>>
        #cut -f1 -d'|' ${tabulatedSampleFilePaths}| sort| uniq >  unique_samples.list
    output {
    	File samplenames = 'samples.list'
        Array[String] array_of_samples = read_lines('samples.list')
    	Array[File] array_of_R1_files = read_lines('R1_fastq.list')
        Array[File] array_of_R2_files = read_lines('R2_fastq.list')
        File unique_samples = 'unique_samples.list'

    }  
}


## cleaning fastq files
task fastp {

String sample_name
File R1_fastq_gz
File R2_fastq_gz
String R1_stripped_basename = basename(R1_fastq_gz, ".fastq.gz")
String R2_stripped_basename = basename(R2_fastq_gz, ".fastq.gz")
String toolpath


command {
    ${toolpath}fastp -i ${R1_fastq_gz} -I ${R2_fastq_gz} -o ${R1_stripped_basename}_cleaned.fastq.gz -O ${R2_stripped_basename}_cleaned.fastq.gz -h report.html -j report.json
}

output {
    File fastq_cleaned_R1 = "${R1_stripped_basename}_cleaned.fastq.gz"
    File fastq_cleaned_R2 = "${R2_stripped_basename}_cleaned.fastq.gz"
    
    }

}



task get_read_group_name{
    File fastq_1
    command <<<
        set -e
        rgname=$(zcat ${fastq_1} |head -n1| cut -f3,4 -d':' |sed -e 's/:/.Lane/g')
        echo $rgname > rgname.txt    
    >>>
    output{
        File rgname = 'rgname.txt'
    }
}


# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  # Command parameters
  String sample_name
  File fastq_1
  File fastq_2
  File rgfile # this expects an one line file containing the readgroupname
  String readgroup_name  = read_lines(rgfile)[0]
  String library_name
  String platform_unit
  String run_date
  String platform_name
  String sequencing_center


  String gatk_jar
  String toolpath

  command {
    
    #rgname=$(zcat ${fastq_1} |head -n1| cut -f3,4 -d':' |sed -e 's/:/.Lane/g')    
    #Flowcell=$(zcat ${fastq_1} |head -n1| cut -f3 -d':' )   
    #cat $readgroup_name > rgname.txt 
    #readgroup_name = $rgname

    java -Xmx8G -jar ${toolpath}${gatk_jar} \
    FastqToSam \
    --FASTQ ${fastq_1} \
    --FASTQ2 ${fastq_2} \
    --OUTPUT ${sample_name}_${readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ${readgroup_name} \
    --SAMPLE_NAME ${sample_name} \
    --LIBRARY_NAME ${library_name} \
    --PLATFORM_UNIT ${platform_unit} \
    --RUN_DATE ${run_date} \
    --PLATFORM ${platform_name} \
    --SEQUENCING_CENTER ${sequencing_center}
  }
#
  output {
    File output_ubam = "${sample_name}_${readgroup_name}.unmapped.bam"
  }
}

task CreateFoFN {
  # Command parameters
  Array[String] array_of_files
  String fofn_name
  #Map [String, String] lista = {"array_of_files":"fofn_name"} 
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list \
   
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }
}

task Creando_inputs {
 Array[String] ubam_path 
 String ubam_name
 

 command {
    mv ${write_lines(ubam_path)}  ${ubam_name}.txt \
   
  }
output {
  File ubamsample = "${ubam_name}.txt"
 }
 }

task Create_inputs_for_preprocesing {
 File bams_sample_names
 File ubams_paths 
# Array[File] = []

command <<<  
python <<CODE 

with open("${bams_sample_names}", "r") as sf:
    samples = sf.readlines()
    samples =[i.strip('\n') for i in samples]
    if samples[-1]=='':
        samples = samples[:-1]
        

with open("${ubams_paths}", "r") as ubf:
    ubams = ubf.readlines()
    ubams =[i.strip('\n') for i in ubams]
    if ubams[-1]=='':
        ubams = ubams[:-1]
      
open_files = []
for i in range(len(samples)):
    sample = samples[i]
    ubam = ubams[i]
    
    filename ='%s.txt'%sample
    if sample not in open_files:
        with open(filename,'w') as f:
            f.write("%s\n"%ubam)
            open_files.append(sample)
        f.close()
    else:
        with open(filename,'a') as f:
            f.write("%s\n"%ubam)
        f.close()

CODE
>>>

output {

    Array[File] ubam_samples = glob("*.txt")
}


}

task path_borrado {
Array[String] path1
Array[String] path2 
String temp1 = "temp1"
String temp2 = "temp2"

command <<<
mv ${write_lines(path1)}  ${temp1}.txt
mv ${write_lines(path2)}  ${temp2}.txt

>>>

output {
File path_borrar1 = "${temp1}.txt"
File path_borrar2 = "${temp2}.txt"

}


}
