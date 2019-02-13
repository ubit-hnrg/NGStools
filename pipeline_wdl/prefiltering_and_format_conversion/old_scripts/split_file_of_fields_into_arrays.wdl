
task read_file_of_fields {
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


task get_rgid_from_fastq {
    File fastqfile
    command <<<
        set -e
        rgid=$(zcat ${R1_fastq_cleaned_gz} |head -n1| cut -f3,4 -d':' |sed -e 's/:/.Lane/g')
    >>>

    output{
        String RGid = rgid
    }
}


#Create a file.list with the content of an array

task write_array {
  # Command parameters
  Array[String] array_to_write
  String file_name
  
  command {
    mv ${write_lines(array_to_write)}  ${fofn_name}.list
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }

}


workflow parse_fields {
    File file_of_fields

    call read_file_of_fields {
        input:
        ffields = file_of_fields
    }
    

  call write_array as write_samples{
    input:
      array_to_write = read_file_of_fields.array_of_samples,
      file_name = 'samples'
  }


  }


}

#    output {
#    	Array[File] array_output = fof_usage_task.array_of_files
#    }

