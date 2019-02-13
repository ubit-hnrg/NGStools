task read_fof {
    File fof
    Array[File] my_files = read_lines(fof)
    #Array[File] my_files=['a','b','c']
    command {
        touch ${fof}.csv
    }
    
    output {
    	Array[File] array_of_files = my_files
    }  
}



# fofN: file of filesNames

task write_FoFN {
  # Command parameters
  Array[String] array_of_files
  String fofn_name
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }

}


workflow fof_wf {
    File file_of_files
    fof_name = 'WRITED'

    call read_fof {
        input:
        fof = file_of_files
    }
    

  #Create a file with a list of the generated ubams
  call write_FoFN {
    input:
      array_of_files = read_fof.array_of_files,
      fofn_name = 'out_test'
  }

}

#    output {
#    	Array[File] array_output = fof_usage_task.array_of_files
#    }

