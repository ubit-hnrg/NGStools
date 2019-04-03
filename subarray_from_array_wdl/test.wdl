workflow test {
    Array[File] files = ['/home/ariel/bla_ooh.txt','/home/ariel/ooh_ooh.txt','/home/ariel/pepe_blabla.esaaaa']
    
    call subset_array_glob_elegante{
        input:
        array_of_files = files,
        pattern = 'bla'
    }
    
}



task subset_array_glob_elegante {
  String pattern = 'bla'
  Array[File] array_of_files

  command{
  }

  output {
    Array[File] subArray = glob('../inputs/*/*${pattern}*')
  }
}




task subset_array_glob_cabeza {
  String pattern = 'bla'
  Array[File] array_of_files

  command{
      for element in ${sep=' ' array_of_files}  ; 
      do
      ln -s $element .
      done;
  }

  output {
    Array[File] subArray = glob('*${pattern}*')
  }
}
