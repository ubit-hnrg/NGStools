#### main anotacion Version 1

import './anotaciones_hnrg_single.wdl' as anotacionesSingle


workflow main_workflow {

File list_of_vcf 
String toolpath 



Array[File] vcf_in = read_lines(list_of_vcf)
Array[String] path = read_lines(list_of_vcf)
Array[Pair[String,File]] vcf_x_path = zip(path, vcf_in)


scatter (pairs in vcf_x_path) {
 
 call anotacionesSingle.FuncionalAnnotationSingle {
        input:
        input_vcf = pairs.left,
        path_vcf_original = pairs.right,
        toolpath = toolpath,
        samplename1 = basename(pairs.left,".vcf"),
        java_heap_memory_initial = "12g"
        
        
      }
  }

}

