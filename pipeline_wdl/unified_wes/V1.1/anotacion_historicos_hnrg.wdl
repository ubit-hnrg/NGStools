
task bplat_annot {

    File input_vcf
    File funcion_lua
    File config_file_vcfanno
    String sample_name
    String fecha_poblacion_hnrg = "hoy"
    String toolpath

    

    command {
        ${toolpath}/vcfanno_linux64 -p 4 -lua ${funcion_lua} ${config_file_vcfanno} ${input_vcf} > ${sample_name}.${fecha_poblacion_hnrg}.vcf

    }
output {
    File out_vcfanno = "${sample_name}.${fecha_poblacion_hnrg}.vcf"
}

}

task get_dir {
String path1

command {
    set -e -o pipefail
    dirname ${path1}

}
output {
   String path = stdout()
}


}

task copy_to_data {
    File output_to_save
    String path_save
    command{
       set -e -o pipefail

       cp -L ${output_to_save} ${path_save}
    }
}



workflow FuncionalAnnotation {

File list_of_vcf 
String toolpath 
#String path_save
File funcion_lua
File config_file_vcfanno



Array[File] vcf_in = read_lines(list_of_vcf)
Array[String] path = read_lines(list_of_vcf)
Array[Pair[File,String]] vcf_x_path = zip(vcf_in, path)


scatter (pairs in vcf_x_path) {


 String sample_name1 = basename(pairs.left, ".vcf")
 #String path_save = basename(vcf,sample_name1+".final_annot.vcf")

 #annotate for bplat format
 call bplat_annot {

 input:
    input_vcf = vcf,
    funcion_lua = funcion_lua, 
    config_file_vcfanno = config_file_vcfanno,
    sample_name = sample_name1,
    toolpath = toolpath

    
 }

 call get_dir {
     input:
     path1 = pairs.right

 }


 call copy_to_data {
        input:
        output_to_save = bplat_annot.out_vcfanno,
        path_save = get_dir.path
    }

 }
}