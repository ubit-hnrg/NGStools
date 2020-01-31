task hnrg_freq {

    File input_vcf
    #File funcion_lua
    File config_file_vcfanno
    String sample_name
    String hnrg = "HNRG"
    String toolpath

    

    command {
        ${toolpath}/vcfanno_linux64 -p 4 ${config_file_vcfanno} ${input_vcf} > ${sample_name}.${hnrg}.vcf

    }
output {
    File out_vcfanno = "${sample_name}.${hnrg}.vcf"
}

}


task copy_to_data {
    File output_to_save
    String path_vcf_original
    command{
       set -e -o pipefail

       path_out=$(dirname ${path_vcf_original})


       cp -L ${output_to_save} $path_out
    }
}



workflow freq_hnrg {

File list_of_vcf 
String toolpath 
#String path_save
File funcion_lua
File config_file_vcfanno



Array[File] vcf_in = read_lines(list_of_vcf)
Array[String] path = read_lines(list_of_vcf)
Array[Pair[String,File]] vcf_x_path = zip(path, vcf_in)


scatter (pairs in vcf_x_path) {


 String sample_name1 = basename(pairs.right, ".vcf")
 #String path_save = basename(vcf,sample_name1+".final_annot.vcf")

 #annotate for bplat format
 call hnrg_freq {

 input:
    input_vcf = pairs.right, 
    config_file_vcfanno = config_file_vcfanno,
    sample_name = sample_name1,
    toolpath = toolpath

    
 }


 call copy_to_data {
        input:
        output_to_save = hnrg_freq.out_vcfanno,
        path_vcf_original = pairs.left
    }

 }
}