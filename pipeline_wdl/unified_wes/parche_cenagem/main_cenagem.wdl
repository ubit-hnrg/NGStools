task annovar{
    File one_sample_vcf
    File annovar_table_pl
    File convert2annovar
    File annotate_variation
    File variants_reduction


    String db_annovar
    String sample
    
    command<<<
        perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol knownGene -operation  g -nastring . -otherinfo
        #perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol refGene,knownGene,intervar_20180118 -operation  g,g,f -nastring . -otherinfo

    >>>

    output {
        File annovar_vcf = '${sample}.hg19_multianno.vcf'
        File annovar_txt = '${sample}.hg19_multianno.txt'
    }

}

 # prepare final multianno.tsv for deliver (do not forget postprocess InterVar)
 #1) localizo las 3 columnas de Otherinfo para volarlas. 
 #2) modifico el header agregando las columnas que faltan del vcf (remplazo x otherinfo)
 #3) joineo al archivo multinanno las columnas de genotipo de las restantes muestras de  la corrida.

task get_tsv_from_annovar {
    File annovar_txt
    File annovar_vcf
    String sample


    command <<<
    #columnas a cortar (localizando Otherinfo column y las 2 siguientes)
    nl0=$(head -n1 ${annovar_txt}|tr '\t' '\n'|nl|grep 'Otherinfo'|cut -f1)
    nl1=$((nl0 + 1))
    nl2=$((nl0 + 2))


    # meto header (dejando el campo 'Otherinfo' que despues va a aser remplazado por las columnas del vcf original)
    head -n1 ${annovar_txt} > ${sample}.hg19_multianno.tsv
    # vuelo las tres columnas de otherinfo
    tail -n+2 ${annovar_txt}|cut -f$nl0,$nl1,$nl2 --complement >>  ${sample}.hg19_multianno.tsv;
    vcf_header=$(grep '#CH' ${annovar_vcf});

    #remplazo el header
    sed -i "s/Otherinfo/$vcf_header/g" ${sample}.hg19_multianno.tsv;

    #join one multianno tsv file AND joint genotyped vcf. This script (join_vcf.py) also postprocess Intervar columns.
    #change dots by tabs.
    sed -i -e "s|\.	|	|g" ${sample}.hg19_multianno.tsv

    >>>
    output {
        File annovar_tsv =  '${sample}.hg19_multianno.tsv'
    }
}

task build_excell_report{
    File annovar_tsv
    File exon_coverage_report
    #String sample
    String samplename2
    #String original_sample
  
     #/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${sample}.output_xlsx

    command{

       /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${samplename2}_variants.xlsx
   
   }    

    output{
        File excell_report = '${samplename2}_variants.xlsx'
    }
}    

task symlink_important_files {
    File output_to_save
    File output_to_save2

    String path_save
    command{
       set -e -o pipefail

       path_out=$(dirname ${path_save}) 
       
       cp -L ${output_to_save} $path_out

       cp -L ${output_to_save2} $path_out
    }
}




workflow parche {
File list_of_vcf 
File list_of_exon_rep



# for annovar prouposes
String db_annovar
File annovar_table_pl

String toolpath 



Array[File] vcf_in = read_lines(list_of_vcf)
Array[String] path = read_lines(list_of_vcf)
Array[File] exon_coverage_reports = read_lines(list_of_exon_rep)


#Array[Pair[String,File]] vcf_x_path = zip(path, vcf_in)


scatter (idx in length(vcf_in)) {

call annovar {
        input:
            one_sample_vcf =  vcf_in[idx],
            sample = basename(vcf_in[idx], ".final_annot_single_V2.vcf"),
            annovar_table_pl = annovar_table_pl,
            db_annovar = db_annovar
        }



        call get_tsv_from_annovar {
            input:
            annovar_txt = annovar.annovar_txt,
            annovar_vcf = annovar.annovar_vcf,
            sample = basename(vcf_in[idx], ".final_annot_single_V2.vcf")
            
        }

       call build_excell_report {
            input:
            annovar_tsv = get_tsv_from_annovar.annovar_tsv,
            samplename2 = basename(vcf_in[idx], ".final_annot_single_V2.vcf"),
            exon_coverage_report = exon_coverage_reports[idx]
            
        }

       call symlink_important_files {
        input:
        output_to_save = annovar.annovar_vcf,
        output_to_save2 = build_excell_report.excell_report,
        path_save = path[idx]

        }
 
    }



}