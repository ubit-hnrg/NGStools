import './processMultisampleVCF4annotation.wdl' as splitVCF
import './anotaciones_hnrg4annotation.wdl' as anotaciones


task copy_important_files {
    File output_to_save
    String path_save
    #ln -s ${output_to_save} ${path_save}
    command{
       cp -L ${output_to_save} ${path_save}
    }
}

task build_excell_report{
    File annovar_tsv
    #File exon_coverage_report
    #String sample
    String samplename2
    #String original_sample
  
     #/home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${sample}.output_xlsx

    command{

       /home/hnrg/NGStools/pipeline_wdl/qualityControl/make_excel_report.py ${annovar_tsv}:Variants  ${samplename2}_variants.xlsx
   
   }    

    output{
        File excell_report = '${samplename2}_variants.xlsx'
    }
}  




workflow functional_annotation {

  ###toolpaths
  String gatk_jar
  String toolpath
  String path_softlink
  
  String db_annovar #path annovar
  File annovar_table_pl #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
  File joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py
    
  String reference_version = "GRCh37.75"

  File TSO 
  #File reportes_by_exon #"${name}_coverage_statistics_by_exon.tsv"
  #File exon_coords
  File tso_bed
  String version = "v0.1.9"

 #Array[File] inputs_scatter_by_exon = read_lines(reportes_by_exon)



#scatter (reportebyexon in inputs_scatter_by_exon){



 call splitVCF.processJointVCF {
     input:
     multisampleVCF = TSO, ###inputs TSO*.vcf.gz
     #array_path_save = mkdir_samplename.path_out_softlink, #### path/softlink + samplename
     version = version,

     toolpath = toolpath, 
     region_padded_bed = tso_bed,##Tso_bed
     path_softlink = path_softlink,

    # for annovar prouposes
    db_annovar = db_annovar,#path annovar
    annovar_table_pl = annovar_table_pl, #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
    joinPY = joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py
 }

Array[String] path_save = processJointVCF.paths_copy


Array[File] Tsv_annovar = processJointVCF.annovar_tsv_out
    scatter (idx in range(length(Tsv_annovar))){
       
       String sample = basename(Tsv_annovar[idx],".multianno_multisample.tsv")
       #String samplename2 = basename(inputs_scatter_by_exon[idx],"_coverage_statistics_by_exon.tsv")
       
       #if(sample==samplename2){
       call build_excell_report {
            input:
            annovar_tsv = Tsv_annovar[idx],
            samplename2 = sample
           # exon_coverage_report = inputs_scatter_by_exon[idx]
            
           }
         # }
      

    }

Array[File?] reporte_variantes = build_excell_report.excell_report
Array[Pair[String,File?]] samples_by_variant = zip (path_save, reporte_variantes)
  scatter (pairs in samples_by_variant) {
    call copy_important_files as build_excell_reportbyvariants{
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
    }
  }



 ####anotaciones funcionales
Array[File] vcf_individuales = processJointVCF.individual_vcfs_annovar
Array[Pair[String,File]] vcf_x_path = zip (path_save, vcf_individuales)
  scatter (vcf in vcf_x_path) {
    call anotaciones.FuncionalAnnotation {
        input:
        input_vcf = vcf.right,
        path_save = vcf.left,
        toolpath = toolpath,
        version1 = version,
        samplename1 = basename(vcf.right,".hg19_multianno.vcf"),
        java_heap_memory_initial = "12g",
        reference_version = reference_version,

      }
    }


####fin workflow
}