##### V38 pipeline UIT marzo 2024
## solo anotacion, toma un vcf ya listo para anotar.
## 

import './only_annot_hnrg_38.wdl' as anotaciones_only #5


########MAIN

workflow main_workflow {

  ###inputs for fastq2ubam workflows
  

  String? pipeline_version = "V2" ###no me lo toma el input como default


  ####metadata
  String run_date                   
  String library_name 
  String platform_name = "Illumina"
  String platformmodel = "NextSeq-500"
  String sequencing_center =  "UIT-HNRG" 
  String readlenght 
  String ubam_list_name = "ubamfiles"

  ###GATK
  String gatk_jar = "gatk-package-4.5.0.0-local.jar"
  String toolpath = "/home/hnrg/HNRG-pipeline-V0.1/tools/"
  String ngs_toolpath = "/home/hnrg/NGStools"


  ###save location
  String path_softlink

  String sample_name = basename(input_vcf, ".V38_restricted.vcf")
  String path = dirname(input_vcf)
  
  String java_heap_memory_initial = "128m"


  ##################anotacion funcional
  String reference_version = "GRCh38.mane.1.2.refseq"


  ##annovar
  ###############agregar path a la 38
  String db_annovar = "/data/new_dbs/annovar/hg38/humandb" #path annovar 
  File input_vcf 

    
    ###esto creo q vuela
    File annovar_table_pl #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
    File joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py

    call anotaciones_only.FuncionalAnnotationSingle {
        input:
        input_vcf = input_vcf, #sin annovar del genotipado , 
        path_save = path_softlink,
        toolpath = toolpath,
        samplename1 = sample_name,
        java_heap_memory_initial = "1g",
        pipeline_version = pipeline_version,
        reference_version = reference_version
        
      }
     


   } 