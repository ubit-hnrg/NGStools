#### main Version 1

import './jointgenotype_single.wdl' as genotypeGVCF
import './anotaciones_hnrg_single.wdl' as anotacionesSingle




task mkdir {
    String path_softlink

    command{
        mkdir -p ${path_softlink}
    }
}

task mkdir_samplename {
    String path_softlink
    String samplename

    command{
        mkdir -p ${path_softlink}${samplename}
    }

    output {
        String path_out_softlink = "${path_softlink}" + "${samplename}"+"/"
}
}

task copy2data {
    File output_to_save
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
    }
}




task symlink_important_files {
    File output_to_save
    String path_save
    #ln -s ${output_to_save} ${path_save}
    command{
       cp -L ${output_to_save} ${path_save}
    }
}

task coord_generator {

  File experiment_lib
  File chromosome_length
  Int padding
  Int merge_tolerance
  String toolpath
  String gatk_jar
  File ref_dict
  String path_save

  
  command <<<
    #!/bin/bash
    set -e
    set -o pipefail
    
    ${toolpath}bedtools2/bin/slopBed -i ${experiment_lib} -g ${chromosome_length} -b ${padding} | sort -k1,1 -k2,2n -V > intervalo_b37_padded_${padding}.bed 

    ###merged
     
    ${toolpath}bedtools2/bin/mergeBed -i intervalo_b37_padded_${padding}.bed -d ${merge_tolerance} > intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed

    java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed -O=intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list -SD=${ref_dict}  

    java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}.bed -O=intervalo_b37_padded_${padding}.interval_list -SD=${ref_dict}
     

    cp -L intervalo_b37_padded_${padding}.bed ${path_save}
    cp -L intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list ${path_save}
    cp -L intervalo_b37_padded_${padding}.interval_list ${path_save}
  >>>

  output {
    File padded_coord = "intervalo_b37_padded_${padding}.bed"
    #File merged_padded_coord = "intervalo_b37_padded_merged_${merge_tolerance}.bed"
    File interval_list = "intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list"
    File eval_interval_list = "intervalo_b37_padded_${padding}.interval_list"
  }

}




workflow main_workflow {

  ###inputs for fastq2ubam workflows
  #File tabulatedSampleFilePaths
  String gatk_jar
  String toolpath
  String ubam_list_name
  File gvcf_list_of_files
  File gvcf_list_of_idx
 
  File generic_exon_coords
  String callset_name
  #File tso_bed

  #inputs for preprocesing ( ubam2gvcf)

  String path_softlink
  String ref_name
  String reference_version


 
    
 
  ########################  command line para el bwa ##################
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 4 -Y"
  
  ########## referencia
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
    
  File dbSNP_vcf
  File dbSNP_vcf_index

    
  ###################### inputs para crear intervalo
  File experiment_lib
  File chromosome_length
  Int padding = "100"
  Int merge_tolerance = "0"

  call mkdir {
    input: 
    path_softlink = path_softlink,
  }



  call coord_generator {
    input:
    experiment_lib = experiment_lib,
    chromosome_length = chromosome_length,
    padding = padding,
    merge_tolerance = merge_tolerance,
    toolpath = toolpath,
    gatk_jar = gatk_jar,
    ref_dict = ref_dict,
    path_save = path_softlink
  }


#   call restrict_to_TSO {
#     input:
#     padded_interval = coord_generator.padded_coord,
#     generic_exon_coords = generic_exon_coords,
#     toolpath = toolpath
#   }


##cantidad de gvcfs
  Int cantidad_gvcf = length(read_lines(gvcf_list_of_files))
  Array[File] array_of_gvcfs = read_lines(gvcf_list_of_files)
  Array[File] array_of_gvcfs_idx = read_lines(gvcf_list_of_idx)

  Array[Pair[File,File]] gvcf_x_idx = zip(array_of_gvcfs,array_of_gvcfs_idx)


 #inputs_bams is an array of files. Each element is a file containing all the aligned and merged bams of a sample.
 scatter (gvcf in gvcf_x_idx)  {


   String sample_name = basename(gvcf.left, ".g.vcf.gz")

   call mkdir_samplename {
    input: 
     path_softlink = path_softlink,
     samplename = sample_name
    }


   ##########################################llamada workflow Joint Genotyping
   call genotypeGVCF.singleGenotypeGVCFs {
    input:
    num_gvcfs= cantidad_gvcf,
    eval_interval_list   = coord_generator.eval_interval_list,
    array_path_save = mkdir_samplename.path_out_softlink,
    dbSNP_vcf = dbSNP_vcf,
    dbSNP_vcf_index = dbSNP_vcf_index,
    callset_name = sample_name,
    ref_fasta = ref_fasta,
    ref_fasta_index =ref_fasta_index,
    ref_dict = ref_dict,
    gatk_jar = gatk_jar,
    toolpath = toolpath,
    sample_names = sample_name,
    input_gvcfs = gvcf.left,
    input_gvcfs_indices = gvcf.right
    #region_padded_bed = coord_generator.padded_coord
    }
 
 call anotacionesSingle.FuncionalAnnotationSingle {
        input:
        input_vcf = singleGenotypeGVCFs.restricted_vcf,
        path_save = mkdir_samplename.path_out_softlink,
        toolpath = toolpath,
        samplename1 = basename(gvcf.left,"_restricted.vcf"),
        java_heap_memory_initial = "12g",
        reference_version = reference_version,
        version_db = "single_V2"

      }


 
  }
    # Outputs that will be retained when execution is complete
 output {

   
   Array[File?] outputvcf = singleGenotypeGVCFs.outputvcf
   Array[File?] outputvcfindex =  singleGenotypeGVCFs.outputvcfindex
   Array[File?] detail_metrics_file =  singleGenotypeGVCFs.metrica1
   Array[File?] summary_metrics_file = singleGenotypeGVCFs.metrica2

  
  }



  
}

