#### realign main version para pedido de las chicas, mapear a una region acotada de la referencia

import './hnrg-fastq2ubam_realing.wdl' as fastq2ubam 
import './ubam2bwa_realign.wdl' as ubam2bwa
import './bam2gvcf_realig.wdl' as bamtogvcf



task borrado_fastp {
File path1
File path2

#command <<<
#mv ${write_lines(path)}  ${t}.txt;
command <<<
while read -r line; do
rm "$line"
done < "${path1}"
while read -r line2; do
rm "$line2"
done < "${path2}"
>>>

}


task borrado {
String archivo_borrar

command <<<
readlink -f ${archivo_borrar} | xargs rm
>>>
}


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


task groupingBams_bysample_glob {
  String pattern 
  Array[File] array_of_files

  command{
  }

  output {
    Array[File] subArray_input_ubam2gvcf = glob('../inputs/*/*${pattern}*')
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

# task coord_generator {

#   File experiment_lib
#   File chromosome_length
#   Int padding
#   Int merge_tolerance
#   String toolpath
#   String gatk_jar
#   File ref_dict
#   String path_save

  
#   #${toolpath}bedtools2/bin/slopBed -i ${experiment_lib} -g ${chromosome_length} -b ${padding} | sort -k1,1 -k2,2n -V > intervalo_b37_padded_${padding}.bed
    
#   ##java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}.bed -O=intervalo_b37_padded_${padding}_merged_preprocessing.interval_list -SD=${ref_dict}

# ##sort -k1,1 -k2,2n intervalo_b37_padded_${padding}.bed | ${toolpath}bedtools2/bin/mergeBed -d ${ merge_tolerance} > intervalo_b37_padded_${padding}_merged_${ merge_tolerance}.bed
#     #java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed -O=intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list -SD=${ref_dict}
   
#   command <<<
#     #!/bin/bash
#     set -e
#     set -o pipefail
    
#     ${toolpath}bedtools2/bin/slopBed -i ${experiment_lib} -g ${chromosome_length} -b ${padding} | sort -k1,1 -k2,2n -V > intervalo_b37_padded_${padding}.bed 

#     ###merged
     
#     ${toolpath}bedtools2/bin/mergeBed -i intervalo_b37_padded_${padding}.bed -d ${merge_tolerance} > intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed

#     java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}_merged_${merge_tolerance}.bed -O=intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list -SD=${ref_dict}  

#     java -jar ${toolpath}${gatk_jar} BedToIntervalList -I=intervalo_b37_padded_${padding}.bed -O=intervalo_b37_padded_${padding}.interval_list -SD=${ref_dict}
     

#     cp -L intervalo_b37_padded_${padding}.bed ${path_save}
#     cp -L intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list ${path_save}
#     cp -L intervalo_b37_padded_${padding}.interval_list ${path_save}
#   >>>

#   output {
#     File padded_coord = "intervalo_b37_padded_${padding}.bed"
#     #File merged_padded_coord = "intervalo_b37_padded_merged_${merge_tolerance}.bed"
#     File interval_list = "intervalo_b37_padded_${padding}_merged_${merge_tolerance}_preprocessing.interval_list"
#     File eval_interval_list = "intervalo_b37_padded_${padding}.interval_list"
#   }

# }


#task restrict_to_TSO {
#  File padded_interval
#  File generic_exon_coords
#  String toolpath

#  command <<<
#   #!/bin/bash
#    set -e
#    set -o pipefail

#  ${toolpath}bedtools2/bin/intersectBed -wa -a ${generic_exon_coords} -b ${padded_interval} | sort -k1,1 -k2,2n -V | uniq > exon_restricted2interval.bed
#  >>>

#  output {

#    File interval_restricted = "exon_restricted2interval.bed"


#  }

#}


workflow main_workflow {

  ###inputs for fastq2ubam workflows
  File tabulatedSampleFilePaths
  String run_date                   
  String library_name 
  #String platform_unit 
  String platform_name 
  String gatk_jar
  String toolpath
  String ubam_list_name
  String sequencing_center =  "UIT-HNRG" 
  String readlenght
  String platformmodel
 
  File generic_exon_coords
  #File tso_bed

  #inputs for preprocesing ( ubam2gvcf)

  String path_softlink
  String ref_name

 
    
 
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
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  ##### parametros trimmeado fastp
  Int trim_front_fastp = "5" 
  Int trim_tail_fastp = "5"
    
 # Int break_bands_at_multiples_of = "1000000"
 # Int haplotype_scatter_count = "2"

  ##################################
  Int compression_level = "1"
  String java_heap_memory_initial = "128m"

  ########optimization flags
  #String gatk_gkl_pairhmm_implementation = "LOGLESS_CACHING"
  #Int gatk_gkl_pairhmm_threads = "1"
    
 # ####optimizacion haplotypecaller
  #String smith_waterman_implementation = "AVX_ENABLED"
  #Float? contamination = "0"
 # String newqual = "true"


  ##################anotacion funcional
  String reference_version = "GRCh37.75"



 #   String db_annovar #path annovar
 #   File annovar_table_pl #/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/table_annovar.pl
#    File joinPY #/home/hnrg/NGStools/pipeline_wdl/process_vcf/join_vcfs.py


    ###################### inputs para crear intervalo
   # File experiment_lib
   # File chromosome_length
   # Int padding = "0"
   # Int merge_tolerance = "0"

  call mkdir {
    input: 
    path_softlink = path_softlink,


  }



  #call coord_generator {
  #  input:
  #  experiment_lib = experiment_lib,
  #  chromosome_length = chromosome_length,
  #  padding = padding,
  #  merge_tolerance = merge_tolerance,
  #  toolpath = toolpath,
  #  gatk_jar = gatk_jar,
  #  ref_dict = ref_dict,
  #  path_save = path_softlink
  #}


# call restrict_to_TSO {
#  input:
#  padded_interval = coord_generator.padded_coord,
#  generic_exon_coords = generic_exon_coords,
#  toolpath = toolpath
#}



  call fastq2ubam.ConvertPairedFastQsToUnmappedBamWf {  
      input: 
      trim_front = trim_front_fastp, 
      trim_tail = trim_tail_fastp,
      tabulatedSampleFilePaths = tabulatedSampleFilePaths,
      run_date = run_date,                   
      library_name = library_name,   
      platform_name = platform_name,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      ubam_list_name = ubam_list_name,
      sequencing_center = sequencing_center,
      read_lenght = readlenght,
      path_softlink = path_softlink,
      platform_model = platformmodel

  }

  call borrado_fastp {
    input:
     path1 = ConvertPairedFastQsToUnmappedBamWf.p_borrar1,
     path2 = ConvertPairedFastQsToUnmappedBamWf.p_borrar2
  }
 
  call ubam2bwa.ubamtobwa {
   input:
    array_unmapped_bams = ConvertPairedFastQsToUnmappedBamWf.output_ubams,
    #flowcell_unmapped_bams_list = sample,
    ref_name = ref_name,
    compression_level = compression_level,
    java_heap_memory_initial = java_heap_memory_initial,
    bwa_commandline = bwa_commandline,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict,
    ref_amb = ref_amb,
    ref_ann = ref_ann,
    ref_bwt = ref_bwt,
    ref_pac = ref_pac,
    ref_sa = ref_sa,
    gatk_jar = gatk_jar,
    toolpath = toolpath
  } 


 Array[File] array_of_samples_txt = ConvertPairedFastQsToUnmappedBamWf.muestras

 #inputs_bams is an array of files. Each element is a file containing all the aligned and merged bams of a sample.
 scatter (sample_txt in array_of_samples_txt)  {
   #ubamtobwa.output_mergedbam_files
  
   #File flowcell_mapped_bams_listfile = sample_txt
   #Array[File] flowcell_mapped_bams = read_lines(flowcell_mapped_bams_listfile)

   String sample_name = basename(sample_txt, ".txt")

   #####subworkflow de fastq2bwa

   call mkdir_samplename {
    input: 
     path_softlink = path_softlink,
     samplename = sample_name
    }

   call groupingBams_bysample_glob {
     input: 
       pattern = sample_name,
       array_of_files = ubamtobwa.output_mergedbam_files

    }

    call bamtogvcf.bam2gvcf {
     input:
      #base_file_name =  ubamtobwa.nombre_base,
      base_file_name =  sample_name,
      #lib_resctricted = coord_generator.padded_coord,#restrict_to_TSO.interval_restricted,
                                          #### ARTEFACTO PARA PROBAR EL WORKFLOW ANTERIOR, ESTO HAY QUE TRABAJRLO Y ARTICULARLO BIEN. 
      #inputs_ubams = ConvertPairedFastQsToUnmappedBamWf.muestras,
      #uniquesample_name = ConvertPairedFastQsToUnmappedBamWf.samplesnames,
      path_save = mkdir_samplename.path_out_softlink,
      #sample_name = sample_name,
      #flowcell_unmapped_bams_list = sample,
      bams_entrada = groupingBams_bysample_glob.subArray_input_ubam2gvcf,
      #bams_entrada = ubamtobwa.output_mergedbam_files[1],                   #### ARTEFACTO PARA PROBAR EL WORKFLOW ANTERIOR, ESTO HAY QUE TRABAJRLO Y ARTICULARLO BIEN. 
      ##bams_entrada = flowcell_mapped_bams,  #array of basm corresponding to ONE sample. 
      #ref_name = ref_name,
      #unmapped_bam_suffix = ".bam",
      #bwa_commandline = bwa_commandline,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      #wes_calling_interval_list = coord_generator.interval_list,
     #break_bands_at_multiples_of = break_bands_at_multiples_of,
      #haplotype_scatter_count = haplotype_scatter_count,
      compression_level = compression_level,
      #gatk_gkl_pairhmm_implementation = gatk_gkl_pairhmm_implementation,
      #gatk_gkl_pairhmm_threads = gatk_gkl_pairhmm_threads,
      #wgs_calling_interval_list = coord_generator.interval_list,
      #wgs_evaluation_interval_list = coord_generator.interval_list,
      gatk_jar = gatk_jar,
      toolpath = toolpath,
      #smith_waterman_implementation = smith_waterman_implementation,
      #contamination = contamination,
      #newqual = newqual,
      java_heap_memory_initial = java_heap_memory_initial,
      #tso_bed = coord_generator.padded_coord
    } 
}

output {
   Array[File] analysis_ready_bam = bam2gvcf.analysis_ready_bam
   Array[File] analysis_ready_bam_index = bam2gvcf.analysis_ready_bam_index
}
}

