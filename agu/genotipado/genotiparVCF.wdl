## This WDL implements the joint discovery and VQSR filtering portion of the GATK 
## Best Practices (June 2016) for germline SNP and Indel discovery in human 
## whole-genome sequencing (WGS) and exome sequencing data.

task ImportGVCFs {
  Array[String] sample_names
  Array[File] input_gvcfs
  Array[File] input_gvcfs_indices
  String interval

  String workspace_dir_name

  String path_herramientas
  
  Int batch_size

  command <<<
    set -e
    set -o pipefail
    
    python << CODE
    gvcfs = ['${sep="','" input_gvcfs}']
    sample_names = ['${sep="','" sample_names}']

    if len(gvcfs)!= len(sample_names):
      exit(1)

    with open("inputs.list", "w") as fi:
      for i in range(len(gvcfs)):
        fi.write(sample_names[i] + "\t" + gvcfs[i] + "\n") 

    CODE

    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.

    #padding (-ip) = 500, yo puse 100
    java -Xmx4g -Xms4g -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${interval} \
    --sample-name-map inputs.list \
    -ip 100

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

  >>>
 
  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}


task GenotypeGVCFs {
  File workspace_tar
  String interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbsnp_vcf
  File dbsnp_vcf_index

  String path_herramientas

  command <<<
    set -e

    tar -xf ${workspace_tar}
    WORKSPACE=$( basename ${workspace_tar} .tar)

    java -Xmx5g -Xms5g -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbsnp_vcf} \
     -G StandardAnnotation \
     --only-output-calls-starting-in-intervals true \
     -new-qual true\
     -V gendb://$WORKSPACE \
     -L ${interval}
  >>>

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename
  String path_herramientas

  command {
    set -e

    java -Xmx3g -Xms3g -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      VariantFiltration \
      -filter "ExcessHet > ${excess_het_threshold}" \
      --filter-name "ExcessHet" \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    java -Xmx3g -Xms3g -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      MakeSitesOnlyVcf \
      -I ${variant_filtered_vcf_filename} \
      -O ${sites_only_vcf_filename}

  }
 
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  }
}




task IndelsVariantRecalibrator {
  String recalibration_filename
  String tranches_filename

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  File sites_only_variant_filtered_vcf
  File sites_only_variant_filtered_vcf_index

  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File dbsnp_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf_index
  String path_herramientas
  

  command {
    java -Xmx12g -Xms12g -jar ${path_herramientas}/gatk-package-4.0.8.1-local.jar \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode INDEL \
      --max-gaussians 4 \
      -resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPoly_resource_vcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbsnp_resource_vcf}
  }

  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
  }
}





######################################################################################################################workflow
workflow JointGenotyping {

  #File unpadded_intervals_file
  String callset_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbsnp_vcf
  File dbsnp_vcf_index

  Array[String] sample_names
  Array[File] input_gvcfs 
  Array[File] input_gvcfs_indices 

  String path_herramientas
  File interval_list

  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69
  
  
  #Array[String] snp_recalibration_tranche_values
  #Array[String] snp_recalibration_annotation_values
  
  #################################### atributos del indels variant recalibration
  #Array[String] indel_recalibration_tranche_values
  #Array[String] indel_recalibration_annotation_values
  #File mills_resource_vcf
  #File mills_resource_vcf_index
  #File axiomPoly_resource_vcf
  #File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf = dbsnp_vcf
  File dbsnp_resource_vcf_index = dbsnp_vcf_index
 
  
  call ImportGVCFs {
      input:
        sample_names = sample_names,
        interval = interval_list,
        #interval = unpadded_intervals[idx],
        workspace_dir_name = "mibasededatos",
        input_gvcfs = input_gvcfs,
        input_gvcfs_indices = input_gvcfs_indices,
        path_herramientas = path_herramientas,
        batch_size = 50
    }

    call GenotypeGVCFs {
      input:
        workspace_tar = ImportGVCFs.output_genomicsdb,
        interval = interval_list,
        #interval = unpadded_intervals[idx],
        output_vcf_filename = "output.vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        path_herramientas = path_herramientas
}

#call HardFilterAndMakeSitesOnlyVcf {
#      input:
#        vcf = GenotypeGVCFs.output_vcf,
#        vcf_index = GenotypeGVCFs.output_vcf_index,
#        excess_het_threshold = excess_het_threshold,
#        #variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
#        #sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
#        variant_filtered_vcf_filename = callset_name + ".variant_filtered.vcf.gz",
#        sites_only_vcf_filename = callset_name + ".sites_only.variant_filtered.vcf.gz",
#        
#        path_herramientas = path_herramientas
#}

#call GatherVcfs as SitesOnlyGatherVcf {
 #   input:
 #     input_vcfs_fofn = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
 #     input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index,
 #     output_vcf_name = callset_name + ".sites_only.vcf.gz",
 #     path_herramientas = path_herramientas
#}

#call IndelsVariantRecalibrator {
#    input:
#      sites_only_variant_filtered_vcf = GenotypeGVCFs.output_vcf,
#      sites_only_variant_filtered_vcf_index = GenotypeGVCFs.output_vcf_index,
#      recalibration_filename = callset_name + ".indels.recal",
#      tranches_filename = callset_name + ".indels.tranches",
#      recalibration_tranche_values = indel_recalibration_tranche_values,
#      recalibration_annotation_values = indel_recalibration_annotation_values,
#      mills_resource_vcf = mills_resource_vcf,
#      mills_resource_vcf_index = mills_resource_vcf_index,
#      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
#      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
#      dbsnp_resource_vcf = dbsnp_resource_vcf,
#      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
#      path_herramientas = path_herramientas
#}


}