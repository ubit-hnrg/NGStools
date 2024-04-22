
## Requirements/expectations :
## - One or more human whole-genome per-sample GVCF files
##


workflow singleGenotypeGVCFs {
  #File unpadded_intervals_file ##lista intervalos
  File eval_interval_list
  String pipeline_version
  
  String array_path_save
  String callset_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String gatk_jar
  String toolpath

  String sample_name
  File input_gvcf 
  File input_gvcf_index

  File dbSNP_vcf
  File dbSNP_vcf_index
  File region_padded_bed
  File exon_coordinates

# for annovar prouposes
    String db_annovar
    File annovar_table_pl
    File joinPY
    File gnomad_plof_db
  


  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excess_het_threshold = 54.69


  call GenotypeGVCFs {
    input:
      #workspace_tar = ImportGVCFs.output_genomicsdb,
    gvcf = input_gvcf,
    gvcf_index = input_gvcf_index,
    #interval = unpadded_intervals[idx],
    output_vcf_filename = "output.vcf.gz",
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict,
    dbSNP_vcf = dbSNP_vcf,
    dbSNP_vcf_index = dbSNP_vcf_index,
    gatk_jar = gatk_jar,
    toolpath = toolpath
    
    }


  call HardFilterAndMakeSitesOnlyVcf {
    input:
    vcf = GenotypeGVCFs.output_vcf,
    vcf_index = GenotypeGVCFs.output_vcf_index,
    excess_het_threshold = excess_het_threshold,
    variant_filtered_vcf_filename = callset_name + ".single.variant_filtered.vcf.gz",
    sites_only_vcf_filename = callset_name + ".single.sites_only.variant_filtered.vcf.gz",
    interval_list = region_padded_bed,
    gatk_jar = gatk_jar,
    toolpath = toolpath
   
  }

  call GatherVcfs as FinalGatherVcf {
    input:
    input_vcfs_fofn  = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
    input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
    output_vcf_name = callset_name+"."+pipeline_version+".vcf.gz",
    gatk_jar = gatk_jar,
    toolpath = toolpath

  }

  call CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
    input:
    input_vcf = FinalGatherVcf.output_vcf,
    input_vcf_index = FinalGatherVcf.output_vcf_index,
    metrics_filename_prefix = callset_name,
    dbSNP_vcf = dbSNP_vcf,
    dbSNP_vcf_index = dbSNP_vcf_index,
    interval_list = eval_interval_list,
    ref_dict = ref_dict,
    gatk_jar = gatk_jar,
    toolpath = toolpath

  }
  

  call restrict_vcf {
    input:
    VCF  = FinalGatherVcf.output_vcf,
    region_padded_bed = region_padded_bed,
    toolpath = toolpath
  }

call exon_distance {
    input:
    vcf_ok = restrict_vcf.VCF_restricted, #step12_clinVar.salida_Snpsift, #final_annot.salida_Snpsift,
    exon_coord = exon_coordinates,
    #exon_coordinates_to_lib = exon_coordinates_to_lib,
    sample_name = sample_name
}
#   call filtro_no_calls {
#     input:
#     gatk_jar = gatk_jar,
#     vcf_in = restrict_vcf.VCF_restricted,
#     sample_name = sample_name,
#     toolpath = toolpath,
#     ref_fasta = ref_fasta,
#     ref_fasta_index = ref_fasta_index,
#     ref_dict = ref_dict,
#     path_save = array_path_save
# }



  call annovar {
            input:
            one_sample_vcf =  restrict_vcf.VCF_restricted, #filtro_no_calls.one_sample_vcf,,#get_individual_vcf.one_sample_vcf,
            sample = sample_name,#idsample.idsample,
            annovar_table_pl = annovar_table_pl,
            toolpath = toolpath,

            db_annovar = db_annovar
         }



        call get_tsv_from_annovar {
            input:
            gnomad_plof = gnomad_plof_db,
            annovar_txt = annovar.annovar_txt,
            annovar_vcf = annovar.annovar_vcf,
            restrictedVCF = restrict_vcf.VCF_restricted,#rename_samples.multisample_vcf_restricted_renamed,
            sample = sample_name,#idsample.idsample,
            #sample1 = sample,
            joinPY = joinPY
        }

  #Array[File] salidas = ["${FinalGatherVcf.output_vcf}","${FinalGatherVcf.output_vcf_index}","${CollectMetricsOnFullVcf.detail_metrics_file}","${CollectMetricsOnFullVcf.summary_metrics_file}"]
  #Array[Pair[String,File]] samples_x_files = cross (array_path_save, salidas)
  #scatter (pairs in samples_x_files) {
  call symlink_important_files {
    input:
    final_gath = FinalGatherVcf.output_vcf,
    final_gath_idx = FinalGatherVcf.output_vcf_index,
    restricted_vcf = restrict_vcf.VCF_restricted,

    metrica1 = CollectMetricsOnFullVcf.detail_metrics_file,
    metrica2 = CollectMetricsOnFullVcf.summary_metrics_file,
    path_save = array_path_save
  }
call symlink_important_files2 {
            input: 
             output_to_save1 = annovar.annovar_vcf,
             #output_to_save2 = get_individual_vcf.one_sample_vcf,
             output_to_save3 = get_tsv_from_annovar.annovar_tsv,
             path_save = array_path_save, 
             sample = sample_name
        }

  output {
    # outputs from the small callset path through the wdl
  
   File vcf_exon_distance = exon_distance.exon_dist 



   File? restricted_vcf = restrict_vcf.VCF_restricted
   File? outputvcf = FinalGatherVcf.output_vcf
   File? outputvcfindex =  FinalGatherVcf.output_vcf_index
   File? metrica1 =  CollectMetricsOnFullVcf.detail_metrics_file
   File? metrica2 = CollectMetricsOnFullVcf.summary_metrics_file
   #File individual_vcfs_notAnnotated = get_individual_vcf.one_sample_vcf
   File individual_vcfs_annovar = annovar.annovar_vcf
   #Array[File] individual_excell_reports = build_excell_report.excell_report
   File annovar_tsv_out = get_tsv_from_annovar.annovar_tsv
   File annovar_gene_list = get_tsv_from_annovar.gene_list
   File gene_plof_file = get_tsv_from_annovar.gene_plof
    # outputs from the large callset path through the wdl
    # (note that we do not list ApplyRecalibration here because it is run in both paths)
    #GatherMetrics.detail_metrics_file
    #GatherMetrics.summary_metrics_file

    # output the interval list generated/used by this run workflow
    #File? inter = DynamicallyCombineIntervals.output_intervals
  }
}

task symlink_important_files {
    File final_gath
    File final_gath_idx
    File restricted_vcf
    
    File metrica1
    File metrica2
    String path_save
    command{
       cp -L ${final_gath} ${path_save}
       cp -L ${final_gath_idx} ${path_save}
       cp -L ${metrica1} ${path_save}
       cp -L ${metrica2} ${path_save}
       cp -L ${restricted_vcf} ${path_save}
    
    }
}

task GetNumberOfSamples {
  File sample_name_map

  command <<<
    wc -l ${sample_name_map} | awk '{print $1}'
  >>>

  output {
    Int sample_count = read_int(stdout())
  }
}

task ImportGVCFs {
  String sample_name
  File input_gvcf
  File input_gvcf_index


  String interval
  String gatk_jar
  String toolpath

  String workspace_dir_name


  Int batch_size

  command <<<
    set -e
    set -o pipefail
    
    python3 << CODE
    gvcf = ['${sep="','" input_gvcf}']
    sample_name = ['${sep="','" sample_name}']

    if len(gvcf)!= len(sample_name):
      exit(1)

    with open("inputs.list", "w") as fi:
      for i in range(len(gvcf)):
       fi.write(sample_name[i] + "\t" + gvcf[i] + "\n") 
    
    CODE
    
    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    java -Xmx4g -Xms4g -jar ${toolpath}${gatk_jar} \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${interval} \
    --sample-name-map inputs.list  \
    --reader-threads 4 \
    -ip 500

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}
  
  >>>

  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {
  #File workspace_tar
  #String interval

  String output_vcf_filename

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index    
 
  String gatk_jar
  String toolpath
  File gvcf
  File gvcf_index

  command <<<
    set -e


    java -Xmx4g -Xms4g -jar ${toolpath}${gatk_jar}\
     GenotypeGVCFs \
     -R ${ref_fasta} \
     -O ${output_vcf_filename} \
     -D ${dbSNP_vcf} \
     -G StandardAnnotation \
     --use-new-qual-calculator \
     --variant ${gvcf}
     
  >>>
#-L ${interval}
  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {
  File vcf
  File vcf_index
  Float? excess_het_threshold
  File interval_list
  String variant_filtered_vcf_filename
  String sites_only_vcf_filename

  String gatk_jar
  String toolpath

###### agrego 
###       --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
#   --filter-name "HardFiltered" \
# -L 
##quito
# --filter-expression "ExcessHet > ${excess_het_threshold}" \
     # --filter-name ExcessHet \
#
#     
  command {
    set -e

    java -Xmx3g -Xms3g -jar ${toolpath}${gatk_jar}\
      VariantFiltration \
     --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
     --filter-name "HardFiltered" \
     -L ${interval_list} \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    java -Xmx3g -Xms3g -jar ${toolpath}${gatk_jar} \
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



task GatherVcfs {
  File input_vcfs_fofn
  File input_vcf_indexes_fofn
  String output_vcf_name
  
  String gatk_jar
  String toolpath


    command <<<
    set -e
    set -o pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    java -Xmx6g -Xms6g -jar ${toolpath}${gatk_jar} \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    -I ${sep=" -I " input_vcfs_fofn} \
    --output ${output_vcf_name}

    java -Xmx6g -Xms6g -jar ${toolpath}${gatk_jar}\
    IndexFeatureFile \
  --input ${output_vcf_name}
  >>>
 
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {
  File input_vcf
  File input_vcf_index

  String metrics_filename_prefix
  File dbSNP_vcf
  File dbSNP_vcf_index
  File interval_list
  File ref_dict

  String gatk_jar
  String toolpath



  command {
    ##### --SEQUENCE_DICTIONARY ${ref_dict}
    java -Xmx4g -Xms4g -jar ${toolpath}${gatk_jar} \
      CollectVariantCallingMetrics \
      --INPUT ${input_vcf} \
      --DBSNP ${dbSNP_vcf} \
      --OUTPUT ${metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${interval_list}
  }
  output {
    File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  }

}

task GatherMetrics {
  File input_details_fofn
  File input_summaries_fofn

  String output_prefix

  String gatk_jar
  String toolpath
  

  command <<<
    set -e
    set -o pipefail

    
    java -Xmx4g -Xms4g -jar ${toolpath}${gatk_jar} \
    AccumulateVariantCallingMetrics \
    --INPUT ${sep=" --INPUT " input_details_fofn} \
    --OUTPUT ${output_prefix}
  >>>

  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}

task DynamicallyCombineIntervals {
  File intervals
  Int merge_count


  command {
    python3 << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ${merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("${intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
  }


  output {
    File output_intervals = "out.intervals"
  }
}

task restrict_vcf{
    File VCF
    File region_padded_bed
    String toolpath
    String base = basename(VCF,'.vcf.gz')
    
    command{

        zless ${VCF} | java -jar ${toolpath}/SnpSift.jar intervals ${region_padded_bed} > ${base}_restricted.vcf
    }

    output {
        File VCF_restricted = "${base}_restricted.vcf"
    }

}

task filtro_no_calls {


File ref_fasta
File ref_fasta_index
File ref_dict

File vcf_in
String sample_name
String toolpath
String gatk_jar

String path_save


    command<<<
        ##split vcf according to $i sample. This file contain all samples but only those relevant for $i one
        cat ${vcf_in} | java -jar ${toolpath}/SnpSift.jar filter "(GEN[${sample_name}].GT!='./.')&(GEN[${sample_name}].GT != '0/0')" >  faceted_one_sample_vcf

        ##these steps remove the remaining samples of the vcf.
        cat <(grep '^##' faceted_one_sample_vcf) <(grep -v '^##' faceted_one_sample_vcf| csvcut -t -c '#CHROM',POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,${sample_name} | csvformat -T) > ${sample_name}_filter.vcf 

        rm faceted_one_sample_vcf

      #   java -Xmx2g -Xms2g -jar ${toolpath}${gatk_jar} \
      #   SelectVariants \
      #   -V ${vcf_in} \
      #   -R ${ref_fasta}\
      #   -output ${sample_name}_filter_agu.vcf \
      #   -invert-select 'vc.getHomRefCount() == vc.getCalledChrCount()/2' \ # outcome: vcf sites where all the samples are not GT = ./. (no call) or GT = 0/0

      #  cp -L ${sample_name}_filter.vcf ${path_save}
      #  cp -L ${sample_name}_filter_agu.vcf ${path_save}

    >>>

output {
    File one_sample_vcf = '${sample_name}_filter.vcf'
   # File filtro_nocall_vcf =  '${sample_name}_filter_agu.vcf'
    }  

}

#### annovar

task annovar{
    File one_sample_vcf
    File annovar_table_pl
    String toolpath

    #File convert2annovar = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/2024/annovar/convert2annovar.pl'
    #File annotate_variation = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/2024/annovar/annotate_variation.pl'
    #File variants_reduction = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/2024/annovar/variants_reduction.pl'
    #File coding_change = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/2024/annovar/coding_change.pl'
    #File retrieve_seq_from_fasta = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/2024/annovar/retrieve_seq_from_fasta.pl'


    String db_annovar #= '/data/new_dbs/annovar/hg38/humandb/'
    String sample 
    
## refGene,intervar_20180118,esp6500siv2_all,1000g2015aug_all,exac03,gnomad312_exome,gnomad312_genome,clinvar_20221231,dbscsnv11,dbnsfp42a,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g
#--slicing_threshold 10bp away from splicesite.
    #perl ${toolpath}annovar/2024/annovar/table_annovar.pl ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg38 -thread 4 -remove -out ${sample} -protocol refGene,intervar_20180118,esp6500siv2_all,1000g2015aug_all,exac03,gnomad40_exome,gnomad40_genome,clinvar_20221231,dbscsnv11,rmsk,cytoBand,wgRna,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation g,f,f,f,f,f,f,f,f,r,r,r,r,r,r,g,g -nastring . --otherinfo -polish -intronhgvs 30

    command<<<
        perl ${toolpath}annovar/table_annovar.pl ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg38 -thread 4 -remove -out ${sample} -protocol refGene,intervar_20180118,ensGene,knownGene -operation g,f,g,g -nastring . --otherinfo -polish -intronhgvs 30
        #refGene,intervar_20180118,esp6500siv2_all,1000g2015aug_all,exac03,gnomad40_exome,gnomad40_genome,clinvar_20221231,dbscsnv11,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g -nastring . -otherinfo --slicing_threshold 30 -polish -intronhgvs

    ###dbnsfp con q anoto? con snpsift?
    >>>

    output {
        File annovar_vcf = '${sample}.hg38_multianno.vcf'
        File annovar_txt = '${sample}.hg38_multianno.txt'
    }

}

 # prepare final multianno.tsv for deliver (do not forget postprocess InterVar)
 #1) localizo las 3 columnas de Otherinfo para volarlas. 
 #2) modifico el header agregando las columnas que faltan del vcf (remplazo x otherinfo)
 #3) joineo al archivo multinanno las columnas de genotipo de las restantes muestras de  la corrida.

task get_tsv_from_annovar {
    File annovar_txt
    File annovar_vcf
    File restrictedVCF
    String sample
    #String sample1
    File joinPY    #this file merge the multianno.tsv file with the original multisample vcf
    File gnomad_plof ###gnomad plof for hnrg -> lo usan en brasil.

    command <<<
    
    #columnas a cortar (localizando Otherinfo column y las 2 siguientes)
    #nl0=$(head -n1 ${annovar_txt}|tr '\t' '\n'|nl|grep 'Otherinfo1'|cut -f1)
    #nl1=$((nl0 + 1))
    #nl2=$((nl0 + 2))
    ###modifico para hacer dinamico el rango de otherinfo cols enero2024

    #header=$(head -n1 ${annovar_txt})

# Utiliza awk para obtener los números de las columnas que contienen "Otherinfo"
#otherinfo_cols=$(echo "$header" | awk '{
#    for(i = 1; i <= NF; i++) {
#        if($i ~ /Otherinfo/) {
#            printf "%s,", i;
#        }
#    }
#}' | sed 's/,$//')



    # meto header (dejando el campo 'Otherinfo' que despues va a ser remplazado por las columnas del vcf original)
    #head -n1 ${annovar_txt} > ${sample}.hg38_multianno.tsv
    # vuelo las tres columnas de otherinfo
    #tail -n+2 ${annovar_txt}|cut -f$nl0,$nl1,$nl2 --complement >>  ${sample}.hg38_multianno.tsv; #modifico para usar el rango dinamco
    #tail -n+2 ${annovar_txt}|cut -f$otherinfo_cols --complement >>  ${sample}.hg38_multianno.tsv;

    #####agu enero 2024
    header=$(head -n1 ${annovar_txt})
    otherinfo_cols=$(echo "$header" | awk '{
    for(i = 1; i <= NF; i++) {
        if($i == "Otherinfo1" || $i == "Otherinfo2" || $i == "Otherinfo3") {
            printf "%s,", i;
        }
    }
}' | sed 's/,$//')

# Convierte los números de columna en un rango para 'cut'
cut_range=$(echo $otherinfo_cols | sed 's/,/,/g')

head -n1 ${annovar_txt} | cut -f$cut_range --complement >  ${sample}.hg38_multianno.tsv
tail -n+2 ${annovar_txt} | cut -f$cut_range --complement >>  ${sample}.hg38_multianno.tsv

vcf_header=$(grep '#CH' ${annovar_vcf});
otherinfo_index=4
vcf_col_index=1

total_vcf_fields=$(echo $vcf_header | awk '{print NF}')

# Bucle para reemplazar Otherinfo4, Otherinfo5, ..., Otherinfo10
for (( ; otherinfo_index<=13 && vcf_col_index<=total_vcf_fields; otherinfo_index++, vcf_col_index++ )); do
    # Extrae el elemento actual del encabezado VCF
    vcf_field=$(echo $vcf_header | awk -v col=$vcf_col_index '{print $col}')

    # Reemplaza "OtherinfoX" con el campo correspondiente del VCF
    sed -i "s/Otherinfo$otherinfo_index/$vcf_field/" ${sample}.hg38_multianno.tsv
done
#####fin agu 

    #remplazo el header
    #sed -i "s/Otherinfo/$vcf_header/g" ${sample}.hg38_multianno.tsv;

    #join one multianno tsv file AND joint genotyped vcf. This script (join_vcf.py) also postprocess Intervar columns.
    python3 ${joinPY} --multianno_tsv=${sample}.hg38_multianno.tsv --vcf_multisample=${restrictedVCF} --output=${sample}.multianno_restrict.tsv
    #change dots by tabs.
    sed -i -e "s|\.	|	|g" ${sample}.multianno_restrict.tsv
    
    ####agrego un awk para buscar los genes de annovar y hacer un archivo con la tabla gnomad_plof para esos genes.
    cat ${sample}.multianno_restrict.tsv | cut -f20 | uniq > ${sample}.gene_list_for_plof.list
    head -n1 ${gnomad_plof} > ${sample}_plof.tsv
    awk 'NR == FNR {gene_list[$1];next} ($1 in gene_list)' ${sample}.gene_list_for_plof.list ${gnomad_plof} >> ${sample}_plof.tsv
 
    >>>
    output{
        File annovar_tsv =  '${sample}.multianno_restrict.tsv'
        File gene_list = '${sample}.gene_list_for_plof.list'
        File gene_plof = '${sample}_plof.tsv'
    }
}

task symlink_important_files2 {
    File output_to_save1
    #File output_to_save2
    File output_to_save3
    String path_save
    String sample
    command{
       cp -L ${output_to_save1} ${path_save} #${sample} 
       cp -L ${output_to_save3} ${path_save} #${sample}
    }
}

task exon_distance {
File vcf_ok
File exon_coord
String sample_name
#File exon_coordinates_to_lib
#String path_save

command{
    #!/bin/bash
    set -e

    grep "^#" ${vcf_ok} > head_vcf.vcf
    grep -v "#" ${vcf_ok}| sort -k1,1 -k2,2n >> head_vcf.vcf
    
    ##test
    ##sort -k1,1V -k2,2n ${exon_coord} >> sorted_exon_bed.bed
    
    bedtools closest -a head_vcf.vcf -b ${exon_coord} -D a | cut -f1,2,12-18  > ${sample_name}.exon_distance.tsv

    rm head_vcf.vcf
    }

output {
File exon_dist = "${sample_name}.exon_distance.tsv"
#File exon_dist_to_lib = "${sample_name}.exon_distance_tolib.tsv"

}
}