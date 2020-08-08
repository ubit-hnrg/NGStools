
## Requirements/expectations :
## - One or more human whole-genome per-sample GVCF files
##


workflow singleGenotypeGVCFs {
  #File unpadded_intervals_file ##lista intervalos
  File eval_interval_list
  String pipeline_v
  
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

# for annovar prouposes
    String db_annovar
    File annovar_table_pl
    File joinPY
  


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
    gatk_jar = gatk_jar,
    toolpath = toolpath
   
  }

  call GatherVcfs as FinalGatherVcf {
    input:
    input_vcfs_fofn  = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf,
    input_vcf_indexes_fofn = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index,
    output_vcf_name = callset_name+"."+pipeline_v+".vcf.gz",
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
  

  call restrict_vcf{
    input:
    VCF  = FinalGatherVcf.output_vcf,
    region_padded_bed = region_padded_bed,
    toolpath = toolpath
  }

  call annovar {
            input:
            one_sample_vcf =  restrict_vcf.VCF_restricted,#get_individual_vcf.one_sample_vcf,
            sample = sample_name,#idsample.idsample,
            annovar_table_pl = annovar_table_pl,
            db_annovar = db_annovar
         }



        call get_tsv_from_annovar {
            input:
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
   File? restricted_vcf = restrict_vcf.VCF_restricted
   File? outputvcf = FinalGatherVcf.output_vcf
   File? outputvcfindex =  FinalGatherVcf.output_vcf_index
   File? metrica1 =  CollectMetricsOnFullVcf.detail_metrics_file
   File? metrica2 = CollectMetricsOnFullVcf.summary_metrics_file
   #File individual_vcfs_notAnnotated = get_individual_vcf.one_sample_vcf
   File individual_vcfs_annovar = annovar.annovar_vcf
   #Array[File] individual_excell_reports = build_excell_report.excell_report
   File annovar_tsv_out = get_tsv_from_annovar.annovar_tsv

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
    
    python << CODE
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
    java -Xmx1g -Xms1g -jar ${toolpath}${gatk_jar} \
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


    java -Xmx1g -Xms1g -jar ${toolpath}${gatk_jar}\
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
  Float excess_het_threshold

  String variant_filtered_vcf_filename
  String sites_only_vcf_filename

  String gatk_jar
  String toolpath

  command {
    set -e

    java -Xmx3g -Xms3g -jar ${toolpath}${gatk_jar}\
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}

    java -Xmx3g -Xms3g -jar ${toolpath}${gatk_jar} \
      MakeSitesOnlyVcf \
     -I=${variant_filtered_vcf_filename} \
     -O=${sites_only_vcf_filename}

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
  --feature-file ${output_vcf_name}
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

    
    java -Xmx2g -Xms2g -jar ${toolpath}${gatk_jar} \
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
    python << CODE
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

#### annovar

task annovar{
    File one_sample_vcf
    File annovar_table_pl
    File convert2annovar = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/convert2annovar.pl'
    File annotate_variation = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/annotate_variation.pl'
    File variants_reduction = '/home/hnrg/HNRG-pipeline-V0.1/tools/annovar/variants_reduction.pl'


    String db_annovar
    String sample
    
    command<<<
        perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol refGene,avsnp150,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_exome,gnomad_genome,clinvar_20180603,intervar_20180118,dbscsnv11,dbnsfp35a,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g -nastring . -otherinfo
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
    File restrictedVCF
    String sample
    #String sample1
    File joinPY    #this file merge the multianno.tsv file with the original multisample vcf


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
    python ${joinPY} --multianno_tsv=${sample}.hg19_multianno.tsv --vcf_multisample=${restrictedVCF} --output=${sample}.multianno_restrict.tsv
    #change dots by tabs.
    sed -i -e "s|\.	|	|g" ${sample}.multianno_restrict.tsv

    >>>
    output{
        File annovar_tsv =  '${sample}.multianno_restrict.tsv'
    }
}

task symlink_important_files2 {
    File output_to_save1
    #File output_to_save2
    File output_to_save3
    String path_save
    String sample
    command{
       cp -L ${output_to_save1} ${path_save}${sample} 
       cp -L ${output_to_save3} ${path_save}${sample}
    }
}