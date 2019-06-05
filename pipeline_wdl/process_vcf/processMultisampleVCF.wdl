## Inputs: This workflow take as main input a multisample VCF (N samples) and a region.paded.bed file. 
## Output: return N single sample vcfs restricted to the provided region.paded.bed file. 
 
workflow processJointVCF {

    File multisampleVCF
    String toolpath
    File region_padded_bed

    # for annovar prouposes
    String db_annovar
    File annovar_table_pl

#    Array[String] path_save

    call restrict_multisample_vcf{
        input:
        multisampleVCF  = multisampleVCF,
        region_padded_bed = region_padded_bed,
        toolpath = toolpath
    }


    
    call rename_samples{
        input:
        multisampleVCF_restricted  = restrict_multisample_vcf.multisampleVCF_restricted
    }




    scatter (sample in rename_samples.renamed_samples){
        
        call get_individual_vcf {
            input:
            multisampleVCF = rename_samples.multisample_vcf_restricted_renamed,
            sample = sample,
            toolpath = toolpath
            }

        call run_annovar {
            input:
            one_sample_vcf =  get_individual_vcf.one_sample_vcf,
            sample = sample,
            annovar_table_pl = annovar_table_pl,
            db_annovar = db_annovar
         }

    }
#    Array[File] individual_vcfs = get_individual_vcf.one_sample_vcf,
#    Array[File] annovar_vcfs = run_annovar.annovar_vcf,
    Array[File] annovar_txt = run_annovar.annovar_txt
    
    






    output{
        File multisampleVCF_restricted = restrict_multisample_vcf.multisampleVCF_restricted
        Array[File] individual_vcfs_notAnnotated = get_individual_vcf.one_sample_vcf
        #File multisampleVCF_restricted_renamed = rename_samples.multisample_vcf_restricted_renamed
    }

}



task restrict_multisample_vcf{
    File multisampleVCF
    File region_padded_bed
    String toolpath
    String base = basename(multisampleVCF,'.vcf.gz')
    
    command{

        zless ${multisampleVCF} | java -jar ${toolpath}/SnpSift.jar intervals ${region_padded_bed} > ${base}_restricted.vcf
    }

    output {
        File multisampleVCF_restricted = "${base}_restricted.vcf"
    }

}


task rename_samples{
    File multisampleVCF_restricted
    String base = basename(multisampleVCF_restricted,'.vcf')
    
    command<<<
        initial_sample_ids=$(cat ${multisampleVCF_restricted} |grep '^#C'|cut -f10- )
        cp ${multisampleVCF_restricted} ${base}'_renamed.vcf'
        for i in $initial_sample_ids;
            do
                # rename numeric id samples
                # prepare id  ## This is only for supportting samplenames starting with numbers. (SnpSift compatibility)
                
                if [[ $i =~ ^[0-9].* ]]
                then
                    id='ID'$i
                    sed -i -e "s/$i/ID$i/g" ${base}'_renamed.vcf';
                else
                    id=$i
                fi
            done
        cat ${base}'_renamed.vcf' |grep '^#C'|cut -f10- |tr '\t' '\n'

    >>>

    output {
        File multisample_vcf_restricted_renamed = "${base}_renamed.vcf"
        Array[String] renamed_samples = read_lines(stdout())
    }

}

task get_individual_vcf{
    File multisampleVCF
    String sample
    String toolpath
    

    command<<<
        ##split vcf according to $i sample. This file contain all samples but only those relevant for $i one
        cat ${multisampleVCF} | java -jar ${toolpath}/SnpSift.jar filter "(GEN[${sample}].GT!='./.')&(GEN[${sample}].GT != '0/0')" >  faceted_one_sample_vcf

        ##these steps remove the remaining samples of the vcf.
        cat <(grep '^##' faceted_one_sample_vcf) <(grep -v '^##' faceted_one_sample_vcf| csvcut -t -c '#CHROM',POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,${sample} | csvformat -T) > ${sample}.vcf
        rm faceted_one_sample_vcf
    >>>

    output{
    File one_sample_vcf = '${sample}.vcf'
    }
}

task run_annovar{
    File one_sample_vcf
    File annovar_table_pl
    String db_annovar
    String sample
    
    command{
        #perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol refGene,avsnp150,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_exome,gnomad_genome,clinvar_20180603,intervar_20180118,dbscsnv11,dbnsfp35a,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g -nastring . -otherinfo
        perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol refGene,knownGene -operation  g,g -nastring . -otherinfo

    }

    output {
        File annovar_vcf = '${sample}.hg19_multianno.vcf'
        File annovar_txt = '${sample}.hg19_multianno.txt'
    }

}
#scatter (fastp in fastp_json_files){
#    call fastp_qual {
#        input:
#        inputs_json_report = fastp
#    }
#}
