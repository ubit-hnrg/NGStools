## Inputs: This workflow take as main input a multisample VCF (N samples) and a region.paded.bed file. 
## Output: return N single sample vcfs restricted to the provided region.paded.bed file. 
 
workflow processJointVCF {

    File multisampleVCF
    String toolpath
    File region_padded_bed ##Tso_bed
    String path_softlink
    Array[String] array_path_save


    # for annovar prouposes
    String db_annovar
    File annovar_table_pl
    File joinPY

    #File build_excell_py
    #Array[File] exon_coverage_reports

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




    scatter (sample in rename_samples.original_samples){

        call idsample{
            input:
            sample = sample
        }

        call get_individual_vcf {
            input:
            multisampleVCF = rename_samples.multisample_vcf_restricted_renamed,
            sample = idsample.idsample,
            toolpath = toolpath,
            original_sample = sample
            }

        call annovar {
            input:
            one_sample_vcf =  get_individual_vcf.one_sample_vcf,
            sample = idsample.idsample,
            annovar_table_pl = annovar_table_pl,
            db_annovar = db_annovar
         }

        call get_tsv_from_annovar {
            input:
            annovar_txt = annovar.annovar_txt,
            annovar_vcf = annovar.annovar_vcf,
            multisampleVCF = rename_samples.multisample_vcf_restricted_renamed,
            sample = idsample.idsample,
            joinPY = joinPY
        }

        #call build_excell_report{
        #    input:
        #    annovar_tsv = get_tsv_from_annovar.[annovar_tsv],
        #    exon_coverage_report =  prefix(sample,exon_coverage_reports)[0],
        #    sample=idsample.idsample,
        #    build_excell_py = build_excell_py
        #}


        call symlink_important_files2 {
            input: 
             output_to_save1 = annovar.annovar_vcf,
             output_to_save2 = get_individual_vcf.one_sample_vcf,
             path_save = path_softlink, 
             sample = sample



        }
        
    }

Array[File] salidas = ["${restrict_multisample_vcf.multisampleVCF_restricted}"]
Array[Pair[String,File]] samples_x_files = cross (array_path_save, salidas)
scatter (pairs in samples_x_files) {
    call symlink_important_files {
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
        #output_to_save = pairs,
        #path_save = array_path_save


    }
}
    
Array[File] salidas1 = ["${get_tsv_from_annovar.annovar_tsv}"]
Array[Pair[String,File]] samples_x_files1 = cross (array_path_save, salidas1)
scatter (pairs in samples_x_files1) {
    call symlink_important_files as links_a_annovar_tsv {
        input:
        output_to_save = pairs.right,
        path_save = pairs.left
        #output_to_save = pairs,
        #path_save = array_path_save


    }
}








    output{
        File multisampleVCF_restricted = restrict_multisample_vcf.multisampleVCF_restricted
        Array[File] individual_vcfs_notAnnotated = get_individual_vcf.one_sample_vcf
        Array[File] individual_vcfs_annovar = annovar.annovar_vcf
        #Array[File] individual_excell_reports = build_excell_report.excell_report

        #File multisampleVCF_restricted_renamed = rename_samples.multisample_vcf_restricted_renamed
    }

}

task symlink_important_files2 {
    File output_to_save1
    File output_to_save2
    String path_save
    String sample
    command{
       ln -s ${output_to_save1} ${path_save}${sample} 
       ln -s ${output_to_save2} ${path_save}${sample}
    }
}

task symlink_important_files {
    File output_to_save
    String path_save
 
    command{
       ln -s ${output_to_save} ${path_save}

    }
}

task idsample{
    String sample
    command<<<
        if [[ ${sample} =~ ^[0-9].* ]]
        then
            echo 'ID'${sample}
        else
            echo ${sample}
        fi
    >>>

    output{
        String idsample = read_string(stdout())
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
        #cat ${base}'_renamed.vcf' |grep '^#C'|cut -f10- |tr '\t' '\n'
        cat ${multisampleVCF_restricted} |grep '^#C'|cut -f10- |tr '\t' '\n'
    >>>

    output {
        File multisample_vcf_restricted_renamed = "${base}_renamed.vcf"
        #Array[String] renamed_samples = read_lines(stdout())
        Array[String] original_samples = read_lines(stdout())
    }

}

task get_individual_vcf{
    File multisampleVCF
    String sample
    String toolpath
    String original_sample
    

    command<<<
        ##split vcf according to $i sample. This file contain all samples but only those relevant for $i one
        cat ${multisampleVCF} | java -jar ${toolpath}/SnpSift.jar filter "(GEN[${sample}].GT!='./.')&(GEN[${sample}].GT != '0/0')" >  faceted_one_sample_vcf

        ##these steps remove the remaining samples of the vcf.
        cat <(grep '^##' faceted_one_sample_vcf) <(grep -v '^##' faceted_one_sample_vcf| csvcut -t -c '#CHROM',POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,${sample} | csvformat -T) > ${original_sample}.vcf
        rm faceted_one_sample_vcf
    >>>

    output{
    File one_sample_vcf = '${original_sample}.vcf'
    }
}

task annovar{
    File one_sample_vcf
    File annovar_table_pl
    File convert2annovar
    File annotate_variation
    File variants_reduction


    String db_annovar
    String sample
    
    command<<<
        #perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol refGene,avsnp150,esp6500siv2_all,1000g2015aug_all,exac03,gnomad_exome,gnomad_genome,clinvar_20180603,intervar_20180118,dbscsnv11,dbnsfp35a,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g -nastring . -otherinfo
        perl ${annovar_table_pl} ${one_sample_vcf} ${db_annovar} -vcfinput -buildver hg19 -remove -out ${sample} -protocol refGene,knownGene,intervar_20180118 -operation  g,g,f -nastring . -otherinfo

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
    File multisampleVCF
    String sample
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
    python ${joinPY} --multianno_tsv=${sample}.hg19_multianno.tsv --vcf_multisample=${multisampleVCF} --output=${sample}.multianno_multisample.tsv
    #change dots by tabs.
    sed -i -e "s|\.	|	|g" ${sample}.multianno_multisample.tsv

    >>>
    output{
        File annovar_tsv =  '${sample}.multianno_multisample.tsv'
    }
}

task build_excell_report{
    File annovar_tsv
    File exon_coverage_report
    String sample
    String original_sample
    File build_excell_py
    
    command{
        ${build_excell_py} ${annovar_tsv}:Variants ${exon_coverage_report}:ExonCoverage ${sample}.output_xlsx
    }     

    output{
        File excell_report = '${sample}.variants.xlsx'
    }
}    

