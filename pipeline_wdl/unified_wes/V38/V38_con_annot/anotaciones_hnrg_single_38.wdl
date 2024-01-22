############################
### pipeline anotacion funcional para hg38
#################



##: bptools 
task bptools {

File input_vcf
String samplename1
String toolpath
String java_heap_memory_initial
String parametros
String nombre_step




command {

set -o pipefail


java -Xmx${java_heap_memory_initial} -jar ${toolpath}bptools.jar ${parametros} ${input_vcf} ${samplename1}.${nombre_step}.vcf
}


output {
    
    File bptools_out = "${samplename1}.${nombre_step}.vcf"
    }

}

task Snpeff {
File input_vcf
String samplename1
String toolpath
String java_heap_memory_initial
String reference_version = "GRCh38.mane.1.2.refseq"

###la version 5 de snpeff no soporta el -t, multihread

#java -Xmx${java_heap_memory_initial} -jar ${toolpath}SnpEff/snpEff/snpEff.jar ${reference_version} -hgvs \
# -lof -noStats -canon -onlyProtein -c ${toolpath}SnpEff/snpEff/snpEff.config \
#${input_vcf} > ${samplename1}.step1_SnpEff.vcf
command {

set -o pipefail
java -Xmx${java_heap_memory_initial} -jar ${toolpath}SnpEff/snpEff/snpEff.jar ${reference_version} -hgvs \
 -lof -noStats -geneID -canon  -c ${toolpath}SnpEff/snpEff/snpEff
 .config \
${input_vcf} > ${samplename1}.step1_SnpEff.vcf 
}

output {
    File step1_snpeff = "${samplename1}.step1_SnpEff.vcf"

}

}



task Snpsift {

File input_vcf
String samplename1
String toolpath
String java_heap_memory_initial
String nombre_step
File database
File database_index
String parametros

command {
set -o pipefail
java -Xmx${java_heap_memory_initial} -jar ${toolpath}SnpEff/snpEff/SnpSift.jar ${parametros} ${database} \
${input_vcf} > ${samplename1}.${nombre_step}.vcf

}

output {
    File salida_Snpsift = "${samplename1}.${nombre_step}.vcf"
}

}


task annovar {

File vcf_in
String toolpath
String out_prefix
String annovar_dbpath 
#dbpath=/home/hnrg/HNRG-pipeline-V0.1/dbs/hg19_annovar/
#annovar=/home/bitgenia/dbs/annovar/table_annovar.pl


command <<<

#perl ${toolpath}annovar/table_annovar.pl ${vcf_in} ${annovar_dbpath} -vcfinput  -buildver hg19 -remove -out ${out_prefix} -protocol dbnsfp35a,gnomad_exome,gnomad_genome,intervar_20180118 -operation f,f,f,f -nastring . 
#perl ${toolpath}annovar/table_annovar.pl ${vcf_in} ${annovar_dbpath} --thread 4 -vcfinput  -buildver hg19 -remove -out ${out_prefix} -protocol dbnsfp35a,gnomad_exome,gnomad_genome,intervar_20180118 -operation f,f,f,f 
perl ${toolpath}annovar/table_annovar.pl ${vcf_in} ${annovar_dbpath} -vcfinput -buildver hg38 -thread 4 -remove -out ${out_prefix} -protocol gnomad40_exome,gnomad40_genome,intervar_20180118 -operation f,f,f -nastring . -slicing_threshold 10 , -polish -intronhgvs
#refGene,intervar_20180118,esp6500siv2_all,1000g2015aug_all,exac03,gnomad40_exome,gnomad40_genome,clinvar_20221231,dbscsnv11,rmsk,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,ensGene,knownGene -operation  g,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,g,g -nastring . -otherinfo --slicing_threshold 10 , -polish -intronhgvs
>>>

output {

 File avinput = "${out_prefix}.avinput"
 File multianno_txt = "${out_prefix}.hg38_multianno.txt"
 File multianno_vcf = "${out_prefix}.hg38_multianno.vcf"
}

}

task intervar_postprocessing {
    File vcfinput
    File multianno_txt
    String srcpath = '/home/hnrg/NGStools/pipeline_wdl/anotacion_funcional'
    String toolpath
    String samplename1
    command<<<
    vcfDB='intervar_sample_BD.vcf'

    # multianno to vcf db file
    python ${srcpath}/create_InterVarDB.py -i=${multianno_txt} -o $vcfDB

    # index
    bgzip $vcfDB
    tabix $vcfDB.gz

    sed -e "s|__vcfDB__|$vcfDB.gz|g" ${srcpath}/intervar_vcfanno_template_fromVCF.tom > config_vcfanno.tom
    ${toolpath}vcfanno_linux64 -p 4 config_vcfanno.tom ${vcfinput} > ${samplename1}_intervar.vcf

    >>>
    output {
        File salida_intervar = "${samplename1}_intervar.vcf"
    }


}

task Snpsift_GWASCat {

File input_vcf
String samplename1
String toolpath
String java_heap_memory_initial
String nombre_step
String parametros
File database

command {
set -o pipefail
java -Xmx${java_heap_memory_initial} -jar ${toolpath}SnpEff/snpEff/SnpSift.jar ${parametros} ${database} \
${input_vcf} > ${samplename1}.${nombre_step}.vcf

}

output {
    File salida_Snpsift = "${samplename1}.${nombre_step}.vcf"
}

}

task Snpsift_nodb {

File input_vcf
String samplename1
String toolpath
String java_heap_memory_initial
String nombre_step
String parametros

command {
set -o pipefail
java -Xmx${java_heap_memory_initial} -jar ${toolpath}SnpEff/snpEff/SnpSift.jar ${parametros} \
${input_vcf} > ${samplename1}.${nombre_step}.vcf

}

output {
    File salida_Snpsift = "${samplename1}.${nombre_step}.vcf"
}

}


task symlink_important_files {
    File output_to_save
    File output_to_save2
    String path_save
    command{
       cp -L ${output_to_save} ${path_save}
       cp -L ${output_to_save2} ${path_save}
    }
}


task hnrg_freq {

    File input_vcf
    File config_file_vcfanno
    String samplename1
    String toolpath
    String nombre_step


    

    command {
        ${toolpath}vcfanno_linux64 -p 4 ${config_file_vcfanno} ${input_vcf} > ${samplename1}.${nombre_step}.vcf

    }
output {
    File out_vcfanno = "${samplename1}.${nombre_step}.vcf"
}

}

task exon_distance {
File vcf_ok
File exon_coord
String sample_name
#File exon_coordinates_to_lib
String path_save

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


workflow FuncionalAnnotationSingle {

File input_vcf 
##String path_herramientas#### cambiaar por toolpath
String toolpath 
String samplename1
String java_heap_memory_initial
String reference_version
String path_save
String pipeline_version
File exon_coordinates
#File exon_coordinates_to_lib









call bptools as step_0_bptools_mma {
    input: 
    samplename1 = samplename1,
    parametros = "-mma",
    input_vcf = input_vcf,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step0_splitMAA"


}

call Snpeff as step_1_Snpeff {
input:
    samplename1 = samplename1,
    input_vcf = step_0_bptools_mma.bptools_out,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    reference_version = reference_version

}

call bptools as step_2_bptools_variant_annotation {
 input: 
    samplename1 = samplename1,
    parametros = "-bpann",
    input_vcf = step_1_Snpeff.step1_snpeff,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "Step2_VariantAnnotator"


}

#Step 3: Annotate with dbSNP156"
call Snpsift as step3_dbSNP {
input:
    samplename1 = samplename1,
    parametros = "annotate",
    input_vcf = step_2_bptools_variant_annotation.bptools_out,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step3_dbSNP"
    
}

#Step 4: Annotate with 1000Genomes
call Snpsift as step4_1000Genomes {
input:
    samplename1 = samplename1,
    parametros = "annotate",
    input_vcf = step3_dbSNP.salida_Snpsift,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step4_1000Genomes"
}

#Step 5: Annotate with hapmap
call Snpsift as step5_hapmap{
input:
    samplename1 = samplename1,
    parametros = "annotate -v",
    #input_vcf = step3_dbSNP151.salida_Snpsift,
    input_vcf = step4_1000Genomes.salida_Snpsift,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step5_hapmap"

}

#Step 6: Annotate with GWASCat (with extra fields)
call Snpsift_GWASCat as step6_Snpsift_GWASCat{
input:
    samplename1 = samplename1,
    parametros = "gwasCat -db",
    input_vcf = step5_hapmap.salida_Snpsift,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step6_Snpsift_GWASCat"

}

#Step 7: Annotate with dbNSFP
call Snpsift as step7_dbNSFP{
input:
    samplename1 = samplename1,
    #parametros = 'dbnsfp -v -f "hg18_pos(1-coor),hg38_chr,hg38_pos,genename,Uniprot_acc,Uniprot_id,cds_strand,CADD_raw,CADD_phred,CADD_raw_rankscore,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_rankscore,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,VEST3_score,VEST3_rankscore,GERP++_RS,GERP++_RS_rankscore,COSMIC_ID,COSMIC_CNT,phyloP46way_placental,phyloP100way_vertebrate,phastCons46way_placental,phastCons100way_vertebrate,Ensembl_geneid,Ensembl_transcriptid"  -db',
    #parametros = 'DbNsfp -v -f "pos(1-based),ref,alt,aaref,aaalt,rs_dbSNP150,hg19_chr,hg19_pos(1-based),hg18_chr,hg18_pos(1-based),genename,cds_strand,refcodon,codonpos,codon_degeneracy,Ancestral_allele,AltaiNeandertal,Denisova,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,aapos,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Uniprot_acc_Polyphen2,Uniprot_id_Polyphen2,Uniprot_aapos_Polyphen2,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_UniprotID,MutationAssessor_variant,MutationAssessor_score,MutationAssessor_score_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,Transcript_id_VEST3,Transcript_var_VEST3,VEST3_score,VEST3_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,CADD_raw,CADD_raw_rankscore,CADD_phred,DANN_score,DANN_rankscore,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,Eigen_coding_or_noncoding,Eigen-raw,Eigen-phred,Eigen-PC-raw,Eigen-PC-phred,Eigen-PC-raw_rankscore,GenoCanyon_score,GenoCanyon_score_rankscore,integrated_fitCons_score,integrated_fitCons_score_rankscore,integrated_confidence_value,GM12878_fitCons_score,GM12878_fitCons_score_rankscore,GM12878_confidence_value,H1-hESC_fitCons_score,H1-hESC_fitCons_score_rankscore,H1-hESC_confidence_value,HUVEC_fitCons_score,HUVEC_fitCons_score_rankscore,HUVEC_confidence_value,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP20way_mammalian,phyloP20way_mammalian_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,TWINSUK_AC,TWINSUK_AF,ALSPAC_AC,ALSPAC_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,ExAC_nonTCGA_AC,ExAC_nonTCGA_AF,ExAC_nonTCGA_Adj_AC,ExAC_nonTCGA_Adj_AF,ExAC_nonTCGA_AFR_AC,ExAC_nonTCGA_AFR_AF,ExAC_nonTCGA_AMR_AC,ExAC_nonTCGA_AMR_AF,ExAC_nonTCGA_EAS_AC,ExAC_nonTCGA_EAS_AF,ExAC_nonTCGA_FIN_AC,ExAC_nonTCGA_FIN_AF,ExAC_nonTCGA_NFE_AC,ExAC_nonTCGA_NFE_AF,ExAC_nonTCGA_SAS_AC,ExAC_nonTCGA_SAS_AF,ExAC_nonpsych_AC,ExAC_nonpsych_AF,ExAC_nonpsych_Adj_AC,ExAC_nonpsych_Adj_AF,ExAC_nonpsych_AFR_AC,ExAC_nonpsych_AFR_AF,ExAC_nonpsych_AMR_AC,ExAC_nonpsych_AMR_AF,ExAC_nonpsych_EAS_AC,ExAC_nonpsych_EAS_AF,ExAC_nonpsych_FIN_AC,ExAC_nonpsych_FIN_AF,ExAC_nonpsych_NFE_AC,ExAC_nonpsych_NFE_AF,ExAC_nonpsych_SAS_AC,ExAC_nonpsych_SAS_AF,gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AN,gnomAD_exomes_AFR_AF,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AN,gnomAD_exomes_AMR_AF,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AN,gnomAD_exomes_ASJ_AF,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AN,gnomAD_exomes_FIN_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AN,gnomAD_exomes_NFE_AF,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AN,gnomAD_exomes_SAS_AF,gnomAD_exomes_OTH_AC,gnomAD_exomes_OTH_AN,gnomAD_exomes_OTH_AF,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_AFR_AC,gnomAD_genomes_AFR_AN,gnomAD_genomes_AFR_AF,gnomAD_genomes_AMR_AC,gnomAD_genomes_AMR_AN,gnomAD_genomes_AMR_AF,gnomAD_genomes_ASJ_AC,gnomAD_genomes_ASJ_AN,gnomAD_genomes_ASJ_AF,gnomAD_genomes_EAS_AC,gnomAD_genomes_EAS_AN,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AC,gnomAD_genomes_FIN_AN,gnomAD_genomes_FIN_AF,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AN,gnomAD_genomes_NFE_AF,gnomAD_genomes_OTH_AC,gnomAD_genomes_OTH_AN,gnomAD_genomes_OTH_AF,clinvar_rs,clinvar_clnsig,clinvar_trait,clinvar_golden_stars,Interpro_domain,GTEx_V6p_gene,GTEx_V6p_tissue" -db',
    #parametros = 'DbNsfp -v -f "MPC_score,BayesDel_addAF_pred,BayesDel_addAF_score,BayesDel_noAF_pred,BayesDel_noAF_score,CADD_phred,CADD_raw,CADD_raw_rankscore,DANN_rankscore,DANN_score,Ensembl_geneid,Ensembl_proteinid,Ensembl_transcriptid,FATHMM_converted_rankscore,FATHMM_pred,FATHMM_score,GENCODE_basic,GERP++_RS,GERP++_RS_rankscore,HGVSc_ANNOVAR,HGVSc_VEP,HGVSc_snpEff,HGVSp_ANNOVAR,HGVSp_VEP,HGVSp_snpEff,LRT_converted_rankscore,LRT_pred,LRT_score,M-CAP_pred,M-CAP_rankscore,M-CAP_score,MetaSVM_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_score,PROVEAN_converted_rankscore,PROVEAN_pred,PROVEAN_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_score,REVEL_score,SIFT4G_converted_rankscore,SIFT4G_pred,SIFT4G_score,SIFT_converted_rankscore,SIFT_pred,SIFT_score,TSL,Uniprot_acc,Uniprot_entry,VEP_canonical,VEST4_rankscore,VEST4_score,aaalt,aapos,aaref,cds_strand,codonpos,genename,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,refcodon,rs_dbSNP,Interpro_domain,ExAC_AF,ESP6500_EA_AF,1000Gp3_AF,1000Gp3_AC,ExAC_AC" -db',
    #parametros = "annotate",
    parametros = 'dbnsfp -v -db /data/new_dbs/annot/dbNSFP',
    input_vcf = step6_Snpsift_GWASCat.salida_Snpsift,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step7_dbNSFP"
#
}

#####step8  annovar
call annovar {
input:
#vcf_in = step6_Snpsift_GWASCat.salida_Snpsift,
vcf_in = step7_dbNSFP.salida_Snpsift,

out_prefix = samplename1,
#File out_prefix
toolpath = toolpath

}


call intervar_postprocessing {
    input:
        vcfinput = annovar.multianno_vcf,
        multianno_txt = annovar.multianno_txt,
        toolpath = toolpath,
        samplename1 = samplename1

}

#Step 8: Annotate with VarType
#### no lleva database

call Snpsift_nodb as step8_VarType{
input:
    samplename1 = samplename1,
    parametros = "varType",
    input_vcf = intervar_postprocessing.salida_intervar,
    #input_vcf = step7_dbNSFP.salida_Snpsift, #####salida annovar
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step8_VarType"

}

#Step 9: Annotate with Exome Variant Server
#call Snpsift as step9_EVS{
#input:
#    samplename1 = samplename1,
#    parametros = "annotate -v -info MAF",
#    input_vcf = step8_VarType.salida_Snpsift,
#    toolpath = toolpath,
#    java_heap_memory_initial = java_heap_memory_initial,
#    nombre_step = "step9_EVS"

#}

#Step 10: Annotate with PhastCons
#http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phastCons100way/
#/data/

call Snpsift_nodb as step10_PhastCons{
input:
    samplename1 = samplename1, #### crearlo para la 38, bajando
    parametros = "phastCons -v /data/new_dbs/annot/phastCons/",
    input_vcf = step8_VarType.salida_Snpsift, # step9_EVS.salida_Snpsift,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "step10_PhastCons"

}

#Step 11: Annotate with CADD (Combined Annotation Dependent Depletion) - 1000 Genome variants

#call Snpsift as step11_CADD{
#input:
#    samplename1 = samplename1,
#    parametros = "annotate",
#    input_vcf = step10_PhastCons.salida_Snpsift,
#    toolpath = toolpath,
#    java_heap_memory_initial = java_heap_memory_initial,
#    nombre_step = "step11_CADD"
#
#}

##step 11 bis, new dbnsfp 
#call Snpsift as step11_dbNSFP{
#input:
#   samplename1 = samplename1,
#   parametros = "annotate -info", #-v -info RefSeq,Ensembl,RefSeq_region,RefSeq_gene,RefSeq_functional_consequence,RefSeq_id_c.change_p.change,Ensembl_region,Ensembl_gene,Ensembl_functional_consequence,Ensembl_id_c.change_p.change,ada_score,rf_score",
#   input_vcf = step10_PhastCons.salida_Snpsift,
#   toolpath = toolpath,
#   java_heap_memory_initial = java_heap_memory_initial,
#   nombre_step = "step11_dbNSFP"

#}


#Step 12: Annotate with ClinVar
#call Snpsift as step12_clinVar{
#input:
#    samplename1 = samplename1,
#    #parametros = "annotate -v -info CLNHGVS,CLNALLE,CLNSRC,CLNORIGIN,CLNSRCID,CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNACC",
#    parametros = "annotate",
#    input_vcf = step10_PhastCons.salida_Snpsift,
#    #input_vcf = step11_dbNSFP.salida_Snpsift,
#    toolpath = toolpath,
#    java_heap_memory_initial = java_heap_memory_initial,
#    nombre_step = "step12_clinVar"

#}

#Step 13: Annotate with PharmGKB
#call Snpsift as step13_pharmGKB{
#input:
#    samplename1 = samplename1,
#    parametros = "annotate -v -info PGKB_INDEX,PGKB_GENE,PGKB_DRUG,PGKB_TYPE,PGKB_EVIDENCE,PGKB_DISEASE,PGKB_RACE",
#    input_vcf = step12_clinVar.salida_Snpsift,
#    toolpath = toolpath,
#    java_heap_memory_initial = java_heap_memory_initial,
#    nombre_step = "step13_pharmGKB"




 #call hnrg_freq {

 #input:
 #   input_vcf = step12_clinVar.salida_Snpsift,   # step13_pharmGKB.salida_Snpsift, 
 #   samplename1 = samplename1,
 #   toolpath = toolpath,
 #   nombre_step = "step14_HNRG_FREQ"

    
 #}


#Step 14: Annotate with ExAC
#Step 12: Annotate with ClinVar
call Snpsift as step12_clinVar_final_annot{
input:
    samplename1 = samplename1,
    #parametros = "annotate -v -info CLNHGVS,CLNALLE,CLNSRC,CLNORIGIN,CLNSRCID,CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNACC",
    parametros = "annotate",
    input_vcf = step10_PhastCons.salida_Snpsift,
    #input_vcf = step11_dbNSFP.salida_Snpsift,
    toolpath = toolpath,
    java_heap_memory_initial = java_heap_memory_initial,
    nombre_step = "final_annot_"+pipeline_version

}


#call Snpsift as final_annot{
#input:

#    samplename1 = samplename1,
#    parametros = "annotate -v -info AN_Adj,AC_Adj,AC_Het,AC_Hom,AC_Hemi,POPMAX,VQSLOD,GQ_MEAN,GQ_STDDEV,HWP",
#    input_vcf = step12_clinVar.salida_Snpsift, # hnrg_freq.out_vcfanno,
#    toolpath = toolpath,
#    java_heap_memory_initial = java_heap_memory_initial,
#    nombre_step = "final_annot_"+pipeline_version

#}

call exon_distance {
    input:
    vcf_ok = step12_clinVar_final_annot.salida_Snpsift, #step12_clinVar.salida_Snpsift, #final_annot.salida_Snpsift,
    exon_coord = exon_coordinates,
    #exon_coordinates_to_lib = exon_coordinates_to_lib,
    sample_name = samplename1,
    path_save = path_save
}

 call symlink_important_files {
         input:
        output_to_save = step12_clinVar_final_annot.salida_Snpsift,#final_annot.salida_Snpsift,
        output_to_save2 = exon_distance.exon_dist,
        path_save = path_save
    }


####tsv file to merge excel file.
output {
File vcf_exon_distance = exon_distance.exon_dist 
}
}

