#!/usr/bin/python

import pandas as pd
import os, sys, argparse
from pandas import ExcelWriter

  


def leer_excel(path_excel):
    '''
    lee el excel y guarda las hojas de variantes y cobertura por exon
    '''
    variantes = pd.read_excel(path_excel,sheet_name='Variants',index_col='Gene.knownGene')
    exon_cov = pd.read_excel(path_excel,sheet_name='ExonCoverage')
    samplename = os.path.basename(path_excel).replace('_variants_V2.xlsx','')

    return variantes, exon_cov, samplename

def busqueda_genes_acmg(path_acmg,variantes,exon_cov,name):
    genes_acmg = pd.read_csv(path_acmg, header = 0, sep = '\t')
    genes = genes_acmg[['Gene','Inheritance ','SF List Version']]
    genes_only = genes['Gene']
    
    ###orden de columnas
    nuevas_col=['Func.knownGene','ExonicFunc.refGene','AAChange.refGene','QUAL','DP','dist','exon_id','Func.refGene','Func.ensGene','cytoBand','ExonicFunc.ensGene','ExonicFunc.knownGene','AAChange.knownGene','Gene.ensGene','InterVarEvidence','InterVarVeredict','CLNSIG','CLNREVSTAT','CLNDN','CLNDISDB','#CHROM','REF', 'ALT','POS','ID','FILTER','INFO','FORMAT','GENOTIPO',name,'fqmax_alt','CLNALLELEID','esp6500siv2_all','1000g2015aug_all','ExAC_ALL','ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN', 'ExAC_NFE','ExAC_OTH','ExAC_SAS','gnomAD_exome_ALL','gnomAD_exome_AFR','gnomAD_exome_AMR','gnomAD_exome_ASJ','gnomAD_exome_EAS','gnomAD_exome_FIN','gnomAD_exome_NFE','gnomAD_exome_OTH','gnomAD_exome_SAS','gnomAD_genome_ALL','gnomAD_genome_AFR','gnomAD_genome_AMR','gnomAD_genome_ASJ','gnomAD_genome_EAS','gnomAD_genome_FIN','gnomAD_genome_NFE','gnomAD_genome_OTH', 'Chr','Start','End','Ref','Alt', 'GeneDetail.knownGene', 'GeneDetail.refGene','avsnp150','dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','SIFT_score','SIFT_converted_rankscore','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_rankscore','Polyphen2_HDIV_pred','Polyphen2_HVAR_score','Polyphen2_HVAR_rankscore','Polyphen2_HVAR_pred','LRT_score','LRT_converted_rankscore','LRT_pred','MutationTaster_score','MutationTaster_converted_rankscore','MutationTaster_pred','MutationAssessor_score','MutationAssessor_score_rankscore','MutationAssessor_pred','FATHMM_score','FATHMM_converted_rankscore','FATHMM_pred','PROVEAN_score','PROVEAN_converted_rankscore','PROVEAN_pred','VEST3_score','VEST3_rankscore','MetaSVM_score','MetaSVM_rankscore','MetaSVM_pred','MetaLR_score', 'MetaLR_rankscore', 'MetaLR_pred','M-CAP_score','M-CAP_rankscore','M-CAP_pred','REVEL_score','REVEL_rankscore','MutPred_score','MutPred_rankscore','CADD_raw','CADD_raw_rankscore','CADD_phred','DANN_score','DANN_rankscore','fathmm-MKL_coding_score','fathmm-MKL_coding_rankscore','fathmm-MKL_coding_pred','Eigen_coding_or_noncoding','Eigen-raw','Eigen-PC-raw','GenoCanyon_score','GenoCanyon_score_rankscore','integrated_fitCons_score','integrated_fitCons_score_rankscore','integrated_confidence_value','GERP++_RS','GERP++_RS_rankscore','phyloP100way_vertebrate','phyloP100way_vertebrate_rankscore','phyloP20way_mammalian','phyloP20way_mammalian_rankscore','phastCons100way_vertebrate','phastCons100way_vertebrate_rankscore','phastCons20way_mammalian', 'phastCons20way_mammalian_rankscore','SiPhy_29way_logOdds','SiPhy_29way_logOdds_rankscore','Interpro_domain', 'GTEx_V6p_gene', 'GTEx_V6p_tissue','rmsk', 'tfbsConsSites', 'wgRna', 'targetScanS', 'genomicSuperDups','dgvMerged', 'gwasCatalog', 'GeneDetail.ensGene','ALTERNATIVES']

    variantes_2do = variantes[nuevas_col][variantes.index.isin(genes_only)].sort_index(ascending=True)
    #exones_2do = exon_cov[exon_cov.index.isin(genes_only)]

    merged = pd.merge(exon_cov,genes, how= "inner", left_on="gene", right_on="Gene", right_index = False)[['gene','transcriptID','exonNumber','start','end','strand','IntervalLength','dp>=1','dp>=10', 'dp>=20','dp>=30','dp>=50','dp>=100','Inheritance ','SF List Version']]
    return genes_acmg, variantes_2do, merged

def crear_excel(df1,df2,df3, path_out):

    '''
    guardamos a excel
    otra forma:
    with pd.ExcelWriter('nombreEXCEL.xlsx',engine = 'xlsxwriter') as writer:
    df1.to_excel(writer, sheet_name= 'hoja1',index = False)
    df2.to_excel(writer, sheet_name= 'hoja2',index = False)
    '''
    
    writer = ExcelWriter(path_out)
    
    df1.to_excel(writer,'Variantes_2dos')
    df2.to_excel(writer,'Cobertura_2dos',index=False)
    df3.to_excel(writer,'ACMG_V3',index=False)
    
    writer.save()


def main(argv):

    parser = argparse.ArgumentParser(description= "A partir del excel de variantes se filtra por los genes de hallazgos secundarios ACMG")
    parser.add_argument('-v','--variantes', help='Archivo de variantes from ANNOVAR *_variants*.xlsx')
    parser.add_argument('-a','--acmg', help = 'archivo genes hallazgos secundarios segun ACMG')
    parser.add_argument('-o','--output', help = 'output path')
    args = parser.parse_args()

    if(len(sys.argv) > 1):
        if(not os.path.isfile(args.variantes)):
            print(sys.argv[1],"no se reconoce el primer archivo")
            sys.exit(1)
        elif(not os.path.isfile(args.acmg)):
            print(sys.argv[2],"no se reconoce el 2do archivo")
            sys.exit(1)
    
    excel = args.variantes#sys.argv[1]
    acmg = args.acmg#sys.argv[2]
    out_dir = args.output#sys.argv[3]

    var, exones, nombre_muestra = leer_excel(excel)
    genes_acmg, variantes_2do, exones_2do = busqueda_genes_acmg(acmg,var,exones, nombre_muestra)

    #guardamos
    crear_excel(variantes_2do,exones_2do,genes_acmg,out_dir)



if __name__ == '__main__':
    main(sys.argv[1:])


