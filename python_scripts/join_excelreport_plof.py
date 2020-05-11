#!/usr/bin/python


import pandas as pd 
import numpy as np
import xlsxwriter
import os
import sys
import openpyxl 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p','--plof_file', help='input constrait file from gnomad')
parser.add_argument('-e','--excel_file', help='excel file with variants from pipelineHNRG')
parser.add_argument('-o','--output_file', help='an excel file merged, variants and plof constraint. must be .xlsx extension')
args = parser.parse_args()

plof = args.plof_file
excel = args.excel_file
out = args.output_file

plof_df = pd.read_csv(plof,header =0,sep ='\t')

excel_file_variants = pd.read_excel(excel, sheet_name= 'Variants')
excel_file_coverage = pd.read_excel(excel, sheet_name= 'ExonCoverage')


merged = pd.merge(left = excel_file_variants, right = plof_df ,left_on = 'Gene.knownGene', right_on = 'gene' , how = 'left' )

columnas_variants = len(excel_file_variants.columns)
columnas_merged = len(merged.columns)
#columnas_plof = len(plof_df.columns)

hnrg_soft = {'Plataforma':['Illumina Nextseq500'],'Referencia':['GRCh37.75'],'Software pipeline HNRG':['Gatk-package-4.0.8.1','fastp 0.20.0','Bwa 0.7.15-r1140','bedtools v2.28.0','SAMTOOLS_VERSION 1.9','SnpEff 4.3t','VCFtools (0.1.17)','cromwell-37'],'Anotacion':['dbSNP:All_20180423 v151','indels_site: Mills_and_1000G_gold_standard.indels.b37','hapmap_3.3','1000G_omni2.5.b37','Gwas_catalog_v1.0.2-associations_e96_r2019-04-06','ESP6500SI-V2-SSA137.snps_indels','Clinvar_20190219','PharmGKBvcf_2016','ExAC.r1.sites.vep'],'Annovar':['refGene','avsnp150','esp6500siv2_all','1000g2015aug_all','exac03','gnomad_exome','gnomad_genome','clinvar_20180603','dbscsnv11','dbnsfp35a,rmsk','tfbsConsSites','cytoBand','wgRna','targetScanS','genomicSuperDups','dgvMerged','gwasCatalog','ensGene','knownGene','intervar_20180118']}
hnrg = pd.DataFrame(list(hnrg_soft.values()), index=hnrg_soft.keys()).replace(np.nan, '', regex=True)

merged_ok = merged.iloc[:,np.r_[0:columnas_variants-19,columnas_variants:columnas_merged,132:columnas_variants]]


with pd.ExcelWriter(out,engine = 'xlsxwriter') as writer:
    merged_ok.to_excel(writer, sheet_name= 'Variants',index = False)
    excel_file_coverage.to_excel(writer, sheet_name= 'ExonCoverage', index=False)
    hnrg.to_excel(writer, sheet_name= 'Software', index = True)