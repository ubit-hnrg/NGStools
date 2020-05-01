#!/usr/bin/python


import pandas as pd 
import numpy as np
import xlsxwriter
import os
import sys
import openpyxl 
import argparse

parser = argparse.ArgumentParser()#prog='bam_sex.py',description='get sample sexuality from bam file.', usage='%(prog)s  --bam --output_file')
parser.add_argument('-d','--distance_file', help='input distance file from bedtools closest')
parser.add_argument('-e','--excel_file', help='excel file with variants from pipelineHNRG')
parser.add_argument('-o','--output_file', help='an excel file merged, variants and exon distance. must be .xlsx extension')
args = parser.parse_args()

dist = args.distance_file
excel = args.excel_file
out = args.output_file


dist = pd.read_csv(dist,header =None,sep ='\t')#,usecols = [0,1,2,3,4,5,9])
dist.rename(columns={0: 'CHROM',1:'pos',2:'e_start',3:'e_end',4:'gen_symbol',5:'exon_number',6:'strand',7:'dist'}, inplace=True)
dist['exon_id']=dist[['gen_symbol','strand','exon_number','e_start','e_end']].apply(lambda x: '_'.join(x.values.astype(str)), axis=1)

excel_file_variants = pd.read_excel(excel, sheet_name= 'Variants')
excel_file_coverage = pd.read_excel(excel, sheet_name= 'ExonCoverage')

hnrg_soft = {'Plataforma':['Illumina Nextseq500'],'Referencia':['GRCh37.75'],'Software pipeline HNRG':['Gatk-package-4.0.8.1','fastp 0.20.0','Bwa 0.7.15-r1140','bedtools v2.28.0','SAMTOOLS_VERSION 1.9','SnpEff 4.3t','VCFtools (0.1.17)','cromwell-37'],'Anotacion':['dbSNP:All_20180423 v151','indels_site: Mills_and_1000G_gold_standard.indels.b37','hapmap_3.3','1000G_omni2.5.b37','Gwas_catalog_v1.0.2-associations_e96_r2019-04-06','ESP6500SI-V2-SSA137.snps_indels','Clinvar_20190219','PharmGKBvcf_2016','ExAC.r1.sites.vep'],'Annovar':['refGene','avsnp150','esp6500siv2_all','1000g2015aug_all','exac03','gnomad_exome','gnomad_genome','clinvar_20180603','dbscsnv11','dbnsfp35a,rmsk','tfbsConsSites','cytoBand','wgRna','targetScanS','genomicSuperDups','dgvMerged','gwasCatalog','ensGene','knownGene','intervar_20180118']}
hnrg = pd.DataFrame(list(hnrg_soft.values()), index=hnrg_soft.keys()).replace(np.nan, '', regex=True)

merged = pd.merge(left = excel_file_variants, right = dist, left_on = 'POS', right_on = 'pos' )

columnas_variants = len(excel_file_variants.columns)
columnas_merged = len(merged.columns)
merged_ok =merged.iloc[:,np.r_[0:19,141,143,columnas_merged-2,columnas_merged-1,20:140,142,144:columnas_variants]]

with pd.ExcelWriter(out,engine = 'xlsxwriter') as writer:
    merged_ok.to_excel(writer, sheet_name= 'Variants',index = False)
    excel_file_coverage.to_excel(writer, sheet_name= 'ExonCoverage', index=False)
    hnrg.to_excel(writer, sheet_name= 'Software', index = True)