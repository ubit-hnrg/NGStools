#!/usr/bin/python


import pandas as pd 
import numpy as np
import xlsxwriter
import os
import sys
import openpyxl 
import argparse

parser = argparse.ArgumentParser()#prog='join_annovar_exon_dist.py',description='creates excel file joining distance exon file with annovar variants', usage='%(prog)s  -d -a -o')
parser.add_argument('-d','--distance_file', help='input exon distance file from functional_annot')
parser.add_argument('-a','--annovar_resctricted', help='annovar variants from Annovar')
parser.add_argument('-o','--output_file', help='an excel file merged, variants and exon distance. must be .xlsx extension')
args = parser.parse_args()

path_dist = args.distance_file
path_variants = args.annovar_resctricted
out = args.output_file


###distancia from vcf annotado.
dist = pd.read_csv(path_dist,header =None,sep ='\t')#,usecols = [0,1,2,3,4,5,9])
dist.rename(columns={0: 'CHROM',1:'pos',2:'e_start',3:'e_end',4:'ENSEMBL_ID',5:'gen_symbol',6:'exon_number',7:'strand',8:'dist'}, inplace=True)
dist['exon_id']=dist[['gen_symbol','strand','exon_number','e_start','e_end']].apply(lambda x: '_'.join(x.values.astype(str)), axis=1)

###variants from annovar
variants = pd.read_csv(path_variants, header= 0, sep ='\t')
##merge
merged = pd.merge(left = variants, right = dist, left_on = 'POS', right_on = 'pos', how= 'left' )
columnas_variants = len(variants.columns)
columnas_merged = len(merged.columns)

##rearrange 
#merged_ok =merged.iloc[:,np.r_[0:19,141,143,columnas_merged-2,columnas_merged-1,20:140,142,144:columnas_variants]]
merged_ok =merged.iloc[:,np.r_[0:19,89,93,94,20:88,90:92 ]] #columnas_merged-2,columnas_merged-1,20:140,142,144:columnas_variants]]

merged_ok.to_csv(out, sep = '\t', index = False)
