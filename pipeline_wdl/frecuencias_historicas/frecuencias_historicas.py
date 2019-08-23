#!/usr/bin/python
import pandas as pd 
import numpy as np
import csv
import argparse

parser = argparse.ArgumentParser(prog='frecuencias_historicas.py',description='Extract historic Allele frequency', usage='%(prog)s  --freq_report --output_file')
parser.add_argument('-i','--freq_report', help='input file from VCFTOOLS --freq')
parser.add_argument('-o','--output_file', help='a tsv file with allele frequency sorted by #CHR&POS')
args = parser.parse_args()

path = args.freq_report
out = args.output_file

allele_freqs = pd.read_csv(path)
allele_freqs.columns = ['freqs']
#nuevo dataset fqr, con split se separan los campos para tener nuevas columnas y expand muestra todas las columnas.
fqr = allele_freqs.freqs.str.split('\t',expand = True)
fqr.rename(columns = {0:'CHROM',1:'POS',2:'N_ALLELS',3:'N_CHR'},inplace =True)
allele_cols = fqr.columns.values[4:]


melted = pd.melt(fqr, id_vars=['CHROM','POS'],value_vars=allele_cols)
melted = melted[~melted.value.isnull()]
melted[['alt','freq']] = melted.value.str.split(':',expand =True)

melted.POS = melted.POS.map(int)
melted.sort_values(by=['CHROM','POS'],ascending = [True,True],inplace =True)
melted.drop(['variable','value'],axis=1,inplace =True)

melted.to_csv(out,index=False,sep ='\t')