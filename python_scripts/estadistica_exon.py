#!/usr/bin/python

import pandas as pd 
import numpy as np
import glob
import argparse
import csv


parser = argparse.ArgumentParser(prog='estadistica_exon.py',description='get statistics from *_variants.xlsx exon files.', usage='%(prog)s  --freq_report --output_file')
parser.add_argument('-d','--directory', help='path of the *_variants.xlsx files')
parser.add_argument('-o','--output_file', help='a tsv file with exon statistics ')
args = parser.parse_args()

directory = args.directory
out = args.output_file



def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_


#directory = r'/home/usuario/Escritorio/reporte_exones/'
files = glob.glob(directory + "*_variants.xlsx")
#fields = ["geneSymbol","dp1","dp10","dp20","dp30"]
tablas = []
#df = pd.DataFrame()

for f in files:
    # Every file is put in a temp dataframe, and operations are performed
    temp_df = pd.read_excel(f, sheet_name= "ExonCoverage")
    tablas.append(temp_df)
tablas = pd.concat(tablas)
genes = tablas.groupby(['geneSymbol','exon_number'])
promedios = genes['dp1','dp10','dp20','dp30'].agg([np.min,np.mean,np.std,percentile(25),percentile(75)])

resultado = pd.merge(tablas[['geneSymbol','exon_number','chr','exon_start', 'exon_end']],promedios, on=('geneSymbol','exon_number'), how='left')

resultado.columns=['geneSymbol','exon_number','chr','exon_start','exon_end','dp1_MIN','dp1_MEAN','dp1_STD','dp1_q25','dp1_q75','dp10_MIN','dp10_MEAN','dp10_STD', 'dp10_q25','dp10_q75','dp20_MIN','dp20_MEAN','dp20_STD','dp20_q25', 'dp20_q75','dp30_MIN','dp30_MEAN','dp30_STD','dp30_q25','dp30_q75']
resultado=resultado[['geneSymbol','exon_number','chr','exon_start','exon_end','dp1_MEAN','dp1_STD','dp1_MIN','dp1_q25','dp1_q75','dp10_MEAN','dp10_STD','dp10_MIN','dp10_q25','dp10_q75','dp20_MEAN','dp20_STD','dp20_MIN','dp20_q25', 'dp20_q75','dp30_MEAN','dp30_STD','dp30_MIN','dp30_q25','dp30_q75']]

resultado.to_csv(out,sep = '\t', index=False)

