import pandas as pd
import argparse 
from StringIO import StringIO
import sys
import numpy as np


parser = argparse.ArgumentParser(prog='left_join_multianno_and_multisampleVCF', usage='%(prog)s [options] > outputfile.tsv')

parser.add_argument('--multianno_tsv', help='modified annovar output (with removed unammed fields)')
parser.add_argument('--vcf_multisample', help='multisample vcf file')
parser.add_argument('--output', help='output file')

args = parser.parse_args()


multianno_tsv=args.multianno_tsv
vcf_file=args.vcf_multisample
output=args.output

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df= pd.read_csv(StringIO(''.join(lines))
        #dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               #'QUAL': str, 'FILTER': str, 'INFO': str}
               ,sep='\t')
    return(df)

def process_genotipo(multianno):
    multianno['GENOTIPO']=multianno[sample].str.split(':',expand=True)[0]
    dp = multianno[sample].str.split(':',expand=True)[2]
    multianno['DP'] = dp
    max_alt = multianno[sample].str.split(':',expand=True)[1].str.split(',').max()
    multianno['max_alt'] = max_alt
    multianno['fqmax_alt'] = multianno['max_alt']/multianno['DP']
    multianno['fqmax_alt'] = multianno['fqmax_alt'].round(2)
    return multianno


#This assume that your vcf file contain only one sample. 
multianno=pd.read_table(multianno_tsv)
sample = multianno.columns[-1]

# this vcf is multisample
vcf=read_vcf(vcf_file)
vcf = vcf.iloc[:,~vcf.columns.isin(['QUAL','FILTER','INFO','FORMAT',sample])]
vcf.rename(columns={'ALT':'ALTERNATIVES'},inplace=True)

multianno = process_genotipo(multianno) 
df = pd.merge(multianno,vcf,how='left',on=['#CHROM','POS','ID','REF'])

## reorganize columns
cols = df.columns
firstCols =  ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GENOTIPO','DP','fqmax_alt',sample]
cols =  [c for c in cols if c not in firstCols]
cols = firstCols + cols
df = df.reindex(columns=cols)
df.to_csv(output,sep='\t',index = False)