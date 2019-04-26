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
    dp = multianno[sample].str.split(':',expand=True)[2].apply(pd.to_numeric)
    multianno['DP'] = dp
    max_alt = multianno[sample].str.split(':',expand=True)[1].str.split(',',expand = True).applymap(pd.to_numeric).max(1)
    multianno['fqmax_alt'] = np.round(max_alt/dp,2)
    return multianno



def process_InterVar(multianno):
    multianno['InterVarEvidence'] = multianno.apply(lambda x: list(x.keys()[x==1]) if any(x==1) else np.nan,axis =1)
    multianno.drop(['PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2',
           'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1',
           'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6',
           'BP7'],axis =1,inplace = True)

    multianno.rename(columns={'InterVar_automated':'InterVarVeredict'},inplace = True)
    return(multianno)



#This assume that your vcf file contain only one sample. 
multianno=pd.read_table(multianno_tsv)
sample = multianno.columns[-1]

# this vcf is multisample
vcf=read_vcf(vcf_file)
vcf = vcf.iloc[:,~vcf.columns.isin(['QUAL','FILTER','INFO','FORMAT',sample])]
vcf.rename(columns={'ALT':'ALTERNATIVES'},inplace=True)

multianno = process_genotipo(multianno) 
multianno = process_InterVar(multianno)
df = pd.merge(multianno,vcf,how='left',on=['#CHROM','POS','ID','REF'])

## reorganize columns
cols = df.columns
firstCols =  ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GENOTIPO','DP','fqmax_alt',sample]
cols =  [c for c in cols if c not in firstCols]
cols = firstCols + cols
df = df.reindex(columns=cols)
df.to_csv(output,sep='\t',index = False)