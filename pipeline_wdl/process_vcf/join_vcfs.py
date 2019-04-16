import pandas as pd
import argparse 
from StringIO import StringIO
import sys


parser = argparse.ArgumentParser(prog='left_join_multianno_and_multisampleVCF', usage='%(prog)s [options] > outputfile.tsv')

parser.add_argument('--multianno_tsv', help='modified annovar output (with removed unammed fields)')
parser.add_argument('--vcf_multisample', help='multisample vcf file')

args = parser.parse_args()


multianno_tsv=args.multianno_tsv
vcf_file=args.vcf_multisample

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df= pd.read_csv(StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},sep='\t')
    return(df)


multianno=pd.read_table(multianno_tsv)
sample = multianno.columns[-1]
vcf=read_vcf(vcf_file)
vcf = vcf.iloc[:,~vcf.columns.isin(['QUAL','FILTER','INFO','FORMAT',sample])]

df = pd.merge(multianno,vcf,how='left',on=['#CHROM','POS','ID','REF','ALT'])
sys.stdout.write(df)
