
import argparse
import pandas as pd
from xlrd import open_workbook
import csv
import os

parser = argparse.ArgumentParser(prog='create_InterVarDB.py',description='Create vcf database for intervar annotation')

parser.add_argument('-i','--multianno_txt', help='raw annovar output in tabulated format')
parser.add_argument('-o','--output_vcf', default='./intervar.vcf')
args =  parser.parse_args()
inputfile = args.multianno_txt
outfile = args.output_vcf



def load_annovar_txt_for_vcf_body(multianno_txt):
    cols1 = pd.read_csv(multianno_txt,sep='\t',header=None,nrows=1).iloc[0,:][:-1]
    cols2 = ['#CHROM','POS','ID','REF','ALT']
    body = pd.read_csv(multianno_txt,sep='\t',header=None,skiprows=1)

    body1 = body.iloc[:,0:len(cols1)]
    body2 = body.iloc[:,(len(cols1)+3):(len(cols1)+8)]
    Body = pd.concat([body1,body2],axis = 1)
    Body.columns = list(cols1) + cols2
    return Body



def main():
    header = ['##fileformat=VCFv4.2','##fileDate=2019-01-23',\
            '##source=HNRG','##reference=GRCh37',\
            '##INFO=<ID=InterVarVeredict,Number=1,Type=String,Description="InterVar_auto">',\
            '##INFO=<ID=InterVarEvidence,Number=1,Type=String,Description="Joined evidence codes ACMG categories">'
            ]

    multianno = load_annovar_txt_for_vcf_body(inputfile)
    multianno_hnrg = process_InterVar(multianno)
    multianno_hnrg = multianno_hnrg[[u'#CHROM',u'POS',u'ID',u'REF',u'ALT',u'InterVarVeredict',u'InterVarEvidence']]
  
    body = multianno_hnrg.sort_values(by=[u'#CHROM','POS'],ascending = [True,True])
    info = 'InterVarVeredict='+multianno_hnrg.InterVarVeredict+';'+'InterVarEvidence='+multianno_hnrg.InterVarEvidence

    body['QUAL']='.'
    body['FILTER']='PASS'
    body['INFO']=info
    body.drop(['InterVarVeredict','InterVarEvidence'],axis = 1,inplace = True)

with open(outfile,'wb') as f:
    for x in header:
        f.write('%s\n'%x)
    f.close()
body.to_csv(outfile,sep ='\t',mode='a',index = False)    
