
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(prog='create_InterVarDB.py',description='Create vcf database for intervar annotation')
parser.add_argument('-i','--multianno_txt', help='raw annovar output in tabulated format')
parser.add_argument('-o','--output_vcf', default='./intervar.vcf')
args =  parser.parse_args()
inputfile = args.multianno_txt
outfile = args.output_vcf



def load_annovar_txt_for_vcf_body(multianno_txt):
    cols1 = pd.read_csv(multianno_txt,sep='\t',header=None,nrows=1).iloc[0,:][:-13] 
    cols2 = ['#CHROM','POS','ID','REF','ALT']
    body = pd.read_csv(multianno_txt,sep='\t',header=None,skiprows=1)
    body1 = body.iloc[:,0:len(cols1)]
    body2 = body.iloc[:,(len(cols1)+3):(len(cols1)+8)]
    Body = pd.concat([body1,body2],axis = 1)
    Body.columns = list(cols1) + cols2
    return Body

def process_InterVar(multianno):
    evidence_cols = ['PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2',
           'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1',
           'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6',
           'BP7']
    multianno['InterVarEvidence'] =multianno[evidence_cols].apply(lambda x: ','.join(list(x.keys()[x=='1'])) if any(x=='1') else np.nan,axis =1)
    #multianno.drop(evidence_cols,axis =1,inplace = True)

    veredict = multianno['InterVar_automated']
    multianno.drop(['InterVar_automated'],inplace = True,axis=1)
    multianno['InterVarVeredict'] = veredict
    return(multianno)


def main():
    header = ['##fileformat=VCFv4.2','##fileDate=2024-01-23',\
            '##source=HNRG','##reference=GRCh38.p14',\
            '##INFO=<ID=InterVarVeredict,Number=1,Type=String,Description="InterVar_auto">',\
            '##INFO=<ID=InterVarEvidence,Number=1,Type=String,Description="Joined evidence codes ACMG categories">'
            ]

    multianno = load_annovar_txt_for_vcf_body(inputfile)
    multianno_hnrg = process_InterVar(multianno)
    multianno_hnrg = multianno_hnrg[[u'#CHROM',u'POS',u'ID',u'REF',u'ALT',u'InterVarVeredict',u'InterVarEvidence']]
  
    body = multianno_hnrg.sort_values(by=[u'#CHROM','POS'],ascending = [True,True])
    info = 'InterVarVeredict='+multianno_hnrg.InterVarVeredict+';'+'InterVarEvidence='+multianno_hnrg.InterVarEvidence

    body['QUAL'] = '.'
    body['FILTER'] = 'PASS'
    body['INFO'] = info
    body['INFO'] = body['INFO'].str.replace('.','-')
    body.drop(['InterVarVeredict','InterVarEvidence'],axis = 1,inplace = True)
    
   # with open(outfile,'wb') as f:
   #     for x in header:
   #         f.write('%s\n'%x)
   #     f.close()
   # body.to_csv(outfile,sep ='\t',mode='a',index = False)    
   
    with open(outfile, 'w') as f:  
        for x in header:
            f.write('%s\n' % x)
    body.to_csv(outfile, sep='\t', mode='a', index=False)

if __name__ == "__main__":
    main()
