import pandas as pd
import argparse
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(prog='left_join_multianno_and_multisampleVCF', usage='%(prog)s [options] > outputfile.tsv')

    parser.add_argument('--multianno_txt', help='modified annovar output (with removed unammed fields)')
    parser.add_argument('--output', help='output file')

    args = parser.parse_args()
    return args

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

def load_annovar_txt(multianno_txt):
    cols = pd.read_csv(multianno_txt,sep='\t',header=None,nrows=1).iloc[0,:][:-1]
    body = pd.read_csv(multianno_txt,sep='\t',header=None,skiprows=1).iloc[:,0:len(cols)]
    body.columns = cols
    return body

def main():
    args = parse_args()
    multianno_txt=args.multianno_txt
    output=args.output
    multianno=load_annovar_txt(multianno_txt)
    multianno = process_InterVar(multianno)
    #print multianno.columns
    cols = ['Chr','Start','End','Ref','Alt','InterVarVeredict','InterVarEvidence']
    multianno_db = multianno[cols]
    multianno_db.to_csv(output,index=False)

if __name__ == "__main__":
    main()