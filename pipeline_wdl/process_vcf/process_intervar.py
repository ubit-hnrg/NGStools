import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(prog='left_join_multianno_and_multisampleVCF', usage='%(prog)s [options] > outputfile.tsv')

    parser.add_argument('--multianno_tsv', help='modified annovar output (with removed unammed fields)')
    parser.add_argument('--output', help='output file')

    args = parser.parse_args()
    return args

def process_InterVar(multianno):
    evidence_cols = ['PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2',
           'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1',
           'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6',
           'BP7']
    multianno['InterVarEvidence'] =multianno[evidence_cols].apply(lambda x: ','.join(list(x.keys()[x=='1'])) if any(x=='1') else np.nan,axis =1)
    multianno.drop(evidence_cols,axis =1,inplace = True)

    veredict = multianno['InterVar_automated']
    multianno.drop(['InterVar_automated'],inplace = True,axis=1)
    multianno['InterVarVeredict'] = veredict
    return(multianno)

def main():
    args = parse_args()
    multianno_tsv=args.multianno_tsv
    output=args.output
    multianno=pd.read_table(multianno_tsv)
    
    multianno = process_InterVar(multianno)
    db = multianno[['#CHROM','POS','ID','REF','ALT','InterVarVeredict','InterVarEvidence']]
    db.to_csv(output,sep='\t',index=False)

if __name__ == "__main__":
    main()