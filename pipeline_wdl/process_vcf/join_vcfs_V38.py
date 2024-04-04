import pandas as pd
import argparse
from io import StringIO  # Actualizado para Python 3
import sys
import numpy as np

parser = argparse.ArgumentParser(prog='left_join_multianno_and_multisampleVCF', usage='%(prog)s [options] > outputfile.tsv')

parser.add_argument('--multianno_tsv', help='modified annovar output (with removed unnamed fields)')
parser.add_argument('--vcf_multisample', help='multisample vcf file')
parser.add_argument('--output', help='output file')

args = parser.parse_args()

multianno_tsv = args.multianno_tsv
vcf_file = args.vcf_multisample
output = args.output

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')  # Actualizado para Python 3
    return df

def process_genotipo(multianno):
    sample = multianno.columns[-1]  # Mover la definición de sample aquí si es específica de esta función
    multianno['GENOTIPO'] = multianno[sample].str.split(':', expand=True)[0]
    dp = multianno[sample].str.split(':', expand=True)[2].apply(pd.to_numeric)
    multianno['DP'] = dp
    max_alt = multianno[sample].str.split(':', expand=True)[1].str.split(',', expand=True).applymap(pd.to_numeric).max(1)
    multianno['fqmax_alt'] = np.round(max_alt / dp, 2)
    return multianno

def process_InterVar(multianno):
    evidence_cols = ['PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2',
                     'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1',
                     'BS1', 'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6',
                     'BP7']
    multianno['InterVarEvidence'] = multianno[evidence_cols].apply(lambda x: ','.join(list(x.keys()[x == '1'])) if any(x == '1') else np.nan, axis=1)
    multianno.drop(evidence_cols, axis=1, inplace=True)

    veredict = multianno['InterVar_automated']
    multianno.drop(['InterVar_automated'], inplace=True, axis=1)
    multianno['InterVarVeredict'] = veredict
    return multianno

# Asumiendo que tu archivo VCF contiene solo una muestra.
multianno = pd.read_csv(multianno_tsv, sep='\t')  # Actualizado para usar read_csv con sep='\t'
sample = multianno.columns[-1]

# Este VCF es multisample
vcf = read_vcf(vcf_file)
vcf = vcf.iloc[:, ~vcf.columns.isin(['QUAL', 'FILTER', 'INFO', 'FORMAT', sample])]
vcf.rename(columns={'ALT': 'ALTERNATIVES'}, inplace=True)
vcf['#CHROM'] = vcf['#CHROM'].astype(str).str.strip()  # fuerzo a string la columna de cromosoma para evitar mistmaching entre enteros y chr de un mismo cromosoma

multianno = process_genotipo(multianno)
multianno = process_InterVar(multianno)
multianno['#CHROM'] = multianno['#CHROM'].astype(str).str.strip()  # igual que arriba

df = pd.merge(multianno, vcf, how='left', on=['#CHROM', 'POS', 'ID', 'REF'])

# Reorganizar columnas
cols = df.columns
firstCols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'GENOTIPO', 'DP', 'fqmax_alt', sample]
cols = [c for c in cols if c not in firstCols]
cols = firstCols + cols
df = df.reindex(columns=cols)
df.to_csv(output, sep='\t', index=False)
