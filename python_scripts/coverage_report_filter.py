#!/usr/bin/python

import pandas as pd
import os
import argparse

###############################################################################
#del archivo de cobertura del bam sin restricciones, se eliminan las lecturas con 0 cobertura. 
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-r','--cov_report', help=' "reporte de cobertura sin restricciones intersectbed with TSO bed')
parser.add_argument('-l','--library', help=' "library bed: Ej: TSO')

parser.add_argument('-o1','--output_not_in_tso')
parser.add_argument('-o2','--output_offtarget')

#parser.add_argument('-s','--samplename')

args =  parser.parse_args()

#ffile = sys.argv[1]
#ffile = args.ensembl_coverage_histogram  # input file
input_file = args.cov_report
lib = args.library
#sample_name = args.samplename
exones_tso_no_cubiertos = args.output_not_in_tso
exones_offtarget = args.output_offtarget

    
 
    
def main():
    not_in_tso = pd.read_csv(input_file, header = None, sep = '\t')
    not_in_tso.rename(columns={0:'chr', 1:'start', 2:'end',3:'gene', 4:'transcriptID', 5:'exonNumber', 6:'strand', 7:'IntervalLength', 8:'dp>=1', 9:'dp>=10', 10:'dp>=20', 11:'dp>=30', 12:'dp>=50',13:'dp>=100'}, inplace = True)

    tso_lib = pd.read_csv(lib, header = None, sep = '\t')
    tso_lib.rename(columns={0:'CHROM',1:'AMP_START',2:'AMP_END',3:'Name'},inplace = True)
    tso_lib[['GENE','chr','start','end']] = tso_lib['Name'].str.split('.',expand= True)
    tso_lib.chr.fillna(tso_lib.CHROM, inplace=True)
    tso_lib.start.fillna(tso_lib.AMP_START, inplace=True)
    tso_lib.end.fillna(tso_lib.AMP_END, inplace=True)
    tso_lib.start = tso_lib.start.astype(int)

    offtarget = not_in_tso[(not_in_tso['dp>=10'] != 0.0) | (not_in_tso['dp>=20'] != 0.0) | (not_in_tso['dp>=30'] != 0.0) | (not_in_tso['dp>=50'] != 0.0) | (not_in_tso['dp>=100']!= 0.0)]
    
    ###lecturas que tienen genes en el TSO pero se cubren los exones cob >0

    offtarget[offtarget.gene.isin(tso_lib.GENE)].to_csv(exones_tso_no_cubiertos,sep = '\t',index = False)
    ### lecturas que no estan en el tso con cob >0
    offtarget[~offtarget.gene.isin(tso_lib.GENE)].to_csv(exones_offtarget,sep = '\t',index = False)

   

if __name__ == "__main__":
    main()
