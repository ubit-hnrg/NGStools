#!/usr/bin/python

import pandas as pd
import os, sys, argparse
from pandas import ExcelWriter

  


def leer_excel(path_excel):
    '''
    lee el excel y guarda las hojas de variantes y cobertura por exon
    '''
    variantes = pd.read_excel(path_excel,sheet_name='Variants',index_col='Gene.knownGene')
    exon_cov = pd.read_excel(path_excel,sheet_name='ExonCoverage')

    return variantes, exon_cov

def busqueda_genes_acmg(path_acmg,variantes,exon_cov):
    genes_acmg = pd.read_csv(path_acmg, header = 0, sep = '\t')
    genes = genes_acmg[['Gene','Inheritance ','SF List Version']]
    genes_only = genes['Gene']

    variantes_2do = variantes[variantes.index.isin(genes_only)]
    #exones_2do = exon_cov[exon_cov.index.isin(genes_only)]

    merged = pd.merge(exon_cov,genes, how= "inner", left_on="gene", right_on="Gene", right_index = False)[['gene','transcriptID','exonNumber','start','end','strand','IntervalLength','dp>=1','dp>=10', 'dp>=20','dp>=30','dp>=50','dp>=100','Inheritance ','SF List Version']]
    return genes_acmg, variantes_2do, merged

def crear_excel(df1,df2,df3, path_out):

    '''
    guardamos a excel
    '''
    
    writer = ExcelWriter(path_out)
    
    df1.to_excel(writer,'Variantes_2dos')
    df2.to_excel(writer,'Cobertura_2dos')
    df3.to_excel(writer,'ACMG_V3')
    
    writer.save()


def main(argv):

    if(len(sys.argv) > 1):
        if(not os.path.isfile(sys.argv[1])):
            print(sys.argv[1],"no se reconoce el primer archivo")
            sys.exit(1)
        elif(not os.path.isfile(sys.argv[2])):
            print(sys.argv[2],"no se reconoce el 2do archivo")
            sys.exit(1)
    
    excel = sys.argv[1]
    acmg = sys.argv[2]
    out_dir = sys.argv[3]

    var, exones = leer_excel(excel)
    genes_acmg, variantes_2do, exones_2do = busqueda_genes_acmg(acmg,var,exones)

    #guardamos
    crear_excel(variantes_2do,exones_2do,genes_acmg,out_dir)



if __name__ == '__main__':
    main(sys.argv)


