#!/usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np  
import csv, sys, os, argparse


def leer_estadistico(path, aux_var):

    '''
    Lee el archivo de estadistica generado en el pipeline de analisis de datos NGS 
    y agrega varias columnas.
    '''

    with open(path, 'rt') as files:
        for linea in files:
            
            input_excel = linea.rstrip('\n')
            archi_estadistica = pd.read_excel(input_excel, sheet_name = 'ExonCoverage', index_col = 'gene')
            nombre_experimento = os.path.dirname(linea).split('/')[4]
            
            archi_estadistica['run'] = nombre_experimento
            aux_var.append(archi_estadistica)
        
        variantes = pd.concat(aux_var)
        variantes = variantes.reset_index()
        variantes = variantes.rename(columns={'index':'run'})

    return variantes

def agrupado_y_media(input_file):
    cobertura_historica_variantes = input_file.groupby(['gene','transcriptID','exonNumber'])['start','end','strand','IntervalLength','dp>=1','dp>=10','dp>=20','dp>=30','dp>=50','dp>=100'].mean()

    return cobertura_historica_variantes
    
def main(argv):
    aux_variantes = []

    if(len(sys.argv) > 1):
        if(not os.path.isfile(sys.argv[1])):
            print(sys.argv[1], "no se reconoce el archivo")
            sys.exit(1)
    inputs = sys.argv[1]
    out_dir = sys.argv[2]

   
    estadisticas_concatenadas = leer_estadistico(inputs, aux_variantes)
    
    res = agrupado_y_media(estadisticas_concatenadas)
    res.to_csv(out_dir,sep='\t',index=True)
   

if __name__ == '__main__':
    main(sys.argv)

