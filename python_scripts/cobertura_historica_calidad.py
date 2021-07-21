#!/usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np  
import sys, os, argparse
import statsmodels.stats.api as sms


def leer_estadistico(path, aux_var):

    '''
    Lee el archivo de estadistica generado en el pipeline de analisis de datos NGS 
    y agrega varias columnas. para qual_report
    '''

    with open(path, 'rt') as files:
        for linea in files:
            
            input_excel = linea.rstrip('\n')
            archi_calidad = pd.read_excel(input_excel, sheet_name = 'Alineamiento', index_col = 'Unnamed: 0').transpose()
            
            try:
                archi_calidad = archi_calidad.apply(lambda x: pd.to_numeric(x.str.strip(' %')))
            except: 
                archi_calidad = archi_calidad
            
            nombre_experimento = os.path.basename(linea).split('_')[0]
            
            archi_calidad['run'] = nombre_experimento
            aux_var.append(archi_estadistica)
        
        calidad = pd.concat(aux_var)
        calidad = calidad.reset_index()
        calidad = calidad.rename(columns={'index':'sample'})

    return calidad

def eficiencia(input_file):
    #cobertura_historica_variantes = input_file.groupby(['gene','transcriptID','exonNumber'])['start','end','strand','IntervalLength','dp>=1','dp>=10','dp>=20','dp>=30','dp>=50','dp>=100'].mean()
    historico_eficiencia = input_file.groupby(['run'])['efficienciy-in-Library'].mean()
    ci95_eficiencia = sms.DescrStatsW(historico_eficiencia).tconfint_mean(0.05)

    media_eficiencia = round(input_file.groupby(['run'])['efficienciy-in-Library'].mean().mean(),3)
    CI95 = str(round(ci95_eficiencia[0],3)) + '-' + str(round(ci95_eficiencia[1],3))

    min_efi = round(input_file.groupby(['run'])['efficienciy-in-Library'].min().mean(),3)
    max_efi = round(input_file.groupby(['run'])['efficienciy-in-Library'].max().mean(),3)
    rango = str(min_efi) + '-' + str(max_efi)

    filtering_ratio = round(input_file.groupby(['run'])['filtering_ratio'].mean().mean(),3)
    on_target_ratio = round(input_file.groupby(['run'])['on_target_ratio-in-Library'].mean().mean(),3)
    unduplicated_ratio = round(input_file.groupby(['run'])['unduplicated_ratio-in-Library'].mean().mean(),3)

    historico_gbases_raw = input_file.groupby(['run'])['GBases_before_trimm'].sum()
    gbases_raw = sms.DescrStatsW(historico_gbases_raw).tconfint_mean(0.05)
    CI95_gbases = str(round(gbases_raw[0],3)) + '-' + str(round(gbases_raw[1],3))
    media_gbases_raw = round(input_file.groupby(['run'])['GBases_before_trimm'].sum().mean(),3)

    eficiencia = {'Media':[media_eficiencia],
              'CI95':[CI95],
              'Rango':[rango],
              'CI95_Gbases_raw':[CI95_gbases],
              'mean_Gbases':[media_gbases_raw],
              'Filtering_ratio':[filtering_ratio],
              'On_target_ratio':[on_target_ratio],
              'Unduplicated_ratio':[unduplicated_ratio],
              'N_samples':[input_file.shape[0]]
             }

    eficiencia_out = pd.DataFrame(eficiencia, columns =['Media','CI95','Rango','Filtering_ratio','On_target_ratio','Unduplicated_ratio','mean_Gbases','CI95_Gbases_raw','N_samples'], index=['TSO', 'T1','AG20200820','T1_MET','T1_ENDO','T1_DERMATO','IN_HE'] ).rename_axis('Eficiencia')


    return eficiencia_out
    
def main(argv):
    aux_variantes = []

    if(len(sys.argv) > 1):
        if(not os.path.isfile(sys.argv[1])):
            print(sys.argv[1], "no se reconoce el archivo")
            sys.exit(1)
    inputs = sys.argv[1]
    out_dir = sys.argv[2]

   
    estadisticas_concatenadas = leer_estadistico(inputs, aux_variantes)
    
    res = eficiencia(estadisticas_concatenadas)
    res.to_csv(out_dir,sep='\t',index=True)
   

if __name__ == '__main__':
    main(sys.argv)

