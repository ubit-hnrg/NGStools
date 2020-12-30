#!/usr/bin/python

import pandas as pd
import numpy as np
import os
import csv

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l','--lista_genes', help='archivo con genes, uno por linea')
parser.add_argument('-ic','--intervalo_de_captura', help='intervalo de captura acotado a ensembl, esta en la carpeta de resultados de cada experimento')
parser.add_argument('-o','--output', help='bed file con genes virtuales')
parser.add_argument('-gne','--out_no_encontrado', help='lista de genes no encontrados')
parser.add_argument('-ge','--out_encontrados', help='lista de genes encontrados')



args =  parser.parse_args()

lista_genes_path = args.lista_genes
int_cap_path = args.intervalo_de_captura
out = args.output
out_no_enc = args.out_no_encontrado
out_enc = args.out_encontrados


def lista_genes(path):
    '''
    lee archivo y genera una lista de genes
    '''

    with open(path) as f:
        genes = f.read().splitlines()
        genes_out = [item.strip() for item in genes if item]
        
    return genes_out

def no_encontrados(lista,intervalo_captura):
    '''
    devuelve los genes no encontrados en la lista de genes suministrada x el usuario
    '''
    aux_no = []
    aux_si = []
    for i in lista:
        if i not in intervalo_captura:
            aux_no.append(i)
        else: 
            aux_si.append(i)
    return aux_si, aux_no

def intervalo_captura(lista_genes,intervalo_experimento):

    experiment_lib = pd.read_csv(intervalo_experimento, sep = '\t',header = None)
    experiment_lib.rename(columns={0:'CHR',1:'START',2:'END',3:'ENSEMBL_ID',4:'GENE_NAME',5:'EXON_NUMBER',6:'STRAND'},inplace = True)


    panel_digital = experiment_lib[experiment_lib.GENE_NAME.isin(lista_genes)]
    genes_lista = panel_digital['GENE_NAME'].drop_duplicates().tolist()
    bases_panel_digital = sum(panel_digital.END - panel_digital.START)
    n_genes = panel_digital.GENE_NAME.drop_duplicates().count()
    n_exones = panel_digital.shape[0]
    
    return panel_digital, genes_lista, bases_panel_digital, n_genes, n_exones
    

def main():
    #if len(args) != 3:
    #    raise SystemExit('Uso adecuado: %s archivo_lista_genes archivo_intervalo_ensembl_2_lib' % args[0])

 
    genes_panel = lista_genes(lista_genes_path) ### 
    panel_digital, genes_lista, bases_panel_digital, n_genes, n_exones = intervalo_captura(genes_panel,int_cap_path)
    encontrado, no_found = no_encontrados(genes_panel,genes_lista)
    out_nofound = pd.DataFrame(no_found,columns=['genes_no_encontrados_en_lista(RECHECK)'])
    out_genes = pd.DataFrame(encontrado, columns=['genes en el panel digital'])

    #print(f'El panel digital {os.path.splitext(os.path.basename(lista_genes_path))[0]} tiene: {n_genes} genes, {n_exones} exones, {bases_panel_digital} bases.')
    #print(f'Los genes: {no_found} de la lista ingresada no se encuentran en el intervalo')
    
    panel_digital.to_csv(out,sep='\t', index = False, header = None)
    out_nofound.to_csv(out_no_enc,sep='\t', index = False)
    out_genes.to_csv(out_enc,sep='\t', index = False)


if __name__ == '__main__':
    #import sys
    #main(sys.argv)
    main()

                        