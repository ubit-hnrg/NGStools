#!/usr/bin/python
import sys
import pandas as pd
import os
import argparse
import numpy as np
import allel
import csv



parser = argparse.ArgumentParser()
parser.add_argument('-i','--freq_file', help='HNRG frequency file')
parser.add_argument('-o','--output_report') 

args =  parser.parse_args()
path = args.freq_file
outfile = args.output_report

## exploro numero maximo de alternativos
alternatives = pd.read_csv(path,comment='#',header =None,sep ='\t',usecols = [4])
maxalt = alternatives[4].str.split(',',expand = True).shape[1]
vcf = allel.read_vcf(path,alt_number=maxalt)

def analyze_alts(vcf, fila_vcf):
    alternativos = vcf['variants/ALT']
    
    a = alternativos[fila_vcf,:]
    a = a[a!='']
    alts = np.arange(len(a))+1
    
    genotipado = allel.GenotypeArray(vcf['calldata/GT'])
    samples_name = vcf['samples']
    res = []
    #freq_alts = []
   
    for alt in alts:
        A0 = genotipado[fila_vcf,:,0]
        A1 = genotipado[fila_vcf,:,1]
        S0 = samples_name[A0 == alt]
        S1 = samples_name[A1 == alt]
        samplelist = np.union1d(S0,S1) 
        res.append(samplelist)
        #alt_samples = ','.join(samplelist)
        #res.append(alt_samples) 
        
    return(res)

def pathology_freq(res):
    sampleSummary = []
    sampleDetail = []
    for i in range(len(res)):
        N1 = pd.Series(res[i])
        sdetail = ','.join(N1)
        sampleDetail.append(sdetail)                       
                                
        nombre = N1.str.extract('([a-zA-Z]+)([^a-zA-Z]+)', expand=True) ###match a group of letters: ([a-zA-Z]+) followed by a group of non letters: ([^a-zA-Z]+)
        n = nombre[0].value_counts().rename_axis('pathology_group').reset_index(name='counting')##nombre[0] xq solo me interesa el nombre de la patologia
        n['joined'] = '('+ n.counting.map(str) + ')' + n.pathology_group.map(str)
        alt_count = '-'.join(n.joined.tolist())
        sampleSummary.append(alt_count)
    return(sampleSummary,sampleDetail)

def refactor_alt_column_skitallel(variant_alt,joinpattern=','):
    return(pd.DataFrame(variant_alt).apply(lambda x: joinpattern.join(x[x!=u'']),1))

def join_sampleSummary_list(variant_alt,joinpattern='|'):
    return(pd.Series(variant_alt).apply(lambda x: str(joinpattern.join(x))))


detailArray = []
detail = []
summary = [] 
#for i in range(vcf['variants/ALT'].shape):
for i in np.arange(1000):
    detailArray.append(analyze_alts(vcf=vcf,fila_vcf=i))
    sSummary, sDetail = pathology_freq(detailArray[-1])
    summary.append(sSummary)
    detail.append(sDetail)

alt = refactor_alt_column_skitallel(vcf['variants/ALT'])
sampleSummary = join_sampleSummary_list(summary,joinpattern='|')
sampleDetail = join_sampleSummary_list(detail,joinpattern='|')

vcf2write = dict()
vcf2write.update({'CHROM':vcf['variants/CHROM']})
vcf2write.update({'POS':vcf['variants/POS']})
vcf2write.update({'REF':vcf['variants/REF']})
vcf2write.update({'ALT':alt})
vcf2write.update({'SamplesDetail':sampleDetail})
vcf2write.update({'SampleSummary':sampleSummary})

df = pd.DataFrame(vcf2write)[['CHROM','POS','REF','ALT','SampleSummary','SamplesDetail']]

df.to_csv(outfile,sep = '\t', index=False)