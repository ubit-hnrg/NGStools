import pandas as pd
import glob
import os
import subprocess
import argparse


def get_args():
    parser = argparse.ArgumentParser(description=  'This script takes a bounch of vcfs and compute the observed frequecy for each observed variant')
    parser.add_argument('-i','--inputpath',required=True,help='input dir containig uncompresed vcf files')
    parser.add_argument('-o','--outpath',required=True)    
    args = parser.parse_args()
    return args



def quit_header_and_read_vcf(f,usecols = ['#CHROM','POS','REF','ALT','QUAL']):
    with open(f,'rb') as h:
        counter = 0
        for line in h:
            counter = counter +1
            if line.startswith('#CHROM'):
                break
    h.close()

    df = pd.read_csv(f,sep = '\t',skiprows=counter-1,usecols=usecols)
    df.rename(columns = {'#CHROM':'CHROM'},inplace = True)
    return df

def get_files(inputpath,pattern = '*0056*final*.vcf'):
    files = glob.glob(inputpath+pattern)
    return files


def store_vcfs(vcf_files):
    dflist = []
    for vcf in vcf_files:
        df = quit_header_and_read_vcf(vcf)
        df.drop_duplicates(subset=['CHROM','POS','REF','ALT'],inplace = True)
        dflist.append(df)
    return pd.concat(dflist)

def apply_filters(concatenated_vcfs,Q=1):
    ind  = concatenated_vcfs.QUAL>=Q
    return concatenated_vcfs[ind]
    

def variant_count(concatenated_vcfs,nsamples):
    grouped = concatenated_vcfs.groupby(['CHROM','POS','REF','ALT'])['QUAL'].count()
    grouped = grouped.reset_index().rename(columns={'QUAL':'vcount'})
    
    fmax = nsamples
    grouped['freq'] = grouped.vcount/float(fmax)
    grouped['fgroup'] = pd.cut(grouped.vcount,[0,fmax/float(20),fmax/float(10),fmax/float(2),fmax],labels = ['(2-5)%','[5-10)%','[10-50)%','[50-100]%'])
    return grouped



if __name__ ==  '__main__':
    args = get_args()
    inputpath, opath = args.inputpath, args.outpath
    vcf_files = get_files(inputpath)
    nsamples = len(vcf_files)
    concatenated_vcfs = store_vcfs(vcf_files)
    filtered_vcfs = apply_filters(concatenated_vcfs)
    vcount = variant_count(filtered_vcfs,nsamples)

    ## save output
    if not os.path.exists(opath):
        os.mkdir(opath)
    vcount.to_csv(opath+'local_frequency.csv',sep ='|')
