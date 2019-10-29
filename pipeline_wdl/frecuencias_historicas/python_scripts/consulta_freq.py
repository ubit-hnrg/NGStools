#!/usr/bin/python
import pandas as pd 
import numpy as np
import csv
import gzip
import allel

import sys 
import argparse


parser = argparse.ArgumentParser(prog='consulta_freq.py',description='gets historic frequency from chrom:pos file', usage='%(prog)s --file --coord --output_file')
parser.add_argument('-db','--freq_db', help='file freq_db vcf format')
parser.add_argument('-c','--coord_file', help='chromosome number:position query file')
parser.add_argument('-o','--output_file')

args = parser.parse_args()
db = args.freq_db
coord = args.coord_file
output = args.output_file

open_coord1 = pd.read_csv(coord,header=None ,sep="\t")
open_freq = pd.read_csv(db,sep="\t")
open_freq.CHROM = open_freq.CHROM.map(str)
open_freq.POS = open_freq.POS.map(int)

chrom_pos = open_coord1[0].str.split(':',expand=True)
chrom_pos.columns=['ch_num','pos']
chrom_pos['pos'] = pd.to_numeric(chrom_pos['pos'])
chrom_pos['ch_num'] = chrom_pos['ch_num'].map(str)

ind_bool=(open_freq['CHROM'].isin(chrom_pos["ch_num"]))& (open_freq['POS'].isin(chrom_pos["pos"]))
res = open_freq[ind_bool]   
res.to_csv(output,sep = '\t', index=False)

