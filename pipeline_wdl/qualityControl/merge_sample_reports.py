#!/usr/bin/python
import sys
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--report_files', help='file containing one report file per line')
parser.add_argument('-o','--output_merged_report') 

args =  parser.parse_args()

reports = args.report_files
outfile = args.output_merged_report
merged_report = []

reports = pd.read_csv(reports,header=None).values.flatten()
for rep in reports:
    print rep
    sample = os.path.basename(rep).split('_')[0]
    print sample
    df = pd.read_table(rep,index_col=[0])
    df.columns=[sample]
    merged_report.append(df)

merged_report = pd.concat(merged_report,axis = 1)
merged_report.index.name=None
merged_report.to_csv(outfile,sep ='\t')

