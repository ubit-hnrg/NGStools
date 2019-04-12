#!/usr/bin/python
import sys
import pandas as pd
import os

reports = sys.argv[1]
outfile = sys.argv[2]
merged_report = []

reports = pd.read_csv(reports).to_list()
for rep in reports:
    print rep
    sample = os.path.basename(rep).split('_')[0]
    print sample
    df = pd.read_table(rep,index_col=[0],header =None)
    df.columns=[sample]
    merged_report.append(df)

merged_report = pd.concat(merged_report,axis = 1)
merged_report.index.name=None
merged_report.to_csv(outfile,sep ='\t')

