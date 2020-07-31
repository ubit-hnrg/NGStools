#!/usr/bin/python
from __future__ import division

from statsmodels.stats.weightstats import DescrStatsW
import sys
import pandas as pd
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i','--exon_coverage_regions', help='exon coverage regions in exons')
parser.add_argument('-g','--output_global_report') #default='global_coverage_statistics.tsv'
parser.add_argument('-e','--output_exon_report')# default='./coverage_statistics_by_exon.tsv')
parser.add_argument('-s','--samplename')

args =  parser.parse_args()



#ffile='EB802_exon_filtered_coverage.tsv.gz'
#ffile = sys.argv[1]
ffile = args.exon_coverage_regions  # input file
sample = args.samplename
output_global_coverage = args.output_global_report
output_coverage_by_exon = args.output_exon_report

coverage = pd.read_table(ffile)

signif=2



####aca creas una funcion para hacer la cuenta para eliminar la covertura 0 del experimento no? 
def depth_fraction(coverage,thr=0,ZeroDepth=False):
    if not ZeroDepth:
        condition = coverage['count']>=thr
    else:
        condition = coverage['count']== thr

#    return coverage[condition].count_length.sum()/float(coverage.count_length.sum())
    a = coverage[condition].count_length.sum()
    b = float(coverage.count_length.sum())
    #return coverage[condition].count_length.sum()/float(coverage.count_length.sum()) 
    return a/b if b != 0 else 0


# ensure intervals fall inside library kit
coverage = coverage[(coverage.chr==coverage.exon_chr)&(coverage.start>=coverage.exon_start)&(coverage.end<=coverage.exon_end)]
#sort and remove duplicates
coverage.sort_values(by=['chr','start','end'],ascending=[True,True,False])

coverage.drop_duplicates(inplace = True,subset=['chr','start','geneSymbol','ENSEMBL_ID'])

coverage['count_length']=coverage.end-coverage.start
coverage['exon_length']=coverage.exon_end-coverage.exon_start


# ## Compute global statistics

bases_without_reads = depth_fraction(coverage,ZeroDepth=True)
bases_greater1 = depth_fraction(coverage,thr=1)
bases_greater10 = depth_fraction(coverage,thr=10)
bases_greater20 = depth_fraction(coverage,thr=20)
bases_greater30 = depth_fraction(coverage,thr=30)
# uso estadisticos pesados popr la longitud de cada intervalo de profundidad constante
weighted_stats = DescrStatsW(coverage['count'], weights=coverage.count_length, ddof=0)


global_depth={}
global_depth.update({'mean_DP':round(weighted_stats.mean,signif)})
global_depth.update({'median_DP':weighted_stats.quantile(0.5).values[0]})
global_depth.update({'std_DP':round(weighted_stats.std,signif)})
global_depth.update({'q25_DP':weighted_stats.quantile(0.25).values[0]})
global_depth.update({'q75_DP':weighted_stats.quantile(0.75).values[0]})
global_depth.update({'q95_DP':weighted_stats.quantile(0.95).values[0]})

global_depth.update({'dp>=1':bases_greater1})
global_depth.update({'dp>=10':bases_greater10})
global_depth.update({'dp>=20':bases_greater20})
global_depth.update({'dp>=30':bases_greater30})
global_depth.update({'q95_DP':weighted_stats.quantile(0.95).values[0]})
#global_depth.update({'DP=0':round(bases_without_reads,signif)})

res = pd.Series(global_depth).to_frame()
res.columns=[sample]
global_statistics = res.loc[[u'dp>=1',u'dp>=10', u'dp>=20', u'dp>=30', 
         u'median_DP', u'q25_DP', u'q75_DP', u'q95_DP',
         u'mean_DP', u'std_DP'],:]

# # compute per gene statistics

depth1 = coverage.groupby(['geneSymbol','ENSEMBL_ID','exon_number','strand'])['count','count_length'].apply(lambda x: depth_fraction(x,thr=1))
depth10 = coverage.groupby(['geneSymbol','ENSEMBL_ID','exon_number','strand'])['count','count_length'].apply(lambda x: depth_fraction(x,thr=10))
depth20 = coverage.groupby(['geneSymbol','ENSEMBL_ID','exon_number','strand'])['count','count_length'].apply(lambda x: depth_fraction(x,thr=20))
depth30 = coverage.groupby(['geneSymbol','ENSEMBL_ID','exon_number','strand'])['count','count_length'].apply(lambda x: depth_fraction(x,thr=30))



depth_by_exon = pd.concat([depth1,depth10,depth20,depth30],axis=1)
depth_by_exon.columns=['dp1','dp10','dp20','dp30']
depth_by_exon =depth_by_exon.round(2)

exon_description = coverage[['chr','exon_start','exon_end','geneSymbol','ENSEMBL_ID','exon_number','strand']].drop_duplicates()

depth_by_exon = pd.merge(depth_by_exon,exon_description,how='left',left_index=True,right_on=['geneSymbol','ENSEMBL_ID','exon_number','strand'])
depth_by_exon = depth_by_exon[[u'geneSymbol',u'ENSEMBL_ID', u'exon_number',u'strand', u'dp1', u'dp10', u'dp20', u'dp30', u'chr', u'exon_start', u'exon_end']]
depth_by_exon.sort_values(by=['chr','exon_start','exon_end'],inplace = True)
#depth_by_exon.reset_index(inplace=True)


global_statistics.round(2).to_csv(output_global_coverage,sep='\t')
depth_by_exon.to_csv(output_coverage_by_exon,index = False,sep='\t')
