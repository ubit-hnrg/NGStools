#!/usr/bin/python
from __future__ import division

from statsmodels.stats.weightstats import DescrStatsW
import sys
import pandas as pd
import numpy as np
import os
import argparse

###############################################################################
#coverage_hist command line with bedtools
#bedtools coverage -a ./caputre_kit.bed -b sample.bam -hist > sample.bam.hist
## Header:
#echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i','--ensembl_coverage_histogram', help=' "bedtools coverage -hist" between library_kit and sample bam file')
parser.add_argument('-o','--output_local_report') #default='global_coverage_statistics.tsv'
parser.add_argument('-s','--samplename')

args =  parser.parse_args()

#ffile='EB802_exon_filtered_coverage.tsv.gz'
#ffile = sys.argv[1]
ffile = args.ensembl_coverage_histogram  # input file
sample = args.samplename
output_local_coverage = args.output_local_report
signif = 2 # harcoded
    
    
def depth_fraction(coverage,thr=0,ZeroDepth=False):
    if not ZeroDepth:
        condition = coverage['DP']>=thr
    else:
        condition = coverage['DP']== thr

    b = float(coverage.BPs.sum()) # total BPs

    if b!= 0:
        return coverage[condition].BPs.sum()/b
    else:
        return 0
    
    
def localdepth(coverage_hist):


    coverage_hist['cumsum'] = 1- coverage_hist.frequency.cumsum()
    weighted_stats = DescrStatsW(coverage_hist.DP-1, weights=coverage_hist.BPs, ddof=0)
    
    
   
    local_depth={}
    local_depth.update({'mean_DP':round(weighted_stats.mean,signif)})
    #local_depth.update({'median_DP':weighted_stats.quantile(0.5).values[0]})
    local_depth.update({'std_DP':round(weighted_stats.std,signif)})

    local_depth.update({'dp>=1':(round(depth_fraction(coverage_hist,thr=1),signif))*100})
    local_depth.update({'dp>=5':(round(depth_fraction(coverage_hist,thr=5),signif))*100})
    local_depth.update({'dp>=10':(round(depth_fraction(coverage_hist,thr=10),signif))*100})
    local_depth.update({'dp>=20':(round(depth_fraction(coverage_hist,thr=20),signif))*100})
    local_depth.update({'dp>=30':(round(depth_fraction(coverage_hist,thr=30),signif))*100})
    #local_depth.update({'mean_DP':round(weighted_stats.mean,signif)})

    #local_depth.update({'dp>=50':round(depth_fraction(coverage_hist,thr=50),signif)})
    #local_depth.update({'dp>=100':round(depth_fraction(coverage_hist,thr=100),signif)})
    return pd.Series(local_depth)

def main():
    ENS_coverage_hist = pd.read_csv(ffile,sep ='\t')
    results = ENS_coverage_hist.groupby(['transcriptID','exonNumber']).apply(localdepth)
    results = results.reset_index()
    info = ENS_coverage_hist.drop_duplicates(['transcriptID','exonNumber'])[['gene','exonNumber','transcriptID','IntervalLength']]
    results = pd.merge(results,info,on = ['transcriptID','exonNumber'])
    results = results[['transcriptID','gene','exonNumber','IntervalLength','mean_DP','std_DP','dp>=1','dp>=5','dp>=10','dp>=20','dp>=30']]
    results.sort_values(by=['gene','exonNumber'],inplace = True) ### sort para ordenar los exones en strand - ##agu 8/10
    results.to_csv(output_local_coverage,sep = '\t',index = False)

if __name__ == "__main__":
    main()
