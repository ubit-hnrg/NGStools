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
parser.add_argument('-i','--library_coverage_histogram', help=' "bedtools coverage -hist" between library_kit and sample bam file')
parser.add_argument('-o','--output_global_report') #default='global_coverage_statistics.tsv'
parser.add_argument('-s','--samplename')

args =  parser.parse_args()

#ffile='EB802_exon_filtered_coverage.tsv.gz'
#ffile = sys.argv[1]
ffile = args.library_coverage_histogram  # input file
sample = args.samplename
output_global_coverage = args.output_global_report
signif = 2 # harcoded

def globaldepth(coverage_hist):

    coverage_hist['cumsum'] = 1- coverage_hist.frequency.cumsum()
    weighted_stats = DescrStatsW(coverage_hist.DP-1, weights=coverage_hist.BPs, ddof=0)
    
    global_depth={}
    global_depth.update({'mean_DP':round(weighted_stats.mean,signif)})
    global_depth.update({'median_DP':weighted_stats.quantile(0.5).values[0]})
    global_depth.update({'std_DP':round(weighted_stats.std,signif)})
    global_depth.update({'q25_DP':weighted_stats.quantile(0.25).values[0]})
    global_depth.update({'q75_DP':weighted_stats.quantile(0.75).values[0]})
    global_depth.update({'q95_DP':weighted_stats.quantile(0.95).values[0]})
    global_depth.update({'q95_DP':weighted_stats.quantile(0.95).values[0]})

    global_depth.update({'dp>=1':round(depth_fraction(coverage_hist,thr=1),signif)})
    global_depth.update({'dp>=10':round(depth_fraction(coverage_hist,thr=10),signif)})
    global_depth.update({'dp>=20':round(depth_fraction(coverage_hist,thr=20),signif)})
    global_depth.update({'dp>=30':round(depth_fraction(coverage_hist,thr=30),signif)})
    global_depth.update({'dp>=50':round(depth_fraction(coverage_hist,thr=50),signif)})
    global_depth.update({'dp>=100':round(depth_fraction(coverage_hist,thr=100),signif)})
    return(global_depth)
    
    
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
    
def main():
    coverage_hist = pd.read_csv(ffile,sep ='\t')
    globalreport = pd.Series(globaldepth(coverage_hist))
    globalreport.to_csv(output_global_coverage)


if __name__ == "__main__":
    main()

