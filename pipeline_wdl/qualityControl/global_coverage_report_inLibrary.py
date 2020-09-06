#!/usr/bin/python
from __future__ import division

from statsmodels.stats.weightstats import DescrStatsW
import sys
import pandas as pd
import numpy as np
import os
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import LogLocator, AutoLocator, AutoMinorLocator
#

###############################################################################
#coverage_hist command line with bedtools
#bedtools coverage -a ./caputre_kit.bed -b sample.bam -hist > sample.bam.hist
## Header:
#echo -e 'chr\tstart\tend\tgene\tDP\tBPs\tIntervalLength\tfrequency' > header.txt
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-i','--library_coverage_histogram', help=' "bedtools coverage -hist" between library_kit and sample bam file')
parser.add_argument('-o','--output_global_report') #default='global_coverage_statistics.tsv'
parser.add_argument('-op','--output_plot')
parser.add_argument('-s','--samplename')

args =  parser.parse_args()

#ffile='EB802_exon_filtered_coverage.tsv.gz'
#ffile = sys.argv[1]
ffile = args.library_coverage_histogram  # input file
sample = args.samplename
output_global_coverage = args.output_global_report
out_plot = args.output_plot
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
    globalreport = globalreport[['mean_DP', 'median_DP', 'std_DP', 'q25_DP', 'q75_DP', 'q95_DP', 'dp>=1',
       'dp>=10', 'dp>=20', 'dp>=30', 'dp>=50', 'dp>=100']].to_frame()
    globalreport.columns=[sample]
    globalreport.to_csv(output_global_coverage,header = True,sep='\t')

    f = plt.figure(figsize=(14,6))
    ax1 = f.add_subplot(121)
    ax2 = f.add_subplot(122)
    ax1.plot(coverage_hist.DP, coverage_hist.frequency)
    ax1.set_yscale('log')
    ax1.set_xscale('symlog')
    ax1.set_xlabel('DP')
    ax1.set_ylabel('density')
    ax1.xaxis.set_minor_locator(LogLocator(subs=np.arange(2, 10)))
    ax1.grid(True, which="both", ls="--")

    ax2.plot(coverage_hist.DP, coverage_hist['cumsum'].values,'.-')
    ax2.set_ylim((-0.02,1))
    ax2.set_xscale('symlog')
    ax2.set_ylabel('Prob( bp > DP )')
    ax2.set_yticks(np.arange(0, 1., 0.1))
    ax2.set_xlabel('DP')
    ax2.xaxis.set_minor_locator(LogLocator(subs=np.arange(2, 10)))
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.grid(True, which="both", ls="--")
    f.savefig(out_plot, format='eps')

if __name__ == "__main__":
    main()

