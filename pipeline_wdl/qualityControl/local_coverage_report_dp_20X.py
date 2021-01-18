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
parser.add_argument('-i','--digital_panel_coverage_histogram', help=' "bedtools coverage -hist" between digital_panel and sample bam file')
parser.add_argument('-o','--output_coverage_digital_report') #default='global_coverage_statistics.tsv'
parser.add_argument('-pn','--panel_name')

args =  parser.parse_args()


ffile = args.digital_panel_coverage_histogram  # input file
panel_name = args.panel_name
output_local_coverage = args.output_coverage_digital_report
signif = 2 # harcoded
    
########################################################################################################### ARI
def depth_fraction(coverage,thr=0,ZeroDepth=False): ##modifico, retorno la cantidad de bases totales.       #
    if not ZeroDepth:                                                                                       #
        condition = coverage['DP']>=thr                                                                     #
    else:
        condition = coverage['DP']== thr

    b = float(coverage.BPs.sum()) # total BPs
    if b!= 0:
        return b, coverage[condition].BPs.sum(), coverage[condition].BPs.sum()/b 
    else:
        return 0

def localdepth(coverage_hist):
    local_depth={}
    b, bases_20x, depth_20X = depth_fraction(coverage_hist,thr=20)
    local_depth.update({'bases_totales':int(b)})
    local_depth.update({'bases_20X':int(bases_20x)})

    local_depth.update({'dp>=20':round(depth_20X,signif)})
    #local_depth.update({'dp>=20':round(depth_fraction(coverage_hist,thr=30),signif)})
    return pd.Series(local_depth)                                                                  #
####################################################################################################

def main():
    cov_panel_digital = pd.read_csv(ffile,sep ='\t')
    results_20X = cov_panel_digital.groupby(['transcriptID','exonNumber']).apply(localdepth) ###sumo cada 
    res_20X = results_20X[results_20X['dp>=20']>0].copy()
    res_20X = res_20X.reset_index()
    
    ######
    genes_panel = cov_panel_digital.gene.drop_duplicates().count() 
    exones_panel = cov_panel_digital.drop_duplicates(['transcriptID','exonNumber']).shape[0]
    bases_panel = cov_panel_digital.BPs.sum()
    bases_20X = res_20X.bases_20X.sum()
    porcentaje_20x = 100*(bases_20X/bases_panel)
    prof_media = cov_panel_digital.DP.mean()

    ###
    panel_dig = {'N_genes':[genes_panel],'N_exones':[exones_panel],'N_bases_panel':[bases_panel],'bases_20X':int(bases_20X),'bases>=20X(%)':round(porcentaje_20x,2),'prof_media':round(prof_media,2)}
    panel_out = pd.DataFrame(panel_dig, columns = ['N_genes','N_exones','N_bases_panel','bases_20X','bases>=20X(%)','prof_media'],index = [panel_name])


    panel_out.to_csv(output_local_coverage,sep = '\t')#,index = False)

if __name__ == "__main__":
    main()
