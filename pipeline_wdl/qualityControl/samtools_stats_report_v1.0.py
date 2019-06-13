#!/usr/bin/python
import pandas as pd
import sys 
import argparse

parser = argparse.ArgumentParser(prog='samtools_stat_report.py',description='Make aligment report from samtools stat reports', usage='%(prog)s  --samtools_global_report --samtools_library_report --output_file')
parser.add_argument('-N','--total_reads', type=float,help='Number of total reads that passed quality criteria')
parser.add_argument('-l','--samtools_library_report', help='library kit restricted report from samtools stat tool')
parser.add_argument('-o','--output_file')

args = parser.parse_args()
total_reads = args.total_reads
samtools_kit_report_file = args.samtools_library_report
output = args.output_file



#sample=sys.argv[1]
#path = '/home/hnrg/resultsHNRG/%s'%sample
#samtools_global_report_file = path + '/%s_samtools.stats'%sample
#samtools_kit_report_file = path + '/%s_TSO_samtools.stats'%sample
#output = path + '/%s_samtools_report.tsv'%sample

def parse_samtools_report_SN(samtools_stat_file):
    order=[]
    sn = {}
    import re
    stats = open(samtools_stat_file, "r")

    for line in stats:
        if re.match("^SN", line):
            line = line.split('\t')
            d = {line[1].strip(':'):line[2].strip('\n')}
            order.append(line[1].strip(':'))
            sn.update(d)
    samtools_report = pd.Series(sn).astype(float).fillna(0.0)
    return samtools_report


#samtools_global_report = parse_samtools_report_SN(samtools_global_report_file)
samtools_kit_report = parse_samtools_report_SN(samtools_kit_report_file)


#absolute_reads_tso = samtools_kit_report['reads properly paired']
#absolute_reads_tso.rename(columns={'reads properly paired':'Number of reads properly paired'},inplace = True)

percents_tso = 100*samtools_kit_report[['reads properly paired','reads duplicated','reads MQ0']]/float(total_reads)
#percents_tso = samtools_kit_report[['reads properly paired','reads duplicated','reads MQ0']]

percents_tso = percents_tso.append(samtools_kit_report[['average quality','maximum length']])
percents_tso = percents_tso.append(samtools_kit_report[['error rate']]*100)

percents_tso.rename(index={'reads properly paired':'Percent of reads properly paired'},inplace = True)

#percents_tso = percents_tso.append(absolute_reads_tso)


percents_tso = percents_tso.append(pd.Series([total_reads],index=['Number of reads properly paired'])) # !!!! uncomment!!!!
percents_tso.index = percents_tso.index+' in Library'

percents_tso.index = percents_tso.index.str.replace(' ','-')

report = percents_tso.loc[[
#    u'raw-total-sequences',
             u'Number-of-reads-properly-paired-in-Library',
             u'Percent-of-reads-properly-paired-in-Library',
#             u'reads-duplicated',
             u'reads-duplicated-in-Library',
#             u'reads-MQ0',
             u'reads-MQ0-in-Library',
#             u'error-rate',
             u'error-rate-in-Library',
             u'maximum-length-in-Library']]

report.iloc[1:5] = report.iloc[1:5].round(2).map('{:,.2f} %'.format)
report[0] = '{:.0f}'.format(report[0])
report.to_csv(output,sep='\t')

#np.average(covTSO['cov'].values,weights= covTSO.cuentas.values)

