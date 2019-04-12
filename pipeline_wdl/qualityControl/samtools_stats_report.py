#!/usr/bin/python
import pandas as pd
import sys 

sample=sys.argv[1]
path = '/home/hnrg/resultsHNRG/%s'%sample
samtools_global_report_file = path + '/%s.samtools_stats'%(sample,%sample)
samtools_kit_report_file = path + '/%s_TSO.samtools_stats'%(sample,%sample)
output = path + '/%s_samtools_report.tsv'%sample

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


samtools_global_report = parse_samtools_report_SN(samtools_global_report_file)
samtools_kit_report = parse_samtools_report_SN(samtools_kit_report_file)


percents = 100*samtools_global_report[['reads properly paired','reads duplicated','reads MQ0']]/samtools_global_report['raw total sequences']
percents = percents.append(samtools_global_report[['average quality','maximum length']])
percents = percents.append(samtools_global_report[['error rate']]*100)

percents_tso = 100*samtools_kit_report[['reads properly paired','reads duplicated','reads MQ0']]/samtools_global_report['raw total sequences']
#percents_tso = percents_tso.round(2)
percents_tso = percents_tso.append(samtools_kit_report[['average quality']])
percents_tso = percents_tso.append(samtools_kit_report[['error rate']]*100)

percents_tso.index = percents_tso.index+' in Library'
percents_tso.rename(columns={'reads properly paired':'reads properly paired'})

#merge both analysis
percents = percents.append(percents_tso)
percents = percents.append(samtools_global_report[['raw total sequences']])
percents.index = percents.index.str.replace(' ','-')


report = percents.loc[[u'raw-total-sequences',
             u'reads-properly-paired',
             u'reads-properly-paired-in-Library',
             u'reads-duplicated',
             u'reads-duplicated-in-Library',
             u'reads-MQ0',
             u'reads-MQ0-in-Library',
             u'error-rate',
             u'error-rate-in-Library',
             u'maximum-length']]

report.iloc[1:9] = report.iloc[1:9].round(2).map('{:,.2f} %'.format)
report[0] = '{:.0f}'.format(report[0])
report.to_csv(output,sep='\t')

#np.average(covTSO['cov'].values,weights= covTSO.cuentas.values)

