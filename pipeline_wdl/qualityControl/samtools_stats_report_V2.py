#!/usr/bin/python
import pandas as pd
import sys 
import argparse

parser = argparse.ArgumentParser(prog='samtools_stat_report.py',description='Make aligment report from samtools stat reports', usage='%(prog)s  --samtools_global_report --samtools_library_report --output_file')
parser.add_argument('-N','--total_reads', type=float,help='Number of total reads that passed quality criteria')
parser.add_argument('-ba','--total_bases_after', type=int,help='Number of total bases after fastq filtering')
parser.add_argument('-bb','--total_bases_before', type=int,help='Number of total bases before fastq filtering')
parser.add_argument('-l','--samtools_library_report', help='library kit restricted report from samtools stat tool')
parser.add_argument('-o','--output_file')

args = parser.parse_args()
total_reads = args.total_reads
bases_after = args.total_bases_after
bases_before = args.total_bases_before
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


samtools_kit_report = parse_samtools_report_SN(samtools_kit_report_file)
percents_tso = (100*samtools_kit_report[['reads properly paired','reads duplicated','reads MQ0']]/float(total_reads)).astype('float64')
percents_tso = percents_tso.append(samtools_kit_report[['maximum length']]).astype('int64')
percents_tso = percents_tso.append(samtools_kit_report[['error rate']]*100)


####bases
percents_tso = percents_tso.append(samtools_kit_report[['bases mapped (cigar)']])
#cigar_after = (100*samtools_kit_report[['bases mapped (cigar)'][0]]/float(bases_after))
#cigar_before = (100*samtools_kit_report[['bases mapped (cigar)'][0]]/float(bases_before))
cigar_after = (100*samtools_kit_report[['bases inside the target'][0]]/float(bases_after))
cigar_before = (100*samtools_kit_report[['bases inside the target'][0]]/float(bases_before))


#print(samtools_kit_report[['bases mapped (cigar)'][0]])
#print(samtools_kit_report[['bases mapped (cigar)']])
percents_tso = percents_tso.append(pd.Series([cigar_after],index=['On target[%]']))
percents_tso = percents_tso.append(pd.Series([cigar_before],index=['On target_raw[%]']))
percents_tso = percents_tso.append(pd.Series([bases_before],index=['bases_before']))
percents_tso = percents_tso.append(pd.Series([bases_after],index=['bases_after']))
#print(pd.Series([cigar_after],index=['CIGAR after']))
#print(pd.Series([cigar_before],index=['CIGAR before']))

percents_tso = percents_tso.append(samtools_kit_report[['bases inside the target']])
#### ese reads properly paired no deberia ser el de cigar? 
percents_tso = percents_tso.append(pd.Series([total_reads],index=['Number of reads properly paired']))

percents_tso.index = percents_tso.index+' in Library'
percents_tso.rename(index={'Number of reads properly paired in Library':'Number of reads properly paired'},inplace = True)
percents_tso.rename(index={'reads properly paired in Library':'Percent of reads properly paired in Library'},inplace = True)
percents_tso.rename(index={'bases inside the target in Library':'Bases inside the target'},inplace = True)

percents_tso.index = percents_tso.index.str.replace(' ','-')
print(percents_tso)

###renombro ontarge[%]
percents_tso.rename(index={'On-target[%]-in-Library':'On-target[%]'},inplace = True)
percents_tso.rename(index={'On-target_raw[%]-in-Library':'On-target_raw[%]'},inplace = True)
percents_tso.rename(index={'bases_before-in-Library':'Bases_before_trimm'},inplace = True)
percents_tso.rename(index={'bases_after-in-Library':'Bases_after_trimm'},inplace = True)

report = percents_tso.loc[[
#    u'raw-total-sequences',
             u'Number-of-reads-properly-paired',
             u'Percent-of-reads-properly-paired-in-Library',
#             u'reads-duplicated',
             u'reads-duplicated-in-Library',
#             u'reads-MQ0',
             u'reads-MQ0-in-Library',
#             u'error-rate',
             u'error-rate-in-Library',
             u'maximum-length-in-Library',
             u'On-target[%]',
             u'On-target_raw[%]',
             u'Bases_before_trimm',
             u'Bases_after_trimm',
             u'Bases-inside-the-target',
             u'bases-mapped-(cigar)-in-Library']]

report.iloc[1:5] = report.iloc[1:5].round(2).map('{:,.2f} %'.format)
#report.iloc[6:8] = report.iloc[7:8].round(2).map('{:,.2f} %'.format)
#report.iloc[9:10] = report.iloc[9:10].map('{:.2f}'.format)
report[0] = '{:.0f}'.format(report[0])
print(report)
report.to_csv(output,sep='\t')

#np.average(covTSO['cov'].values,weights= covTSO.cuentas.values)

