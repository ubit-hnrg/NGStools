#!/usr/bin/env python3
from __future__ import division
import pandas as pd
import argparse
import re

# Setting up the argument parser
parser = argparse.ArgumentParser(prog='samtools_stat_report.py', description='Make alignment report from samtools stat reports')
parser.add_argument('-N', '--total_reads', type=float, help='Number of total reads that passed quality criteria')
parser.add_argument('-ba', '--total_bases_after', type=int, help='Number of total bases after fastq filtering')
parser.add_argument('-bb', '--total_bases_before', type=int, help='Number of total bases before fastq filtering')
parser.add_argument('-l', '--samtools_library_report', help='Library kit restricted report from samtools stat tool')
parser.add_argument('-d', '--samtools_library_report_without_duplicates', help='Library kit restricted report from samtools stat tool without duplicates')
parser.add_argument('-o', '--output_file', help='Output file path')

args = parser.parse_args()

# Parsing the command line arguments
totalreads = args.total_reads
bases_after = args.total_bases_after
bases_before = args.total_bases_before
samtools_kit_report_file = args.samtools_library_report
samtools_kit_report_file_without_duplicates = args.samtools_library_report_without_duplicates
output_file = args.output_file

def parse_samtools_report_SN(samtools_stat_file):
    order = []
    sn = {}
    with open(samtools_stat_file, "r") as stats:
        for line in stats:
            if re.match("^SN", line):
                parts = line.split('\t')
                key = parts[1].strip(':')
                value = parts[2].strip()
                order.append(key)
                sn[key] = value
    samtools_report = pd.Series(sn).astype(float).fillna(0.0)
    return samtools_report

samtools_kit_report = parse_samtools_report_SN(samtools_kit_report_file)
samtools_kit_report_without_duplicates = parse_samtools_report_SN(samtools_kit_report_file_without_duplicates)

bases_on_target_no_dup = samtools_kit_report_without_duplicates['bases mapped (cigar)']
bases_on_library = samtools_kit_report['bases mapped (cigar)']
bases_on_library_nodup = bases_on_target_no_dup

efficiency = bases_on_library_nodup / bases_before
filtering_ratio = bases_after / bases_before
on_target_ratio = bases_on_library / bases_after
duplicated_ratio = bases_on_library_nodup / bases_on_library

percents_tso = (100 * samtools_kit_report[['reads properly paired', 'reads duplicated', 'reads MQ0']] / totalreads).append(
    samtools_kit_report[['maximum length', 'error rate']].multiply(100)).append(
    samtools_kit_report[['bases mapped (cigar)']].divide(10**9))

Gb_before = bases_before / 10**9
Gb_after = bases_after / 10**9
Gbases_nodup = round(bases_on_library_nodup / 10**9, 6)

percents_tso = percents_tso.append(pd.Series([Gb_before], index=['bases_before']))
percents_tso = percents_tso.append(pd.Series([Gb_after], index=['bases_after']))
percents_tso = percents_tso.append(pd.Series([Gbases_nodup], index=['Gbases without dup mapped']))
percents_tso = percents_tso.append(pd.Series([efficiency], index=['efficiency']))
percents_tso = percents_tso.append(pd.Series([filtering_ratio], index=['filtering_ratio']))
percents_tso = percents_tso.append(pd.Series([on_target_ratio], index=['on_target_ratio']))
percents_tso = percents_tso.append(pd.Series([duplicated_ratio], index=['unduplicated_ratio']))

percents_tso.index = percents_tso.index + ' in Library'
percents_tso.rename(index={'Number of reads properly paired in Library': 'Number of reads properly paired'}, inplace=True)
percents_tso.rename(index={'reads properly paired in Library': 'Percent of reads properly paired in Library'}, inplace=True)
percents_tso.rename(index={'bases inside the target in Library': 'Target size (GBPs)'}, inplace=True)
percents_tso.rename(index={'bases mapped (cigar) in Library': 'Gbases mapped (cigar) in Library'}, inplace=True)

percents_tso.index = percents_tso.index.str.replace(' ', '-')

report = percents_tso.loc[['Number-of-reads-properly-paired', 'Percent-of-reads-properly-paired-in-Library', 'reads-duplicated-in-Library', 'reads-MQ0-in-Library', 'error-rate-in-Library', 'maximum-length-in-Library', 'GBases_before_trimm', 'GBases_after_trimm', 'Gbases-mapped-(cigar)-in-Library', 'Gbases-without-dup-mapped-in-Library', 'efficiency-in-Library', 'filtering_ratio', 'on_target_ratio-in-Library', 'unduplicated_ratio-in-Library', 'Target-size-(GBPs)']]

# Formatting and output to file
report.loc['Number-of-reads-properly-paired'] = '{:.0f}'.format(report['Number-of-reads-properly-paired'])
report.to_csv(output_file, sep='\t', float_format='%.4f')


