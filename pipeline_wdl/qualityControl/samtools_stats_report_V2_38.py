#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
import re

# Definición del analizador de argumentos
parser = argparse.ArgumentParser(prog='samtools_stat_report.py', description='Make alignment report from samtools stat reports')
parser.add_argument('-N', '--total_reads', type=float, help='Number of total reads that passed quality criteria')
parser.add_argument('-ba', '--total_bases_after', type=int, help='Number of total bases after fastq filtering')
parser.add_argument('-bb', '--total_bases_before', type=int, help='Number of total bases before fastq filtering')
parser.add_argument('-l', '--samtools_library_report', help='library kit restricted report from samtools stat tool')
parser.add_argument('-d', '--samtools_library_report_without_duplicates', help='library kit restricted report from samtools stat tool without duplicates')
parser.add_argument('-o', '--output_file')
args = parser.parse_args()


args = parser.parse_args()
totalreads = args.total_reads
bases_after = args.total_bases_after
bases_before = args.total_bases_before
samtools_kit_report_file = args.samtools_library_report
samtools_kit_report_file_without_duplicates = args.samtools_library_report_without_duplicates # ARI 19/09
output = args.output_file

# Función para parsear los reportes de samtools
def parse_samtools_report_SN(samtools_stat_file):
    sn = {}
    with open(samtools_stat_file, "r") as stats:
        for line in stats:
            if line.startswith("SN"):
                parts = line.split('\t')
                key = parts[1].strip(':')
                value = parts[2].strip()
                sn[key] = value
    samtools_report = pd.Series(sn).astype(float).fillna(0.0)
    return samtools_report

# Procesamiento de los reportes
samtools_kit_report = parse_samtools_report_SN(args.samtools_library_report)
samtools_kit_report_without_duplicates = parse_samtools_report_SN(args.samtools_library_report_without_duplicates)
bases_on_target_no_dup = samtools_kit_report_without_duplicates['bases mapped (cigar)']
bases_on_library = samtools_kit_report['bases mapped (cigar)']
bases_on_library_nodup = bases_on_target_no_dup

# Cálculo de métricas
efficiency = bases_on_library_nodup / args.total_bases_before
filtering_ratio = args.total_bases_after / args.total_bases_before
on_target_ratio = bases_on_library / args.total_bases_after
duplicated_ratio = bases_on_library_nodup / bases_on_library

# Preparación de la serie percents_tso con pd.concat
percents_tso_parts = [
    100 * samtools_kit_report[['reads properly paired', 'reads duplicated', 'reads MQ0']] / totalreads,
    samtools_kit_report[['maximum length']],
    samtools_kit_report[['error rate']] * 100,
    samtools_kit_report[['bases mapped (cigar)']] / 1e9,
    pd.Series({
        'bases_before': bases_before / 1e9,
        'bases_after': bases_after / 1e9,
        'Gbases without dup mapped': round(bases_on_library_nodup / 1e9, 6),
        'efficienciy': efficienciy,
        'filtering_ratio': filtering_ratio,
        'on_target_ratio': on_target_ratio,
        'unduplicated_ratio': duplicated_ratio,
        'Number of reads properly paired': totalreads
    })
]

percents_tso = pd.concat(percents_tso_parts)
percents_tso.index += ' in Library'

# Cambio de nombres de índices
index_renames = {
    'Number of reads properly paired in Library': 'Number of reads properly paired',
    'reads properly paired in Library': 'Percent of reads properly paired in Library',
    'bases inside the target in Library': 'Target size (GBPs)',
    'bases mapped (cigar) in Library': 'Gbases mapped (cigar) in Library',
    'On-target[%]-in-Library': 'Bases-On-target[%]',
    'On-target_raw[%]-in-Library': 'Bases-On-target_raw[%]',
    'bases_before-in-Library': 'GBases_before_trimm',
    'bases_after-in-Library': 'GBases_after_trimm',
    'filtering_ratio-in-Library': 'filtering_ratio'
}
percents_tso.index = percents_tso.index.str.replace(' ', '-')
percents_tso.rename(index=index_renames, inplace=True)

report = percents_tso.apply(lambda x: '{:.4f}'.format(x) if isinstance(x, float) else x)

report.to_csv(output, sep='\t')
