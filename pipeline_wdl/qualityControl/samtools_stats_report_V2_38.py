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

# Función para parsear los reportes de samtools
def parse_samtools_report_SN(samtools_stat_file):
    order = []
    sn = {}
    with open(samtools_stat_file, "r") as stats:
        for line in stats:
            if re.match("^SN", line):
                parts = line.split('\t')
                sn[parts[1].strip(':')] = float(parts[2].strip())
                order.append(parts[1].strip(':'))
    return pd.Series(sn).fillna(0.0)

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

# Construcción de reporte de porcentajes
percents_tso = pd.concat([
    samtools_kit_report[['reads properly paired', 'reads duplicated', 'reads MQ0']] / args.total_reads * 100,
    pd.Series({
        'maximum length': samtools_kit_report['maximum length'],
        'error rate': samtools_kit_report['error rate'] * 100,
        'bases mapped (cigar)': samtools_kit_report['bases mapped (cigar)'] / (10**9),
        'bases_before': args.total_bases_before / (10**9),
        'bases_after': args.total_bases_after / (10**9),
        'Gbases without dup mapped': bases_on_library_nodup / (10**9),
        'efficiency': efficiency,
        'filtering_ratio': filtering_ratio,
        'on_target_ratio': on_target_ratio,
        'unduplicated_ratio': duplicated_ratio,
        'bases inside the target': samtools_kit_report.get('bases inside the target', 0) / (10**9),
        'Number of reads properly paired': args.total_reads
    })
])

percents_tso.index = percents_tso.index + ' in Library'
percents_tso.rename(index={
    'Number of reads properly paired in Library': 'Number of reads properly paired',
    'reads properly paired in Library': 'Percent of reads properly paired in Library',
    'bases inside the target in Library': 'Target size (GBPs)',
    'bases mapped (cigar) in Library': 'Gbases mapped (cigar) in Library'
}, inplace=True)

percents_tso.index = percents_tso.index.str.replace(' ', '-')

# Renombre de índices para claridad
percents_tso.rename(index={
    'On-target[%]-in-Library': 'Bases-On-target[%]',
    'On-target_raw[%]-in-Library': 'Bases-On-target_raw[%]',
    'bases_before-in-Library': 'GBases_before_trimm',
    'bases_after-in-Library': 'GBases_after_trimm',
    'filtering_ratio-in-Library': 'filtering_ratio'
}, inplace=True)

# Selección de reporte final y formateo
report = percents_tso.loc[[...]] # Aquí deberías incluir los índices exactos que deseas incluir en tu reporte final

# Formateo y guardado del reporte
report.to_csv(args.output_file, sep='\t')