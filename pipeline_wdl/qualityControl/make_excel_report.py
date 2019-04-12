
#!/usr/bin/python
import pandas as pd
import os
import sys

report_pairs =sys.argv[1:-1]
output = sys.argv[-1]

#report_pairs = ['./result_fastp.tsv:Filtrado',\
#                './TSO20190328_samtools_report.tsv:Alineamiento',\
#                'TSO20190328_coverage_statistics.tsv:Profundidad']
#output = 'out.xlsx'

report_names  = [i.split(':')[1] for i in report_pairs]
report_fnames = [i.split(':')[0] for i in report_pairs]

with pd.ExcelWriter(output) as writer:
    for report_name, report_fname in zip(report_names,report_fnames):
        report_table = pd.read_table(report_fname)
        report_table.to_excel(writer, sheet_name=report_name,index = False)
        

