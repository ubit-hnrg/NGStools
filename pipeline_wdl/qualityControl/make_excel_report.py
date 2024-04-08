#!/usr/bin/python3
import pandas as pd
import os
import sys
import openpyxl
import xlsxwriter

report_pairs =sys.argv[1:-1]
output = sys.argv[-1]
workbook = xlsxwriter.Workbook(output, {'constant_memory': True})  # ROBUSTO PARA GRANDES FILES.


report_names  = [i.split(':')[1] for i in report_pairs]
report_fnames = [i.split(':')[0] for i in report_pairs]


for report_name, report_fname in zip(report_names,report_fnames):
    print report_name
    worksheet = workbook.add_worksheet(report_name)
    report_table = pd.read_table(report_fname,encoding='utf-8')
    
    shape = report_table.shape
    row_max = shape[0]
    col_max = shape[1]
    columns = report_table.columns
    index = report_table.index
    nulls = report_table.isnull()  #para evitar errores escribiendo nulls
    #write sheet
    row = 0
    for col in range(0, col_max):
            cell = columns[col]
            worksheet.write(row, col, cell)

    # write the remaining datas
    for row in range(0, row_max):
        for col in range(0, col_max):
            cell = report_table.iloc[row,col]
            null = nulls.iloc[row,col]
            if not null:
                worksheet.write(row+1, col, cell)
workbook.close()

#with pd.ExcelWriter(output) as writer:
#    for report_name, report_fname in zip(report_names,report_fnames):
#        report_table = pd.read_table(report_fname)
#        report_table.to_excel(writer, sheet_name=report_name,index = False)
        

