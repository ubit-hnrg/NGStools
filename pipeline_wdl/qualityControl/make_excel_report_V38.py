#!/usr/bin/env python3
import pandas as pd
import sys
import xlsxwriter

# Obtener argumentos de la línea de comandos, excluyendo el nombre del script.
report_pairs = sys.argv[1:-1]
output = sys.argv[-1]

# Crear un workbook con la opción de 'constant_memory' para manejar grandes archivos.
workbook = xlsxwriter.Workbook(output, {'constant_memory': True})

# Separar los nombres de los reportes y los nombres de los archivos de los pares de argumentos.
report_names = [i.split(':')[1] for i in report_pairs]
report_fnames = [i.split(':')[0] for i in report_pairs]

for report_name, report_fname in zip(report_names, report_fnames):
    print(report_name)
    worksheet = workbook.add_worksheet(report_name)
    # Usar read_table con la opción 'sep' como '\t' para especificar que es un archivo separado por tabulaciones.
    report_table = pd.read_table(report_fname, sep='\t', encoding='utf-8')
    
    # Obtener dimensiones del DataFrame.
    row_max, col_max = report_table.shape
    columns = report_table.columns

    # Escribir el encabezado del archivo.
    for col_num, column_title in enumerate(columns):
        worksheet.write(0, col_num, column_title)

    # Escribir los datos, omitiendo los valores nulos.
    for row_num in range(row_max):
        for col_num in range(col_max):
            cell_value = report_table.iloc[row_num, col_num]
            if pd.notnull(cell_value):
                worksheet.write(row_num + 1, col_num, cell_value)

# Cerrar el workbook después de escribir todos los datos.
workbook.close()
