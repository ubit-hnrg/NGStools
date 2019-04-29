# encoding=utf8
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import argparse
import pandas as pd
from xlrd import open_workbook
import csv
import os

parser = argparse.ArgumentParser(prog='parse-sheets.py',description='parse first column of all sheets in a xlsx file', usage='%(prog)s [options]')

parser.add_argument('-i','--input_xlsx', help='modified annovar output (with removed unammed fields)')
parser.add_argument('-o','--output_path', default='/home/hnrg/metadataHNRG/',help='path for writing outputs')
args =  parser.parse_args()

inputfile = args.input_xlsx
outpath = args.output_path
wb = open_workbook(inputfile)

for i in range(1, wb.nsheets):
       sheet = wb.sheet_by_index(i)
       sampleID=sheet.name.split('-')[0]

       #make dir
       folder=os.path.join(outpath,sampleID)
       outfile=os.path.join(folder,sampleID+'_genelist.txt')
       if not os.path.exists(folder):
              os.mkdir(folder)
       
       #get column values
       cells = sheet.col(0)[1:]
       values = [x.value for x in cells if x.value is not '']

#       l = len(values)-1;
#       k=0
       with open(outfile, 'w') as f:
              for item in values:
                     if item !='':
                            f.write("%s\t\n" % item)

