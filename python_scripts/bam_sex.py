#!/usr/bin/python

from __future__ import division

import pandas as pd 
import numpy as np
import argparse
import csv
import subprocess
import os
#from plumbum import local

parser = argparse.ArgumentParser()#prog='bam_sex.py',description='get sample sexuality from bam file.', usage='%(prog)s  --bam --output_file')
parser.add_argument('-b','--bam', help='input bam file')
#parser.add_argument('-o','--output_file', help='a file with the sex of the sample')
args = parser.parse_args()

bam = args.bam
#out = args.output_file

#check_output
#x = subprocess.Popen(['samtools', 'view', bam ,'X','-c'],stdout=subprocess.PIPE)#.communicate()[0]#,shell = True).communicate()[0]#,stderr=subprocess.STDOUT)#.communicate()[0]#,stdout=subprocess.PIPE,stdin=subprocess.PIPE)#.communicate()[0]
x = subprocess.check_output(['samtools', 'view', bam ,'X','-c'])#,stdout=subprocess.PIPE)#.communicate()[0]#,shell = True).communicate()[0]#,stderr=subprocess.STDOUT)#.communicate()[0]#,stdout=subprocess.PIPE,stdin=subprocess.PIPE)#.communicate()[0]
y = subprocess.check_output(['samtools', 'view', bam ,'Y','-c'])


#y = subprocess.Popen(['samtools', 'view', bam ,'Y','-c'],stdout=subprocess.PIPE)#,stderr=subprocess.STDOUT)#,stderr=subprocess.STDOUT)#.communicate()[0]#shell = True).communicate()[0]#,stderr=subprocess.STDOUT)#.communicate()[0]#,stdout=subprocess.PIPE,stdin=subprocess.PIPE)#.communicate()[0]

#X = x.stdout.read()

#samtools = local['samtools']
#x = samtools['view',bam, 'X -c']

#y = samtools['view',bam, 'Y -c']
#print "holo"
#print x#.communicate()[0]
#print y#.communicate()[0]
sample = os.path.splitext(os.path.basename(bam))[0]#subprocess.call(["basename","bam",".bam"],shell = True)# ,stdout=subprocess.PIPE).communicate()[0]

#print sample
a = int(y)/int(x) >= 0.07
if a >= 0.07: 
    sex = 'M'
elif (a>= 0.03 or a < 0.07:
    sex = '?'
else:
    sex = 'F'

print ('{} \t {}'.format(sample,sex))


