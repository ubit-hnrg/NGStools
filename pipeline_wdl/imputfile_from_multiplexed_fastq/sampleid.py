import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fqfile', dest='fqfile',  help='fastqfiles separated by |')
parser.add_argument('--outfile', dest='outfile',  help='outfile |')
args = parser.parse_args()


with open(args.outfile, 'w') as outfile:

  with open(args.fqfile, 'r') as fq:
    for line in fq:
      R1=line.split('|')[0]
      sample=os.path.basename(R1).split('_')[0]
      res = '%s|%s|%s'%(sample,sample,line)
      outfile.write('%s'%res)
  fq.close()

outfile.close()
