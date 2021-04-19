#!/usr/bin/python

import argparse
import csv, os, sys, gzip, argparse
from operator import itemgetter



#parser = argparse.ArgumentParser(prog='barcode-freq.py',description='counts barcodes freq. Prints on stdout', usage='%(prog)s  --fasta')
#parser.add_argument('--fasta', help='fasta file in gz for count barcodes')
#args = parser.parse_args()

def leo_fasta(fasta_in):

    barcodes = {}
    with gzip.open(fasta_in) as fastq:
        for line in fastq:
                if not line.startswith(b'@'): continue
                bc = line.decode("utf-8").split(':')[-1].strip()
                if bc not in barcodes:
                        barcodes[bc] = 1
                else:
                        barcodes[bc]+=1
    return barcodes






def main(argv):
    if len(sys.argv) != 2:
        raise SystemExit(f'Uso adecuado: {sys.argv[0]} fasta_file (*.gz)')
    
    #informe_camion(sys.argv[1],sys.argv[2],sys.argv[3])
    aux = leo_fasta(sys.argv[1])
    total = sum(aux.values())
    for k, v in sorted(aux.items(), key=itemgetter(1)):
        print(k, v, round(v/total*100, 2))
    
if __name__ == '__main__':
    main(sys.argv)