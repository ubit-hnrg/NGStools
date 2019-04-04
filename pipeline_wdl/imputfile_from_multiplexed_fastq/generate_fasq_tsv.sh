#!/bin/bash
lsfile=$1
sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/|\n/g' $lsfile|sed '2,${s/|//g;n}'|sed -e ':a' -e 'N' -e '$!ba' -e 's/|\n/|/g' > fastqR1R2.txt
python sampleid.py --fqfile=./fastqR1R2.txt --outfile='samples.txt'
rm fastqR1R2.txt
