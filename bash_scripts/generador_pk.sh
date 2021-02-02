#!/bin/bash

set -e

std_1=$1
patron=$2
out=$3

for f in $std_1; do nombre=$(basename $std_1 $patron'.tsv'); awk -v nom=$nombre 'BEGIN{OFS="\t"} NR>1{$(NF+1)=nom} 1' $f > ./$nombre$out; done
