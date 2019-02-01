#!/bin/bash

## This script take a fastq.gz file and generate the @RG (readgroup) id required for ubam generation
## Positional argument:
## fastq.gz file (mandatory)
## --check (optional flag for checking if the fastq file contains only one flowcell)

if [ ${2:-None} = '--check' ]
then
  numero_flowcells=$(zcat $1 |awk '{ if($1 ~/^@/) { print } }'|cut -f3 -d':' |sort|uniq |wc -l)
else
  numero_flowcells='1'
fi

if [ $numero_flowcells = '1' ]
then
  output=$(zcat $1 |head -n1| cut -f3,4 -d':' |sed -e 's/:/.Lane/g')
else
  output='more than one flowcell found'
fi
echo $output
