#!/bin/bash
numero_flowcells=$(zcat $1 |awk '{ if($1 ~/^@/) { print } }'|cut -f3 -d':' |sort|uniq |wc -l)
#echo $numero_flowcells

if [ $numero_flowcells = '1' ]
then
  output=$(zcat $1 |head -n1| cut -f3,4 -d':' |sed -e 's/:/.Lane/g')
else
  output='more than one flowcell found'
fi
echo $output
