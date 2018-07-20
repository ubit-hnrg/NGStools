#! /bin/bash
#!/bin/bash

while getopts ":i:b:o:s:" opt; do
  case $opt in
    i) ivcf="$OPTARG"
    ;;
    b) ibam="$OPTARG"
    ;;
    o) opath="$OPTARG"
    ;;
    s) samplename="$OPTARG"
    ;;
    g) gap="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&5
    ;;
  esac
done

printf "Argument ivcf is %s\n" "$ivcf"
printf "Argument ibam is %s\n" "$ibam"
printf "Argument opath is %s\n" "$opath"
printf "Argument samplename is %s\n" "$samplename"
printf "Argument delta is %s\n" "$gap"

#grep -v "^#" $ivcf | awk '{if ($2 > $delta) print $1, $2-$delta-1, $2+$delta-1;  else print $1, 0, $2+$delta-1; }' > $opath/$ivcf.bed
#samtools view -bh -L "${ivcf}" "${ibam}" > "${opath}"/"${samplename}"_reduced.bam
#samtools index "${opath}"/"${samplename}"_reduced.bam
